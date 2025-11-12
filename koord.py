# koord.py — УНИВЕРСАЛЬНАЯ ВЕРСИЯ (работает на любой osmnx)

import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point, box
import folium
from folium.plugins import HeatMap
import osmnx as ox
import requests
import time
from shapely.geometry import Polygon

class StoreLocationFinder:
    def __init__(self, stores_data):
        self.stores_data = stores_data
        self.geometry = [Point(xy) for xy in zip(stores_data.longitude, stores_data.latitude)]
        self.gdf_stores = gpd.GeoDataFrame(stores_data, geometry=self.geometry, crs="EPSG:4326")

    def create_search_grid(self, bounds, grid_size=0.0025):
        min_lon, min_lat, max_lon, max_lat = bounds
        num_lon = int(np.ceil((max_lon - min_lon) / grid_size)) + 1
        num_lat = int(np.ceil((max_lat - min_lat) / grid_size)) + 1

        lons = np.linspace(min_lon, max_lon, num_lon)
        lats = np.linspace(min_lat, max_lat, num_lat)

        grid_points = [Point(lon, lat) for lon in lons for lat in lats]
        grid = gpd.GeoDataFrame(geometry=grid_points, crs="EPSG:4326")

        bounds_poly = box(min_lon, min_lat, max_lon, max_lat)
        mask = grid.geometry.within(bounds_poly)
        grid = grid[mask].reset_index(drop=True)

        print(f"Сетка: {len(grid)} точек (внутри границ)")
        return grid

    def find_best_location(self, bounds, grid_size=0.0025, water_gdf=None):
        print("Поиск лучшего места...")
        grid = self.create_search_grid(bounds, grid_size)

        # === ФИЛЬТРАЦИЯ ПО ВОДЕ ===
        if water_gdf is not None and not water_gdf.empty:
            print("Исключаем воду...")
            try:
                water_union = water_gdf.unary_union
                before = len(grid)
                grid = grid[~grid.geometry.within(water_union)]
                print(f"После исключения воды: {len(grid)} из {before}")
            except:
                print("Не удалось объединить водоёмы — пропускаем фильтр")

        if len(grid) == 0:
            print("Все точки на воде! Попробуй уменьшить grid_size.")
            return None, None, None

        def min_distance(p):
            return min(p.distance(s) for s in self.gdf_stores.geometry)

        grid['min_distance'] = grid.geometry.apply(min_distance)
        grid_sorted = grid.sort_values('min_distance', ascending=False)
        best = grid_sorted.iloc[0]
        top10 = grid_sorted.head(10)

        km = best['min_distance'] * 111.0
        print(f"ЛУЧШЕЕ МЕСТО: {best.geometry.y:.6f}, {best.geometry.x:.6f} → {km:.2f} км")
        return best, top10, grid

    def create_map(self, best_point, top_10, grid, bounds, water_gdf=None):
        center_lat = self.stores_data.latitude.mean()
        center_lon = self.stores_data.longitude.mean()
        m = folium.Map(location=[center_lat, center_lon], zoom_start=13, tiles="CartoDB positron")

        # Магазины
        for _, s in self.stores_data.iterrows():
            folium.Marker([s.latitude, s.longitude],
                          popup=s['name'],
                          icon=folium.Icon(color="red", icon="shopping-cart", prefix="fa")).add_to(m)

        # Лучшее место
        if best_point is not None:
            folium.Marker([best_point.geometry.y, best_point.geometry.x],
                          popup="<b>ЛУЧШЕЕ МЕСТО</b>",
                          icon=folium.Icon(color="green", icon="star", prefix="fa")).add_to(m)

        # Топ-10
        for i, (_, row) in enumerate(top_10.iloc[1:].iterrows(), 2):
            folium.CircleMarker([row.geometry.y, row.geometry.x],
                                radius=7, color="blue", fill=True,
                                popup=f"#{i} — {row['min_distance']*111:.2f} км").add_to(m)

        # Вода
        if water_gdf is not None and not water_gdf.empty:
            folium.GeoJson(water_gdf, style_function=lambda x: {
                'fillColor': 'lightblue', 'color': 'blue', 'weight': 1, 'fillOpacity': 0.5
            }).add_to(m)

        # Тепловая карта
        heat_data = [[p.y, p.x, d*111] for p, d in zip(grid.geometry, grid.min_distance)]
        HeatMap(heat_data, radius=12, blur=18).add_to(m)

        # Границы
        min_lon, min_lat, max_lon, max_lat = bounds
        folium.Polygon(locations=[[min_lat, min_lon], [min_lat, max_lon],
                                 [max_lat, max_lon], [max_lat, min_lon]],
                       color="black", weight=3, fill=False, dash_array="5").add_to(m)

        return m


# === УНИВЕРСАЛЬНАЯ ЗАГРУЗКА ВОДЫ ===
def get_water_bodies_safe(bounds):
    min_lon, min_lat, max_lon, max_lat = bounds
    south, west, north, east = min_lat, min_lon, max_lat, max_lon

    print("Попытка через osmnx.geometries_from_bbox...")
    try:
        water = ox.geometries_from_bbox(north, south, east, west, tags={'natural': 'water'})
        water = water[water.geometry.notnull()].to_crs("EPSG:4326")
        print(f"Успех! Найдено {len(water)} водоёмов")
        return clip_to_bounds(water, bounds)
    except:
        pass

    print("osmnx не сработал → используем Overpass API напрямую...")
    overpass_url = "https://overpass-api.de/api/interpreter"
    query = f"""
    [out:json][timeout:60];
    (
      way["natural"="water"]({south},{west},{north},{east});
      relation["natural"="water"]({south},{west},{north},{east});
    );
    out body;
    >;
    out skel qt;
    """
    try:
        response = requests.get(overpass_url, params={'data': query})
        data = response.json()
        print(f"Overpass вернул {len(data['elements'])} элементов")

        geometries = []
        for el in data['elements']:
            if el['type'] == 'way' and 'geometry' in el:
                coords = [(p['lon'], p['lat']) for p in el['geometry']]
                if len(coords) > 2:
                    geometries.append(Polygon(coords))
        if geometries:
            gdf = gpd.GeoDataFrame(geometry=geometries, crs="EPSG:4326")
            return clip_to_bounds(gdf, bounds)
    except Exception as e:
        print(f"Overpass тоже упал: {e}")

    print("не загружена way с полем geometry — Overpass не отдал нужные полигоны")
    return gpd.GeoDataFrame()


def clip_to_bounds(gdf, bounds):
    min_lon, min_lat, max_lon, max_lat = bounds
    clip_poly = box(min_lon, min_lat, max_lon, max_lat)
    clipped = gdf.copy()
    clipped.geometry = gdf.geometry.intersection(clip_poly)
    clipped = clipped[~clipped.geometry.is_empty]
    return clipped


def main():
    stores_df = pd.read_csv('existing_stores.csv')

    bounds = [
        30.1923419,
        59.92964418581805,
        30.2561469,
        59.95714219640026
    ]

    print("=== ПОИСК МАГАЗИНА В ПРЕДЕЛАХ ГРАНИЦ ===")
    water_gdf = get_water_bodies_safe(bounds)

    finder = StoreLocationFinder(stores_df)
    best, top10, grid = finder.find_best_location(bounds, grid_size=0.002, water_gdf=water_gdf)

    if best is None:
        print("Не удалось найти место на суше!")
        return

    # Сохранение
    results = [{'rank': i+1,
                'latitude': row.geometry.y,
                'longitude': row.geometry.x,
                'distance_km': round(row['min_distance']*111, 3)}
               for i, (_, row) in enumerate(top10.iterrows())]
    pd.DataFrame(results).to_csv('top_locations.csv', index=False)
    print("top_locations.csv сохранён")

    # Карта
    m = finder.create_map(best, top10, grid, bounds, water_gdf)
    m.save('best_store_location.html')
    print("Карта: best_store_location.html")

    print("\nТОП-5 наиболее удалённых мест:")
    for i, (_, row) in enumerate(top10.head(5).iterrows()):
        print(f"{i+1}. {row.geometry.y:.6f}, {row.geometry.x:.6f} → {row['min_distance']*111:.2f} км от новой точки до ближайшего уже существующего магазина")


if __name__ == "__main__":
    # Установка, если нет requests
    try:
        import requests
    except:
        print("Устанавливаем requests...")
        import subprocess, sys
        subprocess.check_call([sys.executable, "-m", "pip", "install", "requests"])
    main()
