import osm2geojson
import geopandas as gpd
import requests
import json

# Bounding boxes
lat_bounds = range(-90,90)
long_bounds = range(-180,180,5)

count = 0
# Iterate over each region
for lat_lb in lat_bounds:

    # Lists to store data
    polys = []

    for long_lb in long_bounds:
        count+=1
        lat_ub = lat_lb + 1
        long_ub = long_lb + 5

        bounding_box = f"({lat_lb},{long_lb},{lat_ub},{long_ub})"

        # Overpass QL query for mines and quarries in current grid square
        query = f"""
        [out:json];
        (
          relation["industrial"="mine"]{bounding_box};
          relation["historic"="mine"]{bounding_box};
          relation["landuse"="quarry"]{bounding_box};
          way["natural"="water"]{bounding_box};
          way["historic"="mine"]{bounding_box};
          way["landuse"="quarry"]{bounding_box};
        );
        out body;
        >;
        out skel qt;
        """

        # Execute the query
        print(f"Executing Overpass API query for bounding box (Count = {count}): {bounding_box}...")
        url = "http://overpass-api.de/api/interpreter"
        r = requests.get(url, params={'data': query})
        if r.status_code != 200:
            print(f"Error: Server responded with status code {r.status_code}")
            print(f"Response content:\n{r.content}")
            continue
        try:
            data = r.json()
        except json.JSONDecodeError:
            print("Failed to parse JSON from response.")
            print(f"Response content:\n{r.content}")
            continue
        print("Query executed successfully!")

        # Convert the raw OSM data to GeoJSON
        geojson = osm2geojson.json2geojson(data)

        # Load the GeoJSON data into a GeoDataFrame
        gdf = gpd.read_file(json.dumps(geojson))

        # Only include Polygon and MultiPolygon geometries
        gdf = gdf[gdf.geometry.geom_type.isin(['Polygon', 'MultiPolygon'])]

        # Append the geometries to the polys list
        polys.extend(gdf.geometry)

    # Create a new GeoDataFrame with the combined geometries
    final_gdf = gpd.GeoDataFrame(geometry=polys)

    # Set the CRS to EPSG:4326
    final_gdf.crs = "EPSG:4326"

    # Save to Shapefile
    print("Saving data to Shapefile...")
    if (len(final_gdf.index) != 0):
        final_gdf.to_file(f"./preprocessing_files/osm_download/OSM_minesAndQuarries_{lat_lb}.shp")
        print("Shapefile saved successfully!")
    else:
        print(f"No polygons between latitudes {lat_lb} and {lat_ub}.")
    
