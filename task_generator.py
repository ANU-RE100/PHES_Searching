import os
from shapely.geometry import Point, Polygon
import fiona

# Create a function to check if a point is in a polygon
def point_in_country(lat, lon, polygons):
    point = Point(lon, lat)
    return any(polygon.contains(point) for polygon in polygons)

# Function to convert the string into Polygon object
def create_polygon(data):
    coordinates = []
    split_data = data.split(',')
    for i in range(0, len(split_data), 2):
        coordinates.append((float(split_data[i]), float(split_data[i+1])))
    return Polygon(coordinates)

def list_files_in_folder(folder_path):
    file_list = []
    for filename in os.listdir(folder_path):
        file_list.append(filename)
    return file_list

def get_SRTM_filename(lat, lon):
    if lat < 0:
        lat_string = "s"+f"{-lat:02}"
    else:
        lat_string = "n" + f"{lat:02}"
            
    if lon < 0:
        lon_string = "w" + f"{-lon:03}"
    else:
        lon_string = "e" + f"{lon:03}"

    return lat_string + "_" + lon_string + "_1arc_v3.tif"

def get_FABDEM_filename(lat, lon):
    if lat < 0:
        lat_string = "S"+f"{-lat:02}"
    else:
        lat_string = "N" + f"{lat:02}"
            
    if lon < 0:
        lon_string = "W" + f"{-lon:03}"
    else:
        lon_string = "E" + f"{lon:03}"

    return lat_string + lon_string + "_FABDEM_V1-2.tif"

if __name__ == '__main__':
    DEM_type = input("SRTM or FABDEM task list? ")

    if DEM_type == "SRTM":
        folder_path = './input/DEMs'
    elif DEM_type == "FABDEM":
        folder_path = './input/FABDEMs'
    else:
        print("Invalid DEM type.")
        exit()

    DEM_list = list_files_in_folder(folder_path)

    file = open("task_lists/indonesia_tasks_fabdem.txt", "w+")

    task_type = input("Enter task type (e.g. ocean): ")

    input_type = input("COUNTRY or COORDINATE inputs? ")

    if input_type == "COUNTRY":
        country_name = input("Enter the country name: ")

        # Read the text file and create Polygon for Indonesia
        polygons = []
        with open("./input/countries/countries.txt", "r") as file:
            lines = file.readlines()
            country_found = False
            for line in lines:
                if line.startswith("@"+country_name):
                    country_found = True
                elif country_found and not (line.startswith("@")):
                    data = line.strip()
                    polygon = create_polygon(data)
                    polygons.append(polygon)
                elif country_found and line.startswith("@"):
                    break
        
        lowest_lat, lowest_lon, highest_lat, highest_lon = 90, 180, -90, -180
        count = 0
        for polygon in polygons:
            count += 1
            print(str(count) +" of "+str(len(polygons)))
            min_lat, min_lon, max_lat, max_lon = polygon.bounds
            lowest_lat = min(lowest_lat, min_lat)
            lowest_lon = min(lowest_lon, min_lon)
            highest_lat = max(highest_lat, max_lat)
            highest_lon = max(highest_lon, max_lon)

    elif input_type == "COORDINATE":
        lowest_lat = int(input("Enter lowest latitude: "))
        highest_lat = int(input("Enter highest latitude: "))
        lowest_lon = int(input("Enter lowest longitude: "))
        highest_lon = int(input("Enter highest longitude: "))

    else:
        file.close()
        print("Invalid input type.")
        exit()

    for lat in range(int(lowest_lat), int(highest_lat)+1):
        for lon in range(int(lowest_lon), int(highest_lon+1)):
            if DEM_type == "SRTM":
                DEM_name = get_SRTM_filename(lat,lon) 
            elif DEM_type == "FABDEM":
                DEM_name = get_FABDEM_filename(lat,lon)   
            else:
                print("Invalid DEM type.")
                exit()        

            if DEM_name in DEM_list:
                file.write(task_type+" "+str(lon)+" "+str(lat)+"\n")

    file.close()