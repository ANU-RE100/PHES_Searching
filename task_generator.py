import os

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
    DEM_type = input("SRTM or FABDEM task list?")

    if DEM_type == "SRTM":
        folder_path = './input/DEMs'
    elif DEM_type == "FABDEM":
        folder_path = './input/FABDEMs'
    else:
        print("Invalid DEM type.")
        exit()

    DEM_list = list_files_in_folder(folder_path)

    file = open("task_lists/world_tasks.txt", "w+")

    lowest_lat = int(input("Enter lowest latitude: "))
    highest_lat = int(input("Enter highest latitude: "))
    lowest_lon = int(input("Enter lowest longitude: "))
    highest_lon = int(input("Enter highest longitude: "))

    task_type = input("Enter task type (e.g. ocean): ")

    for lat in range(lowest_lat, highest_lat+1):
        for lon in range(lowest_lon, highest_lon+1):
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