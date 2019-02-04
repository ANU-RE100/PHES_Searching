file = open("tasks", "w+")

lowest_lat = int(input("Enter lowest latitude: "))
highest_lat = int(input("Enter highest latitude: "))
lowest_lon = int(input("Enter lowest longitude: "))
highest_lon = int(input("Enter highest longitude: "))

for lat in range(lowest_lat, highest_lat):
    for lon in range(lowest_lon, highest_lon):
        file.write(str(lon)+" "+str(lat)+"\n")
file.close()
