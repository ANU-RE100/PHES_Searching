# run `pip install pyshp` if shapefile is not installed
import shapefile
from math import floor
from pathlib import Path

def main(path_to_shp, path_to_PHES, output_path):
    shp_path = Path(path_to_shp)
    output_folder = shp_path.parent/output_path
    if not output_folder.exists():
        output_folder.mkdir(parents=True)
    
    # Find existing DEMs
    existing_dems = []
    existing_dems_path = Path(path_to_PHES)/'task_lists/world_tasks.txt'
    with open(existing_dems_path) as file:
        for line in file:
            existing_dems.append(line.strip())

    sf = shapefile.Reader(path_to_shp)
    points = set()
    for shape in sf.shapes():
        for (lat, lon) in shape.points:
            # Add surrounding points
            for i in range(-1, 2):
                for j in range(-1, 2):
                    points.add(str(floor(lat)+i)+" "+str(floor(lon)+j))

    output_file = output_folder/(shp_path.stem+"_dems.txt")
    with open(output_file,'w') as file:
        for point in points:
            if point in existing_dems:            
                file.write(point+"\n")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Find DEMs of input shapefile')
    parser.add_argument('--shp_file_path', metavar='shp_file_path', required=True,
                        help='the path to the shp file (from working directory)')
    parser.add_argument('--PHES_path', metavar='PHES_path', required=True,
                        help='the path to the root of the PHES reop (from working directory)')
    parser.add_argument('--output_path', metavar='output_folder_path', default=".",
                        help='the path to the output folder relative to the shp file')

    args = parser.parse_args()
    main(path_to_shp=args.shp_file_path, path_to_PHES=args.PHES_path, output_path=args.output_path)
