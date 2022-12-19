
"""
To run:

python3 group_kmls.py --kml_file_path kml_file_path [--output_path output_folder_path] [--file_size output_file_size_MB]
"""
import os
from pathlib import Path
import xml.etree.cElementTree as ET
ET.register_namespace("", "http://www.opengis.net/kml/2.2")

def main(path_to_tasks_file, output_path="."):
    """
    Split KML file at path_to_kml_file into <sizeMB chunks.
    Files stored at output_path.
    """
    tasks_file_path = Path(path_to_tasks_file)

    output_folder = tasks_file_path.parent/output_path
    if not output_folder.exists():
        output_folder.mkdir(parents=True)

    # Read in tasks file
    with open(tasks_file_path) as f:
        lines = f.readlines()

    sizes = ['2.0GWh_6h','5.0GWh_18h','15GWh_18h','50GWh_50h','150GWh_50h','150GWh_168h','500GWh_168h','1500GWh_168h']

    for size in sizes:
        output_kml = ET.Element("kml")
        document = ET.SubElement(output_kml, "Document", id="Layers")
        add_styling = True
        for line in lines:
            task = line.split(" ")
            ns = "n" if int(task[-1]) >= 0 else "s"
            ew = "e" if int(task[-2]) >= 0 else "w"
            existing = "existing_" if len(task) == 3 else ""
            folder_name = existing+ns+str(abs(int(task[-1])))+"_"+ew+str(abs(int(task[-2])))
            csv_file_name = tasks_file_path.parent/"final_output_classes"/folder_name/(folder_name+"_"+size+".csv")
            # TODO(check if size is in csv and not 0)

            kml_file_name = tasks_file_path.parent/"final_output_classes"/folder_name/(folder_name+"_"+size+".kml")
            print(kml_file_name)
            if os.path.isfile(kml_file_name):
                print("found it!")
                tree = ET.parse(kml_file_name)
                root = tree.getroot()

                if add_styling:
                    add_styling = False
                    # Add styling
                    for i in range(4):
                        document.insert(i, root[0][i])
                    # Add name
                    ET.SubElement(document, "name").text = "name"

                i = 5

                for folder in root[0]:
                    if folder.tag.split("}")[1] == "Folder":
                        for placemark in folder:
                            if placemark.tag.split("}")[1] == "Placemark":
                                document.insert(i, placemark)
                                i += 1

        output_tree = ET.ElementTree(output_kml)
        
        grouped_file_path = output_folder/(size+"_summary.kml")
        print(grouped_file_path)
        output_tree.write(grouped_file_path, encoding='utf8', method='xml')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Split a KML file into <10MB files')
    parser.add_argument('--task_file_path', metavar='task_file_path', required=True,
                        help='the path to task file')
    parser.add_argument('--output_path', metavar='output_folder_path', default=".",
                        help='the path to the output folder relative to the kml file')
    args = parser.parse_args()
    main(path_to_tasks_file=args.task_file_path, output_path=args.output_path)