
"""
To run:

python3 group_kmls.py --task_file_path task_file_path [--output_path output_folder_path]
"""
import os
from pathlib import Path
import pandas as pd
import xml.etree.cElementTree as ET
ET.register_namespace("", "http://www.opengis.net/kml/2.2")

def get_site(task):
    task = task.split(" ")
    ns = "n" if int(task[-1]) >= 0 else "s"
    ew = "e" if int(task[-2]) >= 0 else "w"
    existing = "existing_" if len(task) == 3 else ""
    # TODO(work for pits as well)
    return existing+ns+str(abs(int(task[-1])))+"_"+ew+str(abs(int(task[-2])))

def main(path_to_tasks_file, output_path="."):
    """
    Group kmls and create summary csvs.
    Files stored at output_path.
    """
    tasks_file_path = Path(path_to_tasks_file)

    output_folder = tasks_file_path.parent/output_path
    if not output_folder.exists():
        output_folder.mkdir(parents=True)

    # Read in tasks file
    with open(tasks_file_path) as f:
        tasks = f.readlines()

    records = []
    sizes = []
    # Create total summary csv
    for task in tasks:
        site = get_site(task)
        total_csv_file_name = tasks_file_path.parent/"final_output_classes"/site/(site+"_total.csv")
        name = site[:-10]
        df = pd.read_csv(total_csv_file_name)
        record = list(df.iloc[-1])
        record.pop(3)
        record.pop(1)
        record = [name] + record
        records.append(record)
        if not sizes:
            for i in range(len(df.index) - 2, -1, -1):
                sizes.append(list(df.iloc[i])[1])
            sizes.reverse()

    grouped = pd.DataFrame(records,columns = ['Reservoir','Grid Identifier','Number of paired sites','Total potential capacity (GWh)'])
    grouped.to_csv(output_folder/("summary.csv"))

    for j, size in enumerate(sizes):
        records = []
        output_kml = ET.Element("kml")
        document = ET.SubElement(output_kml, "Document", id="Layers")
        new_kml_folder = None
        for task in tasks:
            site = get_site(task)
            total_csv_file_name = tasks_file_path.parent/"final_output_classes"/site/(site+"_total.csv")
            total_csv_df = pd.read_csv(total_csv_file_name)

            # Should the file exist (volume is not 0)
            if (list(total_csv_df.iloc[j])[-1]) != 0:
                # Create summary csv
                csv_file_name = tasks_file_path.parent/"final_output_classes"/site/(site+"_"+size+".csv")
                df = pd.read_csv(csv_file_name)
                for i in range(len(df.index) - 1):
                    records.append(list(df.iloc[i]))

                # Create summary kml
                kml_file_name = tasks_file_path.parent/"final_output_classes"/site/(site+"_"+size+".kml")
                tree = ET.parse(kml_file_name)
                root = tree.getroot()

                if new_kml_folder is None:
                    # Add styling
                    for i in range(4):
                        document.insert(i, root[0][i])
                    # Add name
                    ET.SubElement(document, "name").text = size
                    new_kml_folder = ET.SubElement(document, "Folder")

                i = 5
                for folder in root[0]:
                    if folder.tag.split("}")[1] == "Folder":
                        for placemark in folder:
                            if placemark.tag.split("}")[1] == "Placemark":
                                new_kml_folder.insert(i, placemark)
                                i += 1

        grouped = pd.DataFrame(records,columns = ['Pair Identifier','Class','Head (m)','Separation (km)','Slope (%)','Volume (GL)','Energy (GWh)','Storage time (h)','Combined water to rock ratio','Country','Non-overlapping','Upper Identifier','Upper elevation (m)','Upper latitude','Upper longitude','Upper reservoir area (ha)','Upper reservoir volume (GL)','Upper dam height (m)','Upper dam length (m)','Upper dam volume (GL)','Upper water to rock ratio','Upper country','Lower Identifier','Lower elevation (m)','Lower latitude','Lower longitude','Lower reservoir area (ha)','Lower reservoir volume (GL)','Lower dam height (m)','Lower dam length (m)','Lower dam volume (GL)','Lower water to rock ratio','Lower country'])
        grouped.to_csv(output_folder/(size+"_summary.csv"))

        output_tree = ET.ElementTree(output_kml)
        grouped_file_path = output_folder/(size+"_summary.kml")
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