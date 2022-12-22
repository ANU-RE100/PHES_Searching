
"""
To run:

python3 output_summary.py --task_file_path task_file_path --final_output_classes_path final_output_classes_path [--output_path output_folder_path]
"""
from tqdm import tqdm
import csv
from pathlib import Path
import pandas as pd
import xml.etree.cElementTree as ET
ET.register_namespace("", "http://www.opengis.net/kml/2.2")

FILESIZE = 1000 # In MB

def byte_len(string):
    """
    Return length in bytes of input string.
    """
    return len(bytes(string, 'utf8'))

def len_of_attrib(dictionary):
    """
    Returns length of input dictionary in bytes.
    Used for finding the length of attribute in a xml.
    """
    ret = 0
    for elm in dictionary:
        ret += byte_len(dictionary[elm]) + byte_len(elm) + 4 # 4 for "" and = and space
    return ret

def size_of_document(document, num_tabs=1):
    """
    Returns length of an xml document in bytes.
    """
    tag_length = byte_len(document.tag)
    if "}" in document.tag:
        tag_length = byte_len(document.tag.split("}")[1])
    size = len_of_attrib(document.attrib) + 7 + num_tabs + 2*tag_length # tag name is bounded above by 10
    if document.text is not None:
        size += byte_len(document.text)
    for child in document:
        size += size_of_document(child, num_tabs=num_tabs+1)
    return size

def get_site(task):
    task = task.split(" ")
    ns = "n" if int(task[-1]) >= 0 else "s"
    ew = "e" if int(task[-2]) >= 0 else "w"
    existing = "existing_" if len(task) == 3 else ""
    # TODO(work for pits as well)
    return existing+ns+str(abs(int(task[-1]))).zfill(2)+"_"+ew+str(abs(int(task[-2]))).zfill(3)

def main(path_to_tasks_file, path_to_final_output_classes, output_path="."):
    """
    Group kmls and create summary csvs.
    Files stored at output_path.
    """
    tasks_file_path = Path(path_to_tasks_file)

    final_output_classes_path = Path(path_to_final_output_classes)/"final_output_classes"

    output_folder = tasks_file_path.parent/output_path
    if not output_folder.exists():
        output_folder.mkdir(parents=True)

    csvs_folder = output_folder/"csvs"
    if not csvs_folder.exists():
        csvs_folder.mkdir(parents=True)

    kmls_folder = output_folder/"kmls"
    if not kmls_folder.exists():
        kmls_folder.mkdir(parents=True)

    # Read in tasks file
    with open(tasks_file_path) as f:
        tasks = f.readlines()

    data = []
    # Create total summary csv
    print("Generating summary csvs...")
    for task in tqdm(tasks):
        site = get_site(task)
        total_csv_file_name = final_output_classes_path/site/(site+"_total.csv")
        with open(total_csv_file_name, 'r') as f:
            final_row = f.readlines()[-1].strip().split(",")
            final_row.pop(3)
            final_row.pop(1)
            final_row = [site] + final_row
            data.append(final_row)

    with open(output_folder/"summary.csv", 'w+') as f:
        writer = csv.writer(f)
        writer.writerow(['Reservoir','Grid Identifier','Number of paired sites','Total potential capacity (GWh)'])
        for i in range(int(len(data))):
            writer.writerow(data[i])

    df = pd.read_csv(total_csv_file_name)
    sizes = list(df["Reservoir type"])[:-1]

    print("Grouping csvs and kmls...")
    for j, size in enumerate(tqdm(sizes)):
        data = []
        header = None
        output_kml = ET.Element("kml")
        document = ET.SubElement(output_kml, "Document", id="Layers")
        total_size = 0
        kml_file_no = 0
        new_kml_folder = None
        for task in tqdm(tasks, leave=False):
            site = get_site(task)
            total_csv_file_name = final_output_classes_path/site/(site+"_total.csv")
            total_csv_df = pd.read_csv(total_csv_file_name)

            # Should the file exist (volume is not 0)
            if (list(total_csv_df.iloc[j])[-1]) != 0:
                # Create summary csv
                csv_file_name = final_output_classes_path/site/(site+"_"+size+".csv")
                with open(csv_file_name, 'r') as f:
                    f_csv = csv.reader(f)
                    header = next(f_csv)
                    for row in f_csv:
                        data.append(row)

                # Create summary kml
                # Read in kml
                kml_file_name = final_output_classes_path/site/(site+"_"+size+".kml")
                tree = ET.parse(kml_file_name)
                root = tree.getroot()

                # Start new folder
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
                            total_size += size_of_document(placemark, num_tabs=3)
                            if placemark.tag.split("}")[1] == "Placemark":
                                new_kml_folder.insert(i, placemark)
                                i += 1
                
                if total_size > FILESIZE * 1000000:
                    output_tree = ET.ElementTree(output_kml)
                    grouped_file_path = kmls_folder/(size+"_summary_"+str(kml_file_no)+".kml")
                    output_tree.write(grouped_file_path, encoding='utf8', method='xml')
                    kml_file_no += 1
                    total_size = 0
                    new_kml_folder = None
                    output_kml = ET.Element("kml")
                    document = ET.SubElement(output_kml, "Document", id="Layers")

        if header is not None:
            with open(csvs_folder/(size+"_summary.csv"), 'w+') as f:
                writer = csv.writer(f)
                writer.writerow(header)
                for i in range(int(len(data))):
                    writer.writerow(data[i])

            output_tree = ET.ElementTree(output_kml)
            grouped_file_path = kmls_folder/(size+"_summary_"+str(kml_file_no)+".kml")
            output_tree.write(grouped_file_path, encoding='utf8', method='xml')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Split a KML file into <10MB files')
    parser.add_argument('--task_file_path', metavar='task_file_path', required=True,
                        help='the path to task file')    
    parser.add_argument('--final_output_classes_path', metavar='final_output_classes_path', required=True,
                        help='the path to final_output_classes folder')
    parser.add_argument('--output_path', metavar='output_folder_path', default=".",
                        help='the path to the output folder relative to the kml file')
    args = parser.parse_args()
    main(path_to_tasks_file=args.task_file_path, path_to_final_output_classes=args.final_output_classes_path, output_path=args.output_path)