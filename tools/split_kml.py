"""
To run:

python3 split_kml.py --kml_file_path kml_file_path [--output_path output_folder_path] [--file_size output_file_size_MB]
"""

from pathlib import Path
import xml.etree.ElementTree as ET
ET.register_namespace("", "http://www.opengis.net/kml/2.2")

FILESIZE = 10 # In MB

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

def main(path_to_kml_file, output_path=".", size=FILESIZE):
    """
    Split KML file at path_to_kml_file into <sizeMB chunks.
    Files stored at output_path.
    """
    kml_file_path = Path(path_to_kml_file)
    output_folder = kml_file_path.parent/output_path
    if not output_folder.exists():
        output_folder.mkdir(parents=True)
    kml_file_name = kml_file_path.name
    # Read in kml
    tree = ET.parse(path_to_kml_file)

    # Determine splitting
    total_file_size = 0
    file_no = 0
    folder = tree.getroot()[0]
    file_numbers = []
    for document_no in range(len(folder)):
        total_file_size += size_of_document(folder[document_no])
        if total_file_size >= size * 1000000 * 0.9:
            file_no += 1
            total_file_size = 0
        if folder[document_no].tag.split("}")[1] != "Document":
            file_numbers.append(-1)   # We want these lines in every file
        else:
            file_numbers.append(file_no)

    # write to files
    for i in range(file_no+1):
        tree = ET.parse(path_to_kml_file)
        root = tree.getroot()
        folder = root[0]
        for j in range(len(file_numbers)-1, -1, -1):
            if file_numbers[j] != i and file_numbers[j] != -1:
                folder.remove(folder[j])
        tree = ET.ElementTree(root)
        split_file_path = output_folder/(kml_file_name[:-4]+"_"+str(i)+".kml")
        tree.write(split_file_path, encoding='utf8', method='xml')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Split a KML file into <10MB files')
    parser.add_argument('--kml_file_path', metavar='kml_file_path', required=True,
                        help='the path to kml file')
    parser.add_argument('--output_path', metavar='output_folder_path', default=".",
                        help='the path to the output folder relative to the kml file')
    parser.add_argument('--file_size', metavar='output_file_size_(MB)', default=FILESIZE,
                        help='the file size in MB to split the kml into')
    args = parser.parse_args()
    main(path_to_kml_file=args.kml_file_path, output_path=args.output_path, size=int(args.file_size))