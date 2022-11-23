import xml.etree.ElementTree as ET
ET.register_namespace("", "http://www.opengis.net/kml/2.2")

# Enter the file path to split
filename = "splitting_example/5GWh 18h.kml"

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

# Read in kml
tree = ET.parse(filename)

# Determine splitting
total_file_size = 0
file_no = 0
folder = tree.getroot()[0]
file_numbers = []
for document_no in range(len(folder)):
    total_file_size += size_of_document(folder[document_no])
    if total_file_size >= 10000000 * 0.9:
        file_no += 1
        total_file_size = 0
    file_numbers.append(file_no)

# to keep
file_numbers[0] = -1
file_numbers[1] = -1

# write to files
for i in range(file_no+1):
    tree = ET.parse(filename)
    root = tree.getroot()
    folder = root[0]
    for j in range(len(file_numbers)-1, -1, -1):
        if file_numbers[j] != i and file_numbers[j] != -1:
            folder.remove(folder[j])
    tree = ET.ElementTree(root)
    split_filename = filename[:-4]+str(i)+".kml"
    tree.write(split_filename, encoding='utf8', method='xml')

    # Ensure correct namespace
    with open(split_filename,'r',encoding='utf-8') as file:
        data = file.readlines()
    data[1] = "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">\n"
    with open(split_filename, 'w', encoding='utf-8') as file:
        file.writelines(data)
