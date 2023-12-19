import time
import requests

letters = ['E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
numbers = list(range(0,61))

""" letters = ['E']
numbers = [19] """

output_folder = "/media/fileshare/PHES_Searching/input/filters/WORLD_URBAN/"

links = []
for letter in letters:
    for number in numbers:
        number_str = f"{number:02d}"
        link_str = "https://sedac.ciesin.columbia.edu/arcgis/web/gmis-hbase/data-download/" + number_str + letter + "_hbase_human_built_up_and_settlement_extent_geographic_30m.zip"
        links.append(link_str)

# overriding requests.Session.rebuild_auth to mantain headers when redirected 
class SessionWithHeaderRedirection(requests.Session):
 
    AUTH_HOST = 'urs.earthdata.nasa.gov'
 
    def __init__(self, username, password):
 
        super().__init__()
 
        self.auth = (username, password)
 
  
 
   # Overrides from the library to keep headers when redirected to or from
 
   # the NASA auth host.
 
    def rebuild_auth(self, prepared_request, response):
 
        headers = prepared_request.headers
 
        url = prepared_request.url
 
  
 
        if 'Authorization' in headers:
 
            original_parsed = requests.utils.urlparse(response.request.url)
 
            redirect_parsed = requests.utils.urlparse(url)
 
  
 
            if (original_parsed.hostname != redirect_parsed.hostname) and redirect_parsed.hostname != self.AUTH_HOST and original_parsed.hostname != self.AUTH_HOST:
 
                del headers['Authorization']
 
  
 
        return
 
  
 
# create session with the user credentials that will be used to authenticate access to the data
 
username = "timothy.weber@anu.edu.au"
 
password= "Dirtygas%97"
 
session = SessionWithHeaderRedirection(username, password)
 
for link in links:
    while True:
        #print(link)
        try:
            # extract the filename from the URL
            filename = link.split("/")[-1]
            
            response = session.get(link, stream=True)
            if response.status_code == 200:
                with open(output_folder + filename, 'wb') as fd:
 
                    for chunk in response.iter_content(chunk_size=1024*1024):
            
                        fd.write(chunk)
                
                print(f"Successfully downloaded {filename}")
                time.sleep(30)
                break

            elif response.status_code == 404:
                print(f"404 for {filename}")
                time.sleep(5)
                break

            else:
                print("Request code for " + filename + ": " + str(response.status_code))
                time.sleep(30)

            
        except:
            print("Error occurred, waiting 120 seconds")
            time.sleep(120)