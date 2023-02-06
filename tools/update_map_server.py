from pathlib import Path
import json
import os
import requests

DESCRIPTION = "As the proportion of wind and solar photovoltaics (PV) in an electrical grid extends into the 50-100% range a combination of additional long-distance high voltage transmission, demand management and local storage is required for stability. Pumped Hydro Energy Storage (PHES) constitutes 97% of electricity storage worldwide because of its low cost.<p>The <a href='http://re100.eng.anu.edu.au/global'>RE100 Group ANU</a> found about 530,000 potentially feasible PHES sites with storage potential of about 22 million Gigawatt-hours (GWh) by using geographic information system (GIS) analysis. This is about one hundred times greater than required to support a 100% global renewable electricity system. Brownfield sites (existing reservoirs, old mining sites) will be included in a future analysis.</p>This information has been developed by the 100% Renewable Energy group from the Research School of Electrical, Energy and Materials Engineering at the Australia National University.  http://re100.eng.anu.edu.au<p>Potential sites for off-river PHES are identified using GIS algorithms with defined search criteria. The surveyed latitude range is up to 60 degrees north and south. Each identified site comprises an upper and lower reservoir pair plus a hypothetical tunnel route between the reservoirs, and includes data such as latitude, longitude, altitude, head, slope, water volume, water area, rock volume, dam wall length, water/rock ratio, energy storage potential and approximate relative cost (categories A-E).</p><p>Wall heights are adjusted for each reservoir in a pair to yield equal water volumes to achieve the targeted energy storage. Energy (= head * usable volume * g * efficiency) and storage-length combinations are provided in Table 1. The approximate number of people that a reservoirs could service for a 100% renewable electricity grid is listed in the third line.</p><style>  table, th, td {    border: 1px solid black;    padding: 5px;    text-align: center;  }</style><center><table>  <tr>   <th>Energy</th> <th>Duration</th> <th>Millions of people</th> </tr> <tr>   <th>2 GWh</th> <td>6 hours</td> <td>0.1</td>   </tr>  <tr>    <th>5 GWh</th> <td>18 hours</td> <td>0.25</td>   </tr>  <tr>    <th>15 GWh</th> <td>18 hours</td> <td>0.75</td>   </tr>  <tr>    <th>50 GWh</th> <td>50 hours</td> <td>2.5</td>   </tr>  <tr>    <th>150 GWh</th> <td>50 hours</td> <td>7.5</td>   </tr>  <tr>    <th>500 GWh</th> <td>168 hours</td> <td>25</td>   </tr>  <tr>    <th>1500 GWh</th> <td>504 hours</td> <td>75</td>   </tr> </table></center><p>Virtually all upper reservoirs are away from rivers, and none intrude on protected area or urban areas listed in the databases that we use below.  There may be local constraints that prevent use of a particular site that is not reflected in these databases. Please refer to the ANU 100% Renewable Energy website for additional information: http://re100.eng.anu.edu.au/global</p>"

# Custom exceptions.
class GeoserverException(Exception):
    def __init__(self, status, message):
        self.status = status
        self.message = message
        super().__init__(f"Status : {self.status} - {self.message}")

class Geoserver:
    """
    Attributes
    ----------
    service_url : str
        The URL for the GeoServer instance.
    username : str
        Login name for session.
    password: str
        Password for session.
    """

    def __init__(
        self,
        service_url: str = "http://localhost:8080/geoserver",  # default deployment url during installation
        username: str = "admin",  # default username during geoserver installation
        password: str = "geoserver",  # default password during geoserver installation
    ):
        self.service_url = service_url
        self.username = username
        self.password = password

        # private request method to reduce repetition of putting auth(username,password) in all requests call. DRY principle

    def _requests(self, method: str, url: str, **kwargs) -> requests.Response:
        if method == "post":
            return requests.post(url, auth=(self.username, self.password), **kwargs)
        elif method == "get":
            return requests.get(url, auth=(self.username, self.password), **kwargs)
        elif method == "put":
            return requests.put(url, auth=(self.username, self.password), **kwargs)
        elif method == "delete":
            return requests.delete(url, auth=(self.username, self.password), **kwargs)


    def workspace_exists(self, workspace):
        """
        get name  workspace if exist
        Example: curl -v -u admin:password -XGET -H "Accept: text/xml"  http://localhost:8080/geoserver/rest/workspaces/acme.xml
        """
        try:
            payload = {"recurse": "true"}
            url = "{}/rest/workspaces/{}.json".format(self.service_url, workspace)
            r = requests.get(url, auth=(self.username, self.password), params=payload)

            if r.status_code == 200:
                return r.json()
            else:
                return None

        except Exception as e:
            raise Exception(e)

    def create_workspace(self, workspace: str):
        """
        Create a new workspace in geoserver.
        The geoserver workspace url will be same as the name of the workspace.
        """
        try:
            url = "{}/rest/workspaces".format(self.service_url)
            data = "<workspace><name>{}</name></workspace>".format(workspace)
            headers = {"content-type": "text/xml"}
            r = self._requests("post", url, data=data, headers=headers)

            if r.status_code == 201:
                return "{} Workspace {} created!".format(r.status_code, workspace)
            else:
                raise GeoserverException(r.status_code, r.content)

        except Exception as e:
            raise Exception(e)

    def delete_workspace(self, workspace: str):
        """
        Parameters
        ----------
        workspace : str
        """
        try:
            payload = {"recurse": "true"}
            url = "{}/rest/workspaces/{}".format(self.service_url, workspace)
            r = requests.delete(
                url, auth=(self.username, self.password), params=payload
            )

            if r.status_code == 200:
                return "Status code: {}, delete workspace".format(r.status_code)

            else:
                raise GeoserverException(r.status_code, r.content)

        except Exception as e:
            raise Exception(e)

    def store_exists(self, store: str, workspace: str):
        """
        Return the data store in a given workspace.
        """
        try:
            url = "{}/rest/workspaces/{}/datastores/{}".format(
                self.service_url, workspace, store
            )

            r = self._requests("get", url)

            if r.status_code == 200:
                return r.json()
            else:
                return None

        except Exception as e:
            raise Exception(e)


    def create_store(
            self,
            store: str,
            workspace: str,
            db: str = "postgres",
            host: str = "localhost",
            port: int = 5432,
            pg_user: str = "postgres",
            pg_password: str = "admin",
        ):
            """
            Create PostGIS store for connecting postgres with geoserver.
            Parameters
            ----------
            store : str
            workspace : str
            db : str
            host : str
            port : int
            pg_user : str
            pg_password : str
            """
            try:
                url = "{}/rest/workspaces/{}/datastores".format(self.service_url, workspace)

                headers = {"content-type": "text/xml"}

                database_connection = """
                        <dataStore>
                        <name>{}</name>
                        <connectionParameters>
                        <entry key="host">{}</entry>
                        <entry key="port">{}</entry>
                        <entry key="user">{}</entry>
                        <entry key="passwd">{}</entry>
                        <entry key="dbtype">postgis</entry>
                        <entry key="database">{}</entry>
                        </connectionParameters>
                        </dataStore>
                        """.format(
                    store, host, port, pg_user, pg_password, db
                )

                r = self._requests(
                    "post",
                    url,
                    data=database_connection,
                    headers=headers,
                )

                if r.status_code in [200, 201]:
                    return "Featurestore created/updated successfully"
                else:
                    raise GeoserverException(r.status_code, r.content)

            except Exception as e:
                raise Exception(e)

    def layer_exists(self, layer_name: str, workspace: str):
        """
        Returns the layer by layer name.
        """
        try:
            url = "{}/rest/layers/{}".format(self.service_url, layer_name)
            if workspace is not None:
                url = "{}/rest/workspaces/{}/layers/{}".format(
                    self.service_url, workspace, layer_name
                )

            r = self._requests("get", url)
            if r.status_code == 200:
                return r.json()
            else:
                return None

        except Exception as e:
            raise Exception(e)

    def publish_layer(
        self,
        store: str,
        pg_table: str,
        workspace: str,
        *args
    ):
        """
        Parameters
        ----------
        store : str
        pg_table : str
        workspace : str
        Returns
        -------
        Notes
        -----
        Only user for postgis vector data
        input parameters: specify the name of the table in the postgis database to be published, specify the store,workspace name, and  the Geoserver user name, password and URL
        """
        try:
            url = "{}/rest/workspaces/{}/datastores/{}/featuretypes/".format(
                self.service_url, workspace, store
            )

            bounding_box = ""
            if args:
                bounding_box = ("<nativeBoundingBox><minx>{}</minx><maxx>{}</maxx><miny>{}</miny><maxy>{}</maxy></nativeBoundingBox><latLonBoundingBox><minx>{}</minx><maxx>{}</maxx><miny>{}</miny><maxy>{}</maxy></latLonBoundingBox>".format(
                    *args, *args
                )
            )

            layer_xml = (
                "<featureType><name>{}</name><title>{}</title>{}</featureType>".format(
                    pg_table, pg_table, bounding_box
                )
            )

            headers = {"content-type": "text/xml"}

            r = requests.post(
                url,
                data=layer_xml,
                auth=(self.username, self.password),
                headers=headers,
            )

            if r.status_code == 201:
                return r.status_code
            else:
                print("Have you run `write_to_postgres.sh`?")
                raise GeoserverException(r.status_code, r.content)

        except Exception as e:
            raise Exception(e)

    def style_layers(self, layer:str, workspace: str):
        try:
            url = "{}/rest/layers/{}:{}".format(
            self.service_url, workspace, layer
            )
            style = "<layer><defaultStyle><name>phes_generic</name></defaultStyle></layer>"
            headers = {"content-type": "text/xml"}

            r = requests.put(
                url,
                data=style,
                auth=(self.username, self.password),
                headers=headers,
            )

            if r.status_code == 200:
                return "Status code: {}, style layer".format(r.status_code)
            else:
                raise GeoserverException(r.status_code, r.content)
        except Exception as e:
            raise Exception(e)


    def delete_layer(self, layer_name: str, workspace: str):
        """
        Parameters
        ----------
        layer_name : str
        workspace : str
        """
        try:
            payload = {"recurse": "true"}
            url = "{}/rest/workspaces/{}/layers/{}".format(
                self.service_url, workspace, layer_name
            )
            r = self._requests(method="delete", url=url, params=payload)
            if r.status_code == 200:
                return "Status code: {}, delete layer".format(r.status_code)
            else:
                raise GeoserverException(r.status_code, r.content)

        except Exception as e:
            raise Exception(e)

    def extend_bounding_box(self, workspace: str, store: str, layer: str):
        """
        Parameters
        ----------
        workspace : str
        store_name : str
        layer: str
        """
        try:
            url = "{}/rest/workspaces/{}/datastores/{}/featuretypes/{}.json".format(
                self.service_url, workspace, store, layer
            )
            r = requests.get(url, auth=(self.username, self.password))
            if r.status_code == 200:
                bounding_box = r.json()["featureType"]["nativeBoundingBox"]
                self.delete_layer(layer, workspace)
                self.publish_layer(store, layer, workspace, bounding_box["minx"] - 1, bounding_box["maxx"] + 1, bounding_box["miny"] - 1, bounding_box["maxy"] + 1)
                return
            else:
                raise GeoserverException(r.status_code, r.content)

        except Exception as e:
            raise Exception(e)


def get_layer_name(site):
    if "." in site:
        l = site.find(".")
        site = site[:l] + site[l + 2 :]
    return site.lower()


def add_reservoir_description(dict, energy, time):
    dict["content"] = (
        "This information has been developed by the 100% Renewable Energy group from the Research School of Electrical, Energy and Materials Engineering at the Australia National University.  http://re100.eng.anu.edu.au<p>These sites will store "
        + str(energy)
        + " GWh of electricity which will discharge over "
        + str(time)
        + " hours yielding "
        + str(round(energy / time, 1))
        + " GW of power. Each pair of reservoirs would provide sufficient energy for "
        + "{:g}".format(float(energy * 0.05))
        + " million people. "
    )

def create_member(site, country, phes_type, name):
    member = {
        "type": "wms",
        "info": [
            {
                "name": "Description",
            },
            {
                "name": "Disclaimer",
                "content": "None of the PHES sites discussed in this study have been the subject of geological, hydrological, environmental, heritage and other studies, and it is not known whether any particular site would be suitable. The commercial feasibility of developing these sites is unknown. As with all major engineering projects, diligent attention to quality assurance would be required for safety and efficacy.<p>There has been no investigation of land tenure apart from exclusion of some environmental areas and urban areas, and no discussions with land owners and managers. Nothing in this list of potential site locations implies any rights for development of these locations. Accuracy of the sites depends on the accuracy of the source data. There may be sites that fall into local protected areas or urban areas that are not identified by the source data.</p>",
            },
            {
                "name": "Access and acknowledgements",
                "content": "In publications that use this information please acknowledge the RE100 Group, Australian National University, http://re100.eng.anu.edu.au/",
            },
            {
                "name": "Source Data",
                "content": "Digital Terrain Model (DTM) https://earthexplorer.usgs.gov/ <p>World Database Protected Areas: https://www.protectedplanet.net/</p> <p>Urban extent:  HBASE http://sedac.ciesin.columbia.edu/data/set/ulandsat-hbase-v1/data-download</p> ",
            },
        ],
        "infoSectionOrder": [
            "Disclaimer",
            "Description",
            "Access and acknowledgements",
            "Source Data",
        ],
        "opacity": 1,
        "featureInfoTemplate": {
            "template": "<div>  <style>    tr {background-color: transparent; ! important;}  </style>  {{description}}</div>"
        },
        "tileErrorHandlingOptions": {"ignoreUnknownTileErrors": True},
    }
    country = country or "Global"
    energy = int(site.split("_")[0].split("G")[0].split(".")[0])
    time = int(site.split("_")[1].split("h")[0])
    add_reservoir_description(member["info"][0], energy, time)
    member["name"] = (
        (name or country + " " + phes_type) + " " + str(energy) + "GWh " + str(time) + "h"
    )
    member["id"] = country + "_" + phes_type + "_" + site
    member["layers"] = get_layer_name(site)
    member["url"] = (
        "https://re100.anu.edu.au/geoserver/"
        + country.lower()
        + "_"
        + phes_type.lower()
        + "/wms"
    )

    return member


def main(output_summary_path, postgres_pw, admin_pw, re100_un, phes_type, country, name):
    assert phes_type in ["Greenfield", "Bluefield", "Brownfield", "Ocean"], (
        'PHES_type must be "Greenfield", "Bluefield", "Brownfield" or "Ocean" but was'
        + phes_type
    )
    if phes_type in ["Bluefield", "Brownfield"] and not country:
        var = input("You have not specifed a country for a " + phes_type + " search. Would you like to continue? [y/n]\n")
        if var.lower() != "y":
            quit()

    # UPDATE GEOSERVER

    geo = Geoserver("https://re100.anu.edu.au/geoserver", username="admin", password=admin_pw)

    workspace = country.lower() + "_" + phes_type.lower()
    store = workspace + "_store"
    if geo.workspace_exists(workspace=workspace):
        var = input("workspace " + workspace + " already exists and will be replaced. If you do not wish to delete please type \'N\', otherwise type any other key.\n")
        if var.lower() == "n":
            quit()
        geo.delete_workspace(workspace=workspace)
    geo.create_workspace(workspace=workspace)

    assert not geo.store_exists(store=store, workspace=workspace)
    geo.create_store(store=store, workspace=workspace, db=workspace, pg_user='postgres', pg_password=postgres_pw)

    path_to_csvs = Path(output_summary_path) / "csvs"
    for pth in path_to_csvs.iterdir():
        layer = get_layer_name(pth.stem[:-8])
        print("publishing", layer)
        assert not geo.layer_exists(layer, workspace=workspace)
        geo.publish_layer(workspace=workspace, store=store, pg_table=layer)
        geo.extend_bounding_box(workspace=workspace, store=store, layer=layer)
        geo.style_layers(layer=layer, workspace=workspace)


# UPDATE SIMPLE.JSON
    
    # Copy locally
    print("Copying locally:")
    os.system("scp "+re100_un+"@re100.anu.edu.au:/etc/var/www/RE100Map/wwwroot/init/simple.json .")

    with open("./simple.json", "r") as f:
        simple_json = json.load(f)

        group = None
        for g in simple_json["catalog"]:
            if g["id"] == phes_type:
                group = g
                break
        if not group:
            group = {}
            group["id"] = phes_type
            group["type"] = "group"
            group["name"] = name or "Global " + phes_type if not country else phes_type
            group["description"] = DESCRIPTION
            group["members"] = []
            simple_json["catalog"].append(group)
        if country:
            new_group = None
            for g in group["members"]:
                if g["id"] == phes_type + "_" + country:
                    new_group = g
                    break
            if not new_group:  # country does not exist, create new group
                new_group = {}
                new_group["type"] = "group"
                new_group["id"] = phes_type + "_" + country
                new_group["name"] = (name or country + " " + phes_type) + " Sites"
                new_group["description"] = DESCRIPTION
                new_group["members"] = []
                group["members"].append(new_group)
            group = new_group
        if group["members"] != []:
            var = input(
                "members already exist in "
                + group["name"]
                + " and will be replaced. If you do not wish to delete please type \'N\', otherwise type any other key.\n"
            )
            if var.lower() == "n":
                quit()

        members = []
        path_to_csvs = Path(output_summary_path) / "csvs"
        for pth in path_to_csvs.iterdir():
            members.append(create_member(pth.stem[:-8], country, phes_type, name))
        group["members"] = members

    with open("./simple.json", "w") as outfile:
        json.dump(simple_json, outfile, indent=2)

    # Copy to re100
    print("Copying to RE100:")
    os.system("scp ./simple.json "+re100_un+"@re100.anu.edu.au:/etc/var/www/RE100Map/wwwroot/init")

    # Delete local copy
    os.remove("./simple.json")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Write to simple.json")
    parser.add_argument(
        "--output_summary_path",
        metavar="output_summary_path",
        required=True,
        help="the path to the output directory generated from output_summary.py (from working directory)",
    )
    parser.add_argument(
        "--postgres_password",
        metavar="postgres_password",
        required=True,
        help="Password for postgres user",
    )
    parser.add_argument(
        "--admin_password",
        metavar="admin_password",
        required=True,
        help="Password for geoserver admin user",
    )
    parser.add_argument(
        "--re100_username",
        metavar="re100_username",
        required=True,
        help="Password for geoserver admin user",
    )
    parser.add_argument(
        "--PHES_type",
        metavar="PHES_type",
        required=True,
        help='"Greenfield", "Bluefield", "Brownfield" or "Ocean"',
    )
    parser.add_argument(
        "--country",
        metavar="country",
        help="If specified, this will appear as a subsection to the PHES_type",
    )
    parser.add_argument(
        "--name",
        metavar="name",
        help="Name of data to be displayed on map server, if not specified will be `country PHES_type` (or `Global PHES_type` if country not specified)",
    )
    args = parser.parse_args()
    main(
        output_summary_path=args.output_summary_path,
        postgres_pw=args.postgres_password,
        admin_pw=args.admin_password,
        re100_un=args.re100_username,
        phes_type=args.PHES_type,
        country=args.country,
        name=args.name  
    )
