from pathlib import Path
import requests

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


def main(csvs_path, postgres_pw, admin_pw, phes_type, country):
    assert phes_type in ["Greenfield", "Bluefield", "Brownfield", "Ocean"], (
        'PHES_type must be "Greenfield", "Bluefield", "Brownfield" or "Ocean" but was'
        + phes_type
    )

    geo = Geoserver("https://re100.anu.edu.au/geoserver", username="admin", password=admin_pw)

    workspace = country.lower() + "_" + phes_type.lower()
    store = workspace + "_store"
    print(workspace)
    if geo.workspace_exists(workspace=workspace):
        var = input("workspace " + workspace + " already exists and will be deleted. If you do not wish to delete please type \'N\', otherwise type any other key.\n")
        if var.lower() == "n":
            quit()
        geo.delete_workspace(workspace=workspace)
    geo.create_workspace(workspace=workspace)

    assert not geo.store_exists(store=store, workspace=workspace)
    geo.create_store(store=store, workspace=workspace, db=workspace, pg_user='postgres', pg_password=postgres_pw)

    def get_layer_name(site):
        if "." in site:
            l = site.find(".")
            site = site[:l] + site[l + 2 :]
        return site.lower()

    path_to_csvs = Path(csvs_path) / "csvs"
    for pth in path_to_csvs.iterdir():
        layer = get_layer_name(pth.stem[:-8])
        print("publishing", layer)
        assert not geo.layer_exists(layer, workspace=workspace)
        geo.publish_layer(workspace=workspace, store=store, pg_table=layer)
        geo.extend_bounding_box(workspace=workspace, store=store, layer=layer)
        geo.style_layers(layer=layer, workspace=workspace)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create layers")
    parser.add_argument(
        "--csvs_path",
        metavar="csvs_dir_path",
        required=True,
        help="the path to the output dir generated from output_summary containing the csvs dir (from working directory)",
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
    args = parser.parse_args()
    main(
        csvs_path=args.csvs_path,
        postgres_pw=args.postgres_password,
        admin_pw=args.admin_password,
        phes_type=args.PHES_type,
        country=args.country
    )
# ./write_to_postgres.sh path_to_output_folder database postgres_password
