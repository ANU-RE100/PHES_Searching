from pathlib import Path
import json

DESCRIPTION = "As the proportion of wind and solar photovoltaics (PV) in an electrical grid extends into the 50-100% range a combination of additional long-distance high voltage transmission, demand management and local storage is required for stability. Pumped Hydro Energy Storage (PHES) constitutes 97% of electricity storage worldwide because of its low cost.<p>The <a href='http://re100.eng.anu.edu.au/global'>RE100 Group ANU</a> found about 530,000 potentially feasible PHES sites with storage potential of about 22 million Gigawatt-hours (GWh) by using geographic information system (GIS) analysis. This is about one hundred times greater than required to support a 100% global renewable electricity system. Brownfield sites (existing reservoirs, old mining sites) will be included in a future analysis.</p>This information has been developed by the 100% Renewable Energy group from the Research School of Electrical, Energy and Materials Engineering at the Australia National University.  http://re100.eng.anu.edu.au<p>Potential sites for off-river PHES are identified using GIS algorithms with defined search criteria. The surveyed latitude range is up to 60 degrees north and south. Each identified site comprises an upper and lower reservoir pair plus a hypothetical tunnel route between the reservoirs, and includes data such as latitude, longitude, altitude, head, slope, water volume, water area, rock volume, dam wall length, water/rock ratio, energy storage potential and approximate relative cost (categories A-E).</p><p>Wall heights are adjusted for each reservoir in a pair to yield equal water volumes to achieve the targeted energy storage. Energy (= head * usable volume * g * efficiency) and storage-length combinations are provided in Table 1. The approximate number of people that a reservoirs could service for a 100% renewable electricity grid is listed in the third line.</p><style>  table, th, td {    border: 1px solid black;    padding: 5px;    text-align: center;  }</style><center><table>  <tr>   <th>Energy</th> <th>Duration</th> <th>Millions of people</th> </tr> <tr>   <th>2 GWh</th> <td>6 hours</td> <td>0.1</td>   </tr>  <tr>    <th>5 GWh</th> <td>18 hours</td> <td>0.25</td>   </tr>  <tr>    <th>15 GWh</th> <td>18 hours</td> <td>0.75</td>   </tr>  <tr>    <th>50 GWh</th> <td>50 hours</td> <td>2.5</td>   </tr>  <tr>    <th>150 GWh</th> <td>50 hours</td> <td>7.5</td>   </tr>  <tr>    <th>500 GWh</th> <td>168 hours</td> <td>25</td>   </tr>  <tr>    <th>1500 GWh</th> <td>504 hours</td> <td>75</td>   </tr> </table></center><p>Virtually all upper reservoirs are away from rivers, and none intrude on protected area or urban areas listed in the databases that we use below.  There may be local constraints that prevent use of a particular site that is not reflected in these databases. Please refer to the ANU 100% Renewable Energy website for additional information: http://re100.eng.anu.edu.au/global</p>"


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


def create_member(site, country, phes_type):
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
        country + " " + phes_type + " " + str(energy) + "GWh " + str(time) + "h"
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


def main(path_to_simple_json, csvs_path, phes_type, country):
    assert phes_type in ["Greenfield", "Bluefield", "Brownfield", "Ocean"], (
        'PHES_type must be "Greenfield", "Bluefield", "Brownfield" or "Ocean" but was'
        + phes_type
    )
    with open(path_to_simple_json + "/simple.json", "r") as f:
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
            group["name"] = (
                "Global " if phes_type == "Greenfield" or phes_type == "Ocean" else ""
            ) + phes_type
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
                new_group["name"] = country + " " + phes_type + " Sites"
                new_group["description"] = DESCRIPTION
                new_group["members"] = []
                group["members"].append(new_group)
            group = new_group
        if group["members"] != []:
            val = input(
                "members already exist in "
                + group["name"]
                + ' and will be replaced, please type "Y" to contine\n'
            )
            if val.lower() != "y":
                quit()

        members = []
        path_to_csvs = Path(csvs_path) / "csvs"
        for pth in path_to_csvs.iterdir():
            members.append(create_member(pth.stem[:-8], country, phes_type))
        group["members"] = members

        # TODO(change to simple)
    with open(path_to_simple_json + "/simple.json", "w") as outfile:
        json.dump(simple_json, outfile, indent=2)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Write to simple.json")
    parser.add_argument(
        "--simple_json_path",
        metavar="simple_json_path",
        required=True,
        help="the path to simple.json (from working directory)",
    )
    parser.add_argument(
        "--csvs_path",
        metavar="csvs_dir_path",
        required=True,
        help="the path to the output dir generated from output_summary containing the csvs dir (from working directory)",
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
        path_to_simple_json=args.simple_json_path,
        csvs_path=args.csvs_path,
        phes_type=args.PHES_type,
        country=args.country,
    )
