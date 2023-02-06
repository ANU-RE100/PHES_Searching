from pathlib import Path
import json
import os
import requests
from update_map_server import Geoserver

DESCRIPTION = "As the proportion of wind and solar photovoltaics (PV) in an electrical grid extends into the 50-100% range a combination of additional long-distance high voltage transmission, demand management and local storage is required for stability. Pumped Hydro Energy Storage (PHES) constitutes 97% of electricity storage worldwide because of its low cost.<p>The <a href='http://re100.eng.anu.edu.au/global'>RE100 Group ANU</a> found about 530,000 potentially feasible PHES sites with storage potential of about 22 million Gigawatt-hours (GWh) by using geographic information system (GIS) analysis. This is about one hundred times greater than required to support a 100% global renewable electricity system. Brownfield sites (existing reservoirs, old mining sites) will be included in a future analysis.</p>This information has been developed by the 100% Renewable Energy group from the Research School of Electrical, Energy and Materials Engineering at the Australia National University.  http://re100.eng.anu.edu.au<p>Potential sites for off-river PHES are identified using GIS algorithms with defined search criteria. The surveyed latitude range is up to 60 degrees north and south. Each identified site comprises an upper and lower reservoir pair plus a hypothetical tunnel route between the reservoirs, and includes data such as latitude, longitude, altitude, head, slope, water volume, water area, rock volume, dam wall length, water/rock ratio, energy storage potential and approximate relative cost (categories A-E).</p><p>Wall heights are adjusted for each reservoir in a pair to yield equal water volumes to achieve the targeted energy storage. Energy (= head * usable volume * g * efficiency) and storage-length combinations are provided in Table 1. The approximate number of people that a reservoirs could service for a 100% renewable electricity grid is listed in the third line.</p><style>  table, th, td {    border: 1px solid black;    padding: 5px;    text-align: center;  }</style><center><table>  <tr>   <th>Energy</th> <th>Duration</th> <th>Millions of people</th> </tr> <tr>   <th>2 GWh</th> <td>6 hours</td> <td>0.1</td>   </tr>  <tr>    <th>5 GWh</th> <td>18 hours</td> <td>0.25</td>   </tr>  <tr>    <th>15 GWh</th> <td>18 hours</td> <td>0.75</td>   </tr>  <tr>    <th>50 GWh</th> <td>50 hours</td> <td>2.5</td>   </tr>  <tr>    <th>150 GWh</th> <td>50 hours</td> <td>7.5</td>   </tr>  <tr>    <th>500 GWh</th> <td>168 hours</td> <td>25</td>   </tr>  <tr>    <th>1500 GWh</th> <td>504 hours</td> <td>75</td>   </tr> </table></center><p>Virtually all upper reservoirs are away from rivers, and none intrude on protected area or urban areas listed in the databases that we use below.  There may be local constraints that prevent use of a particular site that is not reflected in these databases. Please refer to the ANU 100% Renewable Energy website for additional information: http://re100.eng.anu.edu.au/global</p>"


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
        (name or country + " " + phes_type)
        + " "
        + str(energy)
        + "GWh "
        + str(time)
        + "h"
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


def main(
    admin_pw, re100_un, heatmap_type, heatmap_cost, transmission_type, country, name
):
    assert heatmap_type in ["Solar", "Wind"], (
        'heatmap_type must be "Solar" or "Wind" but was' + heatmap_type
    )
    assert heatmap_cost in ["Low", "Medium", "High"], (
        'heatmap_cost must be "Low", "Medium" or "High" but was' + heatmap_type
    )
    assert transmission_type in ["Overhead", "Underground"], (
        'transmission_type must be "Overhead" or "Underground" but was' + heatmap_type
    )
    if not country:
        var = input(
            "You have not specifed a country for a "
            + heatmap_type
            + " search. Would you like to continue? [y/n]\n"
        )
        if var.lower() != "y":
            quit()

    # UPDATE GEOSERVER
    geo = Geoserver(
        "https://re100.anu.edu.au/geoserver", username="admin", password=admin_pw
    )

    workspace = "heatmaps"
    layer = (
        country.lower()
        + "_"
        + transmission_type.lower()
        + "_"
        + heatmap_type.lower()
        + "_"
        + heatmap_cost.lower()
        + "-cost_heatmap"
    )
    store = layer + "_store"

    assert not geo.store_exists(store=store, workspace=workspace)
    tif_path = f"/mnt/data/Heatmap/{heatmap_type} {transmission_type.lower()} {heatmap_cost.lower()}-cost/Australia Raster/Australia_raster.tif"
    geo.create_raster_store(
        store=store,
        workspace=workspace,
        tif_path=tif_path,
    )
    print("publishing", layer)

    assert not geo.layer_exists(layer, workspace=workspace)
    geo.publish_raster_layer(
        workspace=workspace, store=store, layer_name=layer, tif_path=tif_path
    )
    geo.style_heatmap_layers(layer=layer, workspace=workspace)

    # UPDATE SIMPLE.JSON

    # Copy locally
    print("Copying locally:")
    os.system(
        "scp "
        + re100_un
        + "@re100.anu.edu.au:/etc/var/www/RE100Map/wwwroot/init/simple.json ."
    )

    with open("./simple.json", "r") as f:
        simple_json = json.load(f)

        group = None
        for g in simple_json["catalog"]:
            if g["id"] == heatmap_type:
                group = g
                break
        if not group:
            group = {}
            group["id"] = heatmap_type
            group["type"] = "group"
            group["name"] = (
                name or "Global " + heatmap_type if not country else heatmap_type
            )
            group["description"] = DESCRIPTION
            group["members"] = []
            simple_json["catalog"].append(group)
        if country:
            new_group = None
            for g in group["members"]:
                if g["id"] == heatmap_type + "_" + country:
                    new_group = g
                    break
            if not new_group:  # country does not exist, create new group
                new_group = {}
                new_group["type"] = "group"
                new_group["id"] = heatmap_type + "_" + country
                new_group["name"] = (name or country + " " + heatmap_type) + " Sites"
                new_group["description"] = DESCRIPTION
                new_group["members"] = []
                group["members"].append(new_group)
            group = new_group
        if group["members"] != []:
            var = input(
                "members already exist in "
                + group["name"]
                + " and will be replaced. If you do not wish to delete please type 'N', otherwise type any other key.\n"
            )
            if var.lower() == "n":
                quit()
        exit()

        members = []
        members.append(create_member(pth.stem[:-8], country, heatmap_type, name))
        group["members"] = members

    with open("./simple.json", "w") as outfile:
        json.dump(simple_json, outfile, indent=2)

    # Copy to re100
    print("Copying to RE100:")
    os.system(
        "scp ./simple.json "
        + re100_un
        + "@re100.anu.edu.au:/etc/var/www/RE100Map/wwwroot/init"
    )

    # Delete local copy
    os.remove("./simple.json")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Write to simple.json")
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
        "--heatmap_type",
        metavar="heatmap_type",
        required=True,
        help='"Solar" or "Wind"',
    )
    parser.add_argument(
        "--heatmap_cost",
        metavar="heatmap_cost",
        required=True,
        help='"Low", "Medium" or "High"',
    )
    parser.add_argument(
        "--transmission_type",
        metavar="transmission_type",
        required=True,
        help='"Overhead" or "Underground"',
    )
    parser.add_argument(
        "--country",
        metavar="country",
        help="If specified, this will appear as a subsection",
    )
    parser.add_argument(
        "--name",
        metavar="name",
        help="Name of data to be displayed on map server, if not specified will be `country PHES_type` (or `Global PHES_type` if country not specified)",
    )
    args = parser.parse_args()
    main(
        admin_pw=args.admin_password,
        re100_un=args.re100_username,
        heatmap_type=args.heatmap_type,
        heatmap_cost=args.heatmap_cost,
        transmission_type=args.transmission_type,
        country=args.country,
        name=args.name,
    )
