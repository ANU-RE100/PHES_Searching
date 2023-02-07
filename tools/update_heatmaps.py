from pathlib import Path
import json
import os
import requests
from update_map_server import Geoserver

DESCRIPTION = """
Solar and wind play major roles in decarbonising the energy system. But where are the best sites?

In our high-resolution heat maps, an indicative cost of electricity (in AUD/MWh) is calculated for each pixel (1km x 1km for solar and 250m x 250m for wind), comprising the cost of energy from a solar/wind farm PLUS an associated powerline connecting the solar/wind farm to the existing and planned high voltage transmission network.

Access to transmission is currently the largest constraint for solar/wind farms. Heat maps identify the potential of all prospective areas close to the transmission network and are much more useful than solar and wind data alone.

Our heat maps provide detailed qualitative and quantitative information about prospective places for solar and wind farms. They balance knowledge between developers and landholders. They empower landholders and local Government to identify prospective area and facilitate collective bargaining with developers.

The maps use colours to show the lowest indicative costs for wind or solar. The example below shows the area around Tamworth for the low-cost, overhead transmission line scenario for solar farms. The lowest indicative cost is in red, pink shows the second lowest costs, and the unsuitable areas within urban areas and the surrounding National Parks are shown in green.
"""


def add_description(dict):
    dict[
        "content"
    ] = "This information has been developed by the 100% Renewable Energy group from the Research School of Electrical, Energy and Materials Engineering at the Australia National University.  http://re100.eng.anu.edu.au"


def create_member(country, heatmap_type, heatmap_cost, transmission_type, name):
    member = {
        "type": "wms",
        "info": [
            {
                "name": "Description",
            },
            {
                "name": "Disclaimer",
                "content": """None of the solar or wind sites discussed in this study have been the subject of any other studies, and it is not known whether any particular site would be suitable. The commercial feasibility of developing these sites is unknown. As with all major engineering projects, diligent attention to quality assurance would be required for safety and efficacy. The publicly available datasets used for this work contain inaccuracies – results should be viewed as indicative.

There has been no investigation of land tenure apart from exclusion of some environmental areas and urban areas (which are marked as “Unsuitable” and assigned a value of 0 in the maps), and no discussions with land owners and managers. The heat maps do not imply any rights for development at any location. Accuracy of the sites depends on the accuracy of the source data. The indicative costs are assumed to be constant within a pixel, which is 1km x 1km for solar heat maps and 250m x 250m for wind heat maps. There may be pixels that are not marked as “Unsuitable” in the heat maps but fall into protected areas or urban areas that are not identified by the source data.

Only transmission lines with rating 275kV and higher were included except for Tasmania (220 kV).
Transmission constraints were not included in the study. Substation costs were not included in
indicative costs.""",
            },
            {
                "name": "Access and acknowledgements",
                "content": "In publications that use this information please acknowledge the RE100 Group, Australian National University, http://re100.eng.anu.edu.au/",
            },
            # {
            # "name": "Source Data",
            # "content": "",
            # },
        ],
        "infoSectionOrder": [
            "Disclaimer",
            "Description",
            "Access and acknowledgements",
            # "Source Data",
        ],
        "opacity": 0.7,
        # "featureInfoTemplate": {
        # "template": "<div>  <style>    tr {background-color: transparent; ! important;}  </style>  {{description}}</div>"
        # },
        "tileErrorHandlingOptions": {"ignoreUnknownTileErrors": True},
    }
    country = country or "Global"
    add_description(member["info"][0])
    member["name"] = (
        name
        or country
        + " "
        + heatmap_type.lower()
        + " "
        + transmission_type.lower()
        + " "
        + heatmap_cost.lower()
        + "-cost"
    )

    member["id"] = (
        country.lower()
        + "_"
        + transmission_type.lower()
        + "_"
        + heatmap_type.lower()
        + "_"
        + heatmap_cost.lower()
        + "-cost_heatmap"
    )
    member["layers"] = member["id"]
    member["url"] = "https://re100.anu.edu.au/geoserver/heatmaps/wms"

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
            if g["id"] == "heatmaps":
                group = g
                break
        if not group:
            group = {}
            group["id"] = "heatmaps"
            group["type"] = "group"
            group["name"] = "Heatmaps"
            group["description"] = DESCRIPTION
            group["members"] = []
            simple_json["catalog"].append(group)

        type_group = None
        for g in group["members"]:
            if g["id"] == heatmap_type:
                type_group = g
                break
        if not type_group:
            type_group = {}
            type_group["id"] = heatmap_type.lower()
            type_group["type"] = "group"
            type_group["name"] = heatmap_type
            type_group["description"] = DESCRIPTION
            type_group["members"] = []
            group["members"].append(type_group)

        if country:
            new_group = None
            for g in type_group["members"]:
                if g["id"] == heatmap_type + "_" + country:
                    new_group = g
                    break
            if not new_group:  # country does not exist, create new group
                new_group = {}
                new_group["type"] = "group"
                new_group["id"] = heatmap_type + "_" + country
                new_group["name"] = country
                new_group["description"] = DESCRIPTION
                new_group["members"] = []
                type_group["members"].append(new_group)
            group = new_group
        if group["members"] != []:
            var = input(
                "members already exist in "
                + group["name"]
                + " and will be replaced. If you do not wish to delete please type 'N', otherwise type any other key.\n"
            )
            if var.lower() == "n":
                quit()

        cid = (
            country.lower()
            + "_"
            + transmission_type.lower()
            + "_"
            + heatmap_type.lower()
            + "_"
            + heatmap_cost.lower()
            + "-cost_heatmap"
        )
        group["members"] = list(filter(lambda g: g["id"] == cid, group["members"]))

        group["members"].append(
            create_member(country, heatmap_type, heatmap_cost, transmission_type, name)
        )

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
