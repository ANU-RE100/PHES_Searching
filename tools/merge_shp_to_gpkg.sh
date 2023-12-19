#!/bin/bash
# This Bash script is intended as an example of how to merge the Shapefiles downloaded
# using ./tools/osm_mine_downloader.py into a single GPKG

# create a new GeoPackage
ogr2ogr -f "GPKG" -nlt PROMOTE_TO_MULTI /path/to/output.gpkg /path/to/first/shapefile.shp

# append each additional shapefile as a new layer in the GeoPackage
for file in *.shp; do
if [ "$file" != "/path/to/first/shapefile.shp" ]
then
    ogr2ogr -f "GPKG" -update -append -nlt PROMOTE_TO_MULTI /path/to/output.gpkg $file -nln merged
fi
done

# Note that merged.gpkg will be apprxomiately 13GB in size.