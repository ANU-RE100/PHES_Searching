#!/bin/bash

password="$1"

sizes="2gwh_6h 5gwh_18h 15gwh_18h 50gwh_50h 150gwh_50h 500gwh_168h 1500gwh_504h"
#types="global_bluefield global_greenfield australia_bluefield global_ocean indonesia_bluefield philippines_bluefield bhutan_bluefield asean_bluefield"
types="global_brownfield"
#types="global_greenfield"

for type in $types; do
  for size in $sizes; do
      curl -v -u admin:"$password" \
        "https://re100.anu.edu.au/geoserver/gwc/rest/seed/$type:$size" \
    --data-raw \
          'threadCount=01&type=reseed&gridSetId=WebMercatorQuad&tileFormat=image%2Fpng&zoomStart=00&zoomStop=12&parameter_STYLES=phes_generic&minX=&minY=&maxX=&maxY=&tileFailureRetryCount=1&tileFailureRetryWaitTime=100&totalFailuresBeforeAborting=1000' \
    --compressed
  done
done
