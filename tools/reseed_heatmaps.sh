#!/bin/bash

password="$1"

for htype in solar wind; do
  for hcost in low medium high; do
    for ttype in overhead underground; do
      curl -v -u admin:"$password" \
        "https://re100.anu.edu.au/geoserver/gwc/rest/seed/heatmaps:australia_${ttype}_${htype}_${hcost}-cost_heatmap" \
    --data-raw 'threadCount=01&type=reseed&gridSetId=WebMercatorQuad&tileFormat=image%2Fpng&zoomStart=00&zoomStop=10&parameter_STYLES=heatmaps%3AHeatmap&minX=&minY=&maxX=&maxY=&tileFailureRetryCount=1&tileFailureRetryWaitTime=100&totalFailuresBeforeAborting=1000' \
    --compressed
    done
  done
done
