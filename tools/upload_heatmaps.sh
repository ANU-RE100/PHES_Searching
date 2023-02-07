#!/bin/bash

password="$1"

for htype in Solar Wind; do
  for ttype in Overhead Underground; do
    for hcost in Low Medium High; do
      python update_heatmaps.py --admin_password "$password" --re100_username u6311272 --heatmap_type $htype --heatmap_cost $hcost --transmission_type $ttype --country Australia
    done
  done
done
