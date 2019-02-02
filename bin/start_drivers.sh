#!/bin/bash

ndrivers=$1
me=`hostname`

set -e
for ((task=0;task<ndrivers;task++))
do
    bin/search_driver &
done

echo Started $ndrivers searchers on $me