#!/bin/bash

export PHES_BINDIR=bin

ndrivers=$1
tasks_file=$2
processes_file=$3
me=`hostname`
echo $me

set -e
for ((task=0;task<ndrivers;task++))
do
    $PHES_BINDIR/search_driver $tasks_file $processes_file &
done

echo Started $ndrivers searchers on $me

    