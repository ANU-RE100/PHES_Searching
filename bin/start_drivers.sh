#!/bin/bash

# file looks like
# hostname1 nproc
# hostname2 nproc
# ..

export PHES_BINDIR=bin

process_count_file=$1
tasks_file=$2

#myhostid=$1
#ndrivers=$1
me=`hostname`
echo $me

ndrivers=`grep $me $process_count_file | cut -d' ' -f 2`
task_start=`awk -v mh=$me '{ if (done==0) { if (index($1, mh)) done=1; else sum+=$2}}END{print sum}' $process_count_file`
total_searchers=`awk '{sum+=$2}END{print sum}' $process_count_file`

echo total_searchers=$total_searchers

# should cleanup search management directories (lock files dir, task done dir etc) but 

for ((task=0;task<ndrivers;task++))
do
    taskid=$((task + task_start))
    $PHES_BINDIR/search_driver $taskid $total_searchers $tasks_file &  #other args 
done

echo Started $ndrivers searchers on $me

    