# Single Run
To run a single process, use either ./bin/<process> <lon> <lat> for greenfield searching of a DEM square or ./bin/<process> “<existing reservoir name>” for a brownfield site. <lon> and <lat> should be specified as integers between -180 to 180 and -90 to 90 respectively. DEMs and existing reservoirs must be run in the order screening -> pairing -> pretty_set -> constructor.

### To run existing reservoir
Run pairing both with 
Bin/pairing “RESERVOIR NAME” 1
Bin/pairing <long> <lat> 1

<!-- TODO move to single run -->

<!-- TODO update for cmake -->

<!-- TODO add debugging instructions -->

# Bulk Run
__NOTE: Before a bulk run, delete <storage location>/driver_files and <storage location>/debug. If in <storage location> directory, can use make clean command instead (If make file present).__

To start a bulk run, ensure variables and tasks/processes files (see below) are set correctly, then use ./bin/search_driver to start one driver or make run n=<num processes> to start several drivers in parallel. The drivers will automatically go through the tasks and processes, storing any error information (See debug below) and will output done when all tasks are complete.

## Tasks File 
The tasks file (filename specified in variables) contains a list of all the tasks to complete (DEM squares or existing reservoir names). This is a text file with one task per line in the format <lon> <lat> for a DEM square or “<reservoir name>” for an existing reservoir. For example:

<!-- TODO add image -->

## Processes File

The processes file (filename specified in variables) contains a list of all the processes to complete (Out of screening, pairing, pretty_set and constructor). This is a text file with one process per line. For example:

<!-- TODO add image -->
