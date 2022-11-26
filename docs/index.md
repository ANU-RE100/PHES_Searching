# How to use PHES searching software

[Input](./input.md)

## Variables 

All variables that can be adjusted are in the variables file in the home directory of <storage location> and are briefly outlined in the table below. Each line is in the format <variable> = <value>; and // is treated as a comment (Similar to C++ formatting).

<!-- TODO add table -->

## Generating Shapefile Tiles 
Ensure filters to tile are set as variables using filter_to_tile = <filename>; in the variables file at <storage location>. Running shapefile_tiling (use ./bin/shapefile_tiling in bash) will then generate shapefile tiles (Subsets of the polygons) for each of the cells in the task list.

## Single Run
To run a single process, use either ./bin/<process> <lon> <lat> for greenfield searching of a DEM square or ./bin/<process> “<existing reservoir name>” for a brownfield site. <lon> and <lat> should be specified as integers between -180 to 180 and -90 to 90 respectively. DEMs and existing reservoirs must be run in the order screening -> pairing -> pretty_set -> constructor.

<!-- TODO update for cmake -->

## Bulk Run
__NOTE: Before a bulk run, delete <storage location>/driver_files and <storage location>/debug. If in <storage location> directory, can use make clean command instead (If make file present).__

To start a bulk run, ensure variables and tasks/processes files (see below) are set correctly, then use ./bin/search_driver to start one driver or make run n=<num processes> to start several drivers in parallel. The drivers will automatically go through the tasks and processes, storing any error information (See debug below) and will output done when all tasks are complete.

### Tasks File 
The tasks file (filename specified in variables) contains a list of all the tasks to complete (DEM squares or existing reservoir names). This is a text file with one task per line in the format <lon> <lat> for a DEM square or “<reservoir name>” for an existing reservoir. For example:

<!-- TODO add image -->

### Processes File

The processes file (filename specified in variables) contains a list of all the processes to complete (Out of screening, pairing, pretty_set and constructor). This is a text file with one process per line. For example:

<!-- TODO add image -->
