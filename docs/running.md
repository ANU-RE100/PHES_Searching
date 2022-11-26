# Compiling
Building this code depends on shapelib, gdal and boost and uses a CMake build system.

On Ubuntu, these dependencies can be installed using:

```bash
sudo apt-get install -y libgdal-dev libshp-dev libboost-all-dev gdal-bin cmake
```

Then, compile using CMake. For example

```
mkdir build
cd build
cmake ..
make -j $nproc
make install -j $nproc
cd ..
```

# Single Run
## Greenfield
To run a single process on a DEM square, use 
```
./bin/<process> <lon> <lat> 
```
where `<lon>` and `<lat>` should be specified as integers between -180 to 180 and -90 to 90 respectively. 

`<process>` refers to `screening`,  `pairing`, `pretty_set` or `constructor` which must be run in the order `screening` -> `pairing` -> `pretty_set` -> `constructor`.

## Brownfield
To run a single process on an existing reservoir use
```
./bin/<process> <existing reservoir name>
```
To run an existing reservoir, `screening` should be run with both the existing reservoir and the corresponding DEM square. `pairing`, `pretty_set` and `constructor` only need to be run with the existing reservoir name. 

## Bluefield
To run a single process with the ocean as a lower reservoir use
```
./bin/<process> ocean <lon> <lat>
```

## Debug mode
To output debugging statements, a `1` should follow the arguments provided in the above commands, for example
```
./bin/screening 144 -17 1
```
If no arguement is given or a number other than `1` is provided then only errors will be output. 

# Bulk Run
__NOTE: Before a bulk run, delete `<storage location>/driver_files` and `<storage location>/debug`. If in `<storage location> directory`, can use `make clean` command instead (if make file present).__

To start a bulk run, ensure variables and tasks/processes files (see below) are set correctly, then use `./bin/search_driver` to start one driver or `make run n=<num processes>` to start several drivers in parallel. The drivers will automatically go through the tasks and processes, storing any error information (see debug below) and will output done when all tasks are complete.

## Tasks File 
The tasks file (filename specified in variables) contains a list of all the tasks to complete (DEM squares or existing reservoir names). This is a text file with one task per line in the format `<lon>` `<lat>` for a DEM square or `<reservoir name>` for an existing reservoir. For example:

```
"Cooloolabin Dam"
152 -27
```

## Processes File

The processes file (filename specified in variables) contains a list of all the processes to complete (out of `screening`, `pairing`, `pretty_set` and `constructor`). This is a text file with one process per line. For example:

```
screening
pairing
pretty_set
constructor
```
