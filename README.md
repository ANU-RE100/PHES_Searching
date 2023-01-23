# ANU RE100 PHES Searching Software

## Installation

Building this code depends on shapelib, gdal and boost and uses a CMake build system.

On Ubuntu, these dependencies can be installed using:

```bash
sudo apt-get install -y libgdal-dev libshp-dev libboost-all-dev gdal-bin cmake g++
```

Then, compile using CMake. For example

```bash
mkdir build
cd build
cmake ..
make -j
make install -j
cd ..
```

or a convenient script has been set up to perform this with the command

```bash
./build.sh
```

Alternately, there is a docker file for use on other systems, this can be run with

```
docker build . -f docker/Dockerfile -t phes
docker run -v {path_to_storage_location} -it --rm phes
```

## Usage

Place all input files in `PHES_Searching/input/`

Edit variables in `PHES_Searching/variables`

To run:

```bash
   ./bin/screening <long> <lat> 1
   ./bin/pairing <long> <lat> 1
   ./bin/pretty_set <long> <lat> 1
   ./bin/constructor <long> <lat> 1
```

To run large batch:

```bash
	./bin/start_drivers.sh <num_processes>
```
