CXXFLAGS = -std=c++2a -O3 -Wall -Wextra -g
LIBS = -lgdal -lshp -lboost_system -lboost_filesystem
GEN_OBJS = src/reservoir.o src/coordinates.o src/phes_base.o  src/variable_parser.o src/kml.o src/csv.o src/polygons.o src/constructor_helpers.o
OBJS0 = src/shapefile_tiling.o $(GEN_OBJS)
OBJS1 = src/screening.o $(GEN_OBJS)
OBJS2 = src/pairing.o $(GEN_OBJS)
OBJS3 = src/pretty_set.o $(GEN_OBJS)
OBJS4 = src/constructor.o $(GEN_OBJS)
OBJS5 = src/search_driver.o $(GEN_OBJS)
DIRS = bin input output processing_files driver_files
INCDIRS = -Iinclude

utils: $(shell mkdir -p $(DIRS)) bin/shapefile_tiling bin/screening bin/pairing bin/pretty_set bin/constructor bin/search_driver bin/reservoir_constructor bin/depression_volume_finding

bin/shapefile_tiling: $(OBJS0)
	g++ $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

bin/screening: $(OBJS1)
	g++ $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

bin/pairing: $(OBJS2)
	g++ $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

bin/pretty_set: $(OBJS3)
	g++  $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

bin/constructor: $(OBJS4)
	g++  $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

bin/search_driver: $(OBJS5)
	g++  $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

bin/reservoir_constructor: src/reservoir_constructor.o $(GEN_OBJS)
	g++  $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

bin/depression_volume_finding: src/depression_volume_finding.o $(GEN_OBJS)
	g++  $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

.c.o:
	g++ $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm -f src/*.o bin/screening bin/pairing bin/pretty_set bin/constructor bin/search_driver bin/reservoir_constructor bin/depression_volume_finding bin/shapefile_tiling

clear:
	rm -r -f processing_files driver_files && mkdir -p $(DIRS)

run:
	bin/start_drivers.sh $(n)

compile_commands:
	compiledb make

install:
	sudo apt-get install libgdal-dev libshp-dev libboost-all-dev gdal-bin -y
