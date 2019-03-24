
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -g
LIBS = -lgdal -lshp -lboost_system -lboost_filesystem
GEN_OBJS = src/reservoir.o src/coordinates.o src/phes_base.o  src/variable_parser.o src/kml.o src/csv.o
OBJS0 = src/shapefile_tiling.o $(GEN_OBJS)
OBJS1 = src/screening.o $(GEN_OBJS)
OBJS2 = src/pairing.o $(GEN_OBJS)
OBJS3 = src/pretty_set.o $(GEN_OBJS)
OBJS4 = src/constructor.o $(GEN_OBJS)
OBJS5 = src/search_driver.o $(GEN_OBJS)
DIRS = bin input output processing_files driver_files
INCDIRS = -Iinclude

utils: $(shell mkdir -p $(DIRS)) bin/shapefile_tiling bin/screening bin/pairing bin/pretty_set bin/constructor bin/search_driver

bin/shapefile_tiling: $(OBJS0)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJS0) $(LIBS) -o $@

bin/screening: $(OBJS1)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJS1) $(LIBS) -o $@

bin/pairing: $(OBJS2)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJS2) $(LIBS) -o $@

bin/pretty_set: $(OBJS3)
	g++  $(CXXFLAGS) $(LDFLAGS) $(OBJS3) $(LIBS) -o $@

bin/constructor: $(OBJS4)
	g++  $(CXXFLAGS) $(LDFLAGS) $(OBJS4) $(LIBS) -o $@

bin/search_driver: $(OBJS5)
	g++  $(CXXFLAGS) $(LDFLAGS) $(OBJS5) $(LIBS) -o $@

.c.o:
	g++ $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm -f $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) bin/screening bin/pairing bin/pretty_set bin/constructor bin/search_driver

clear:
	rm -r -f processing_files driver_files && mkdir -p $(DIRS)

run:
	bin/start_drivers.sh $(n)

install:
	sudo apt-get install libgdal-dev && sudo apt-get install libshp-dev && sudo apt-get install libboost-all-dev && sudo apt-get install gdal-bin