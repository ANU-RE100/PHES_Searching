#ifndef PHES_BASE_H
#define PHES_BASE_H

#include "gdal/gdal_priv.h"
#include "shapefil.h"
#include <gdal/cpl_conv.h>
#include <gdal/cpl_string.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <bits/stdc++.h>

#include "search_config.hpp"

using namespace std;

void parse_variables(char *filename);

const double EPS = 1.0e-6;
const double INF = 1.0e18;
const double J_GWh_conversion = 3.6e12;
const double cubic_metres_GL_conversion = 1.0e6;
const double resolution = 30.87;
const double meters_per_hectare = 10000;

// Search Driver
extern string tasks_file;     // File with list of cells to do line by line in
                              // format <lon> <lat>
extern string processes_file; // File with list of processes to complete
extern string existing_reservoirs_csv;
extern string existing_reservoirs_shp;
extern string existing_reservoirs_shp_names;

// General
extern string file_storage_location; // Where to look for input files and store
                                     // output files
extern int border;       // Number of cells to add as border around DEM square
extern double dambatter; // Slope on sides of dam
extern double cwidth;    // Width of top of dam
extern double freeboard; // Freeboard on dam

// Shapefile tiling
extern vector<string>
    filter_filenames_to_tile; // Shapefiles to split into tiles

// Screening
extern double min_watershed_area; // Minimum watershed area in hectares to be
                                  // consisered a stream
extern int stream_threshold;      // Number of cells required to reach minimum
                                  // watershed area
extern int contour_height; // Contour interval along streams for finding dam
                           // sites to test

extern double min_reservoir_volume; // Minimum reservoir volume (GL) at maximum
                                    // dam wall height
extern double min_reservoir_water_rock; // Minimum reservoir water to rock ratio
                                        // at optimal dam wall height
extern double
    min_max_dam_height; // Minimum maximum dam height (m) (Before overlapping
                        // filters) to be considered a potential reservoir

extern vector<string> filter_filenames;
extern vector<double> dam_wall_heights; //  Wall heights to test and export

// Pairing
extern int min_head; // Minimum head (m) to be considered a potential pair
extern int max_head; // Maximum head (m) to be considered a potential pair
extern double min_pair_water_rock; // Minimum pair water to rock ratio based on
                                   // interpolated values
extern double min_slope;    // Minimum slope based on interpolated nearest point
                            // seperation between two reservoirs
extern double min_pp_slope; // Minimum slope based on pourpoint seperation
                            // between two reservoirs
extern int max_lowers_per_upper; // Maximum number of lower reservoirs to keep
                                 // per upper reservoir
extern double tolerance_on_FOM;
extern double max_head_variability;   // Maximum amount the head can vary during
                                      // water transfer (Default 0.35)
extern int num_altitude_volume_pairs; // Number of altitude-volume pairs
                                      // provided with an existing pit
extern int pit_height_resolution; // Height resolution of top and bottom of pit
                                  // in metres

// Common
extern double gravity;               // Acceleration due to gravity (m/s/s)
extern double generation_efficiency; // Efficiency of generation
extern double usable_volume;         // Usable volume of reservoir
extern double water_density;         // Density of water (kg/m^3)
extern int max_wall_height;

// Output
extern int good_colour[4];
extern int bad_colour[4];
extern string upper_colour;
extern string lower_colour;
extern double volume_accuracy;            // Maximum ratio error on final volume
extern double dam_wall_height_resolution; // Resolution of dam wall height (m)
extern double minimum_dam_height;

// FOM Calculations
extern double powerhouse_coeff;
extern double power_exp;
extern double head_exp;
extern double power_slope_factor;
extern double slope_int;
extern double head_coeff;
extern double power_offset;
extern double tunnel_fixed;
extern double dam_cost;

// Ocean FOM Calculations
extern double lining_cost;
extern double sea_power_scaling;
extern double ref_marine_cost;
extern double ref_power;
extern double ref_head;

struct Test {
  double energy_capacity;
  int storage_time;
  bool operator<(const Test &o) const {
    if (abs(energy_capacity - o.energy_capacity) < EPS) {
      return storage_time > o.storage_time;
    }
    return energy_capacity > o.energy_capacity;
  }
};

extern vector<Test> tests;

struct Direction {
  int row, col;
};

const array<Direction, 8> directions = {
    {{0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}}};

struct CategoryCutoff {
  char category;
  double power_cost;
  double storage_cost;
  bool operator<(const CategoryCutoff &o) const {
    return category > o.category;
  }
};

extern vector<CategoryCutoff> category_cutoffs;

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef AVG
#define AVG(x, y) 0.5 * (x + y)
#endif
#ifndef SQ
#define SQ(a) ((a) * (a))
#endif

#define RADIANS(x) (0.01745329251994329576 * (x))
#define COS cos
#define SQRT sqrt
#define FLOOR floor

bool file_exists(char *name);
bool file_exists(string name);


#include "coordinates.h"

struct BigModel {
  GridSquare neighbors[9];
  Model<short> *DEM;
  Model<char> *flow_directions[9];
};



#include "model2D.h"
#include "polygons.h"
#include "reservoir.h"

#include "csv.h"

int convert_to_int(double f);
double max(vector<double> a);
double convert_to_dam_volume(int height, double length);
double convert_to_dam_volume(int height, double length);
double linear_interpolate(double value, vector<double> x_values,
                          vector<double> y_values);
string str(int i);
unsigned long walltime_usec();
double find_required_volume(int energy, int head);
char *convert_string(string str);
void write_to_csv_file(FILE *csv_file, vector<string> cols);
vector<string> read_from_csv_file(string line);
string dtos(double f, int nd);
Model<short> *read_DEM_with_borders(GridSquare sq, int border);
BigModel BigModel_init(GridSquare sc);
void set_FOM(Pair *pair);
string str(Test test);
string energy_capacity_to_string(double energy_capacity);
GeographicCoordinate get_origin(double latitude, double longitude, int border);
ExistingReservoir get_existing_reservoir(string name);
vector<ExistingReservoir> get_existing_reservoirs(GridSquare grid_square);
RoughBfieldReservoir existing_reservoir_to_rough_reservoir(ExistingReservoir r);
vector<ExistingPit> get_pit_details(GridSquare grid_square);
void depression_volume_finding(GridSquare grid_square);

#endif
