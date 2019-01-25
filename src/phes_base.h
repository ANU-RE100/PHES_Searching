/*
 * Some common parameters and base types for PHES utilities
 */

#ifndef PHES_BASE_H
#define PHES_BASE_H

//Being risky (for devel only)
#include <bits/stdc++.h>
using namespace std;

#include "model2D.h"

const double resolution = 30.87;         // Approx. 30 m for 1 arc-second DEM
const double EPS = 0.000001;             // Value added in flat sections
const double INF = 1000000000000000000;  
const int min_watershed_area = 10;      // Minimum watershed area in hectares
const int stream_threshold = (int)(11.1*min_watershed_area);        // Flow accumulation to be defined a creek
const double contour_height = 10.0;      // Contour interval for finding dam sites to test
const double freeboard = 1.5;            // Freeboard on dam

const int MIN_HEAD = 100;
const double min_reservoir_volume = 1.0;
const double min_reservoir_water_rock = 3.0;
const double min_max_dam_height = 5.0;

const int dambatter = 3;                      // Slope on sides of dam
const int cwidth = 10;                        // Width of top of dam
const int border = 600;                       // Number of cells to add as border around square

const double gravity = 9.8; 				  
const double generation_efficiency = 0.9;
const double usable_volume = 0.85;
const double J_GWh_conversion = 3.6e12;
const double water_density = 1000;
const double cubic_metres_GL_conversion = 1.0e6;

const int NWALL_HEIGHTS = 10;
const int MAX_WALL_HEIGHT = 100;
const double dam_wall_heights[NWALL_HEIGHTS] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}; //  Wall heights to test and export
const double all_wall_heights[NWALL_HEIGHTS+1] = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

struct Test{
	int energy_capacity;
	int storage_time;
	int slope_coefficient;
	int min_FOM;
	bool operator<(const Test &o) const
	    {
	    	if(energy_capacity==o.energy_capacity){
	    		return storage_time > o.storage_time;
	    	}
			return energy_capacity > o.energy_capacity;
	    }
};

const array<Test, 8> tests = {{
	{2, 6, 0, 0},
	{5, 18, 0, 0},
	{5, 6, 300, 0},
	{15, 18, 100, 0},
	{15, 6, 300, 0},
	{50, 18, 100, 0},
	{50, 6, 300, 0},
	{150, 18, 100, 0}}};

const int ndirections = 8;

struct direction{
	int row, col;
	int val;
};

const direction directions[ndirections] ={{0, 1, 1}, 
				      {1, 1, 2},
				      {1, 0, 4},
				      {1, -1, 8},
				      {0, -1, 16},
				      {-1, -1, 32},
				      {-1, 0, 64},
				      {-1, 1, 128}};

#define ARCSECLEN 30.87
#define DEGREELEN 111132.00



#ifndef MIN
#define MIN(x, y)  (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y)  (((x) > (y)) ? (x) : (y))
#endif
#ifndef MAX
#define AVG(x, y)  0.5*(x+y)
#endif

#define SQ(a) ((a)*(a))

#define RADIANS(x) (0.01745329251994329576*(x))
#define inv_rt2 0.70710678118654752440
#define inv3600 2.777777777777778e-4
#define COS cos
#define SQRT sqrt
#define FLOOR floor

struct GeographicCoordinate{
	double lat,lon;
};

struct ArrayCoordinate{
	int row, col;
	GeographicCoordinate origin;
};

struct GridSquare{
	int lat,lon;
};

struct ArrayCoordinateWithHeight {
	int row, col;
	double h;
	bool operator<(const ArrayCoordinateWithHeight &o) const
	    {
		return h > o.h;
	    }
};

// extreme points of shape in the 8 directions
struct Shape_bound {
	ArrayCoordinate bound[ndirections];
};


// In coordinates
GeographicCoordinate GeographicCoordinate_init(double latitude, double longitude);
ArrayCoordinate ArrayCoordinate_init(int row, int col);
ArrayCoordinate ArrayCoordinate_init(int row, int col, GeographicCoordinate origin);
GridSquare GridSquare_init(int latitude, int longitude);
ArrayCoordinateWithHeight ArrayCoordinateWithHeight_init(int row, int col, double h);

GeographicCoordinate get_origin(GridSquare square, int border);
bool check_within(ArrayCoordinateWithHeight c, int shape[2]);
bool check_within(ArrayCoordinate c, int shape[2]);
double find_slope(ArrayCoordinate c1, ArrayCoordinate c2, Model_double *DEM);
string str(GridSquare square);
int find_lowest_neighbor(ArrayCoordinate c, Model_double *DEM);
int find_lowest_neighbor(ArrayCoordinate c, Model_double *DEM, double coslat);

double find_area(ArrayCoordinate c);

double find_distance(ArrayCoordinate c1, ArrayCoordinate c2);
double find_distance(ArrayCoordinate c1, ArrayCoordinate c2, double coslat);
double find_distance_sqd(ArrayCoordinate c1, ArrayCoordinate c2);
double find_distance_sqd(ArrayCoordinate c1, ArrayCoordinate c2, double coslat);
double find_distance(GeographicCoordinate c1, GeographicCoordinate c2);
double find_distance(GeographicCoordinate c1, GeographicCoordinate c2, double coslat);
double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2);
double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2, double coslat);

int flows_to(ArrayCoordinate c1, ArrayCoordinate c2, Model_int16 *flow_directions);

GeographicCoordinate convert_coordinates(ArrayCoordinate c);
ArrayCoordinate convert_coordinates(GeographicCoordinate c, GeographicCoordinate origin);
double find_orthogonal_nn_distance(ArrayCoordinate c1, ArrayCoordinate c2);

struct RoughReservoir{
	string identifier;
	double latitude;
	double longitude;
	int elevation;
	ArrayCoordinate pour_point;
	double volumes[NWALL_HEIGHTS];
	double dam_volumes[NWALL_HEIGHTS];
	double areas[NWALL_HEIGHTS];
	double water_rocks[NWALL_HEIGHTS];
	double watershed_area;
	double max_dam_height;
	ArrayCoordinate shape_bound[NWALL_HEIGHTS][ndirections];
	bool operator<(const RoughReservoir &o) const
	    {
		return elevation > o.elevation;
	    }
};

struct Reservoir{
	string identifier;
	double latitude;
	double longitude;
	int elevation;
	ArrayCoordinate pour_point;
	double volume;
	double dam_volume;
	double dam_length;
	double area;
	double water_rock;
	double watershed_area;
	double average_water_depth;
	double dam_height;
	double max_dam_height;
	ArrayCoordinate shape_bound[ndirections];
	bool operator<(const Reservoir &o) const
	    {
		return elevation > o.elevation;
	    }
};

struct Pair {
	Reservoir upper;
	Reservoir lower;
	string identifier;
	double distance;
	double pp_distance;
	double slope;
	double required_volume;
	double volume;
	double FOM;
	double water_rock;
	int energy_capacity;
	int storage_time;
	int head;
	bool operator<(const Pair &o) const
	    {
		return FOM > o.FOM;
	    }
};

// In reservoir
RoughReservoir RoughReservoir_init(ArrayCoordinate pour_point, int elevation);
Reservoir Reservoir_init(ArrayCoordinate pour_point, int elevation);
void write_rough_reservoir_csv_header(FILE *csv_file);
void write_rough_reservoir_data_header(FILE *csv_file);
void write_rough_reservoir_csv(FILE *csv_file, RoughReservoir reservoir);
void write_rough_reservoir_data(FILE *csv_file, RoughReservoir reservoir);
vector<RoughReservoir> read_rough_reservoir_data(char* filename);
void write_rough_pair_csv_header(FILE *csv_file);
void write_rough_pair_data_header(FILE *csv_file);
void write_rough_pair_csv(FILE *csv_file, Pair *pair);
void write_rough_pair_data(FILE *csv_file, Pair *pair);
array<vector<Pair>,tests.size()> read_rough_pair_data(char* filename);
void update_reservoir_boundary(ArrayCoordinate (*dam_shape_bounds)[ndirections], ArrayCoordinate point, int elevation_above_pp);
void write_pair_csv_header(FILE *csv_file);
void write_pair_csv(FILE *csv_file, Pair *pair);

struct Models{
	GridSquare neighbors[9];
	Model_int16* DEMs[9];
	Model_int16* flow_directions[9];
	GeographicCoordinate origin;
};

// In phes_base
int convert_to_int(double f);
double max_over_wall_heights(double *a);
double convert_to_dam_volume(int height, double length);
double convert_to_dam_volume(int height, double length);
double linear_interpolate(double value, const double *x_values, const double *y_values);
string str(int i);
unsigned long walltime_usec();
double find_required_volume(int energy, int head);
char* convert_string(string str);
void write_to_csv_file(FILE *csv_file, vector<string> cols);
vector<string> read_from_csv_file(string line);
string dtos(double f, int nd);
Model_int16* read_DEM_with_borders(GridSquare sq);
Models Models_init(GridSquare sc);
void set_FOM(Pair* pair);
string str(Test test);


#endif