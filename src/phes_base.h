#ifndef PHES_BASE_H
#define PHES_BASE_H

#include <sys/time.h>
#include <sys/stat.h> 
#include "shapefil.h"
#include <gdal/gdal.h>
#include <gdal/cpl_conv.h>
#include <gdal/cpl_string.h>

#include <bits/stdc++.h>
using namespace std;

void parse_variables(char* filename);

const double EPS = 1.0e-6;
const double INF = 1.0e18;  

extern double resolution;
extern double min_watershed_area;
extern int stream_threshold;
extern double contour_height;
extern double freeboard;
extern int MIN_HEAD;
extern double min_reservoir_volume;
extern double min_reservoir_water_rock;
extern double min_max_dam_height;
extern double dambatter;
extern double cwidth;
extern int border;
extern double gravity; 				  
extern double generation_efficiency;
extern double usable_volume;
extern double J_GWh_conversion;
extern double water_density;
extern double cubic_metres_GL_conversion;
extern int MAX_WALL_HEIGHT;
extern vector<double> dam_wall_heights;
extern vector<string> filter_filenames;
extern double powerhouse_coeff;
extern double power_exp;
extern double head_exp;
extern double power_slope_factor;
extern double slope_int;
extern double head_coeff;
extern double power_offset;
extern double tunnel_fixed;
extern double dam_cost;
extern int display;

struct Test{
	int energy_capacity;
	int storage_time;
	int min_FOM;
	bool operator<(const Test &o) const
	    {
	    	if(energy_capacity==o.energy_capacity){
	    		return storage_time > o.storage_time;
	    	}
			return energy_capacity > o.energy_capacity;
	    }
};

extern vector<Test> tests;

struct Direction{
	int row, col;
	int val;
};

const array<Direction, 8> directions =
	{{{ 0,  1, 1  }, 
      { 1,  1, 2  },
      { 1,  0, 4  },
      { 1, -1, 8  },
      { 0, -1, 16 },
      {-1, -1, 32 },
      {-1,  0, 64 },
      {-1,  1, 128}}};


#ifndef MIN
#define MIN(x, y)  (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y)  (((x) > (y)) ? (x) : (y))
#endif
#ifndef AVG
#define AVG(x, y)  0.5*(x+y)
#endif
#ifndef SQ
#define SQ(a) ((a)*(a))
#endif

#define RADIANS(x) (0.01745329251994329576*(x))
#define COS cos
#define SQRT sqrt
#define FLOOR floor

bool file_exists(char* name);

#include "model2D.h"
#include "TIFF_IO.h"
#include "coordinates.h"
#include "reservoir.h"
#include "csv.h"

struct Models{
	GridSquare neighbors[9];
	Model_int16* DEMs[9];
	Model_int16* flow_directions[9];
	GeographicCoordinate origin;
};

int convert_to_int(double f);
double max(vector<double> a);
double convert_to_dam_volume(int height, double length);
double convert_to_dam_volume(int height, double length);
double linear_interpolate(double value, vector<double> x_values, vector<double> y_values);
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
