#include "phes_base.h"

// Search Driver
string tasks_file;					// File with list of cells to do line by line in format <lon> <lat>
string processes_file;				// File with list of processes to complete
string existing_reservoirs_csv;
string existing_reservoirs_shp;
string existing_reservoirs_shp_names;
bool use_tiled_bluefield;

// General
string file_storage_location;		// Where to look for input files and store output files
int border;							// Number of cells to add as border around DEM square
double dambatter;					// Slope on sides of dam
double cwidth;						// Width of top of dam
double freeboard;            		// Freeboard on dam
string dem_type;         // The digital elevation model to be used (SRTM or FABDEM)

// Shapefile tiling
vector<string> filter_filenames_to_tile; // Shapefiles to split into tiles

// Screening
double min_watershed_area;			// Minimum watershed area in hectares to be consisered a stream
int stream_threshold;				// Number of cells required to reach minimum watershed area
int contour_height;				// Contour interval along streams for finding dam sites to test

double min_reservoir_volume;		// Minimum reservoir volume (GL) at maximum dam wall height
double min_reservoir_water_rock;	// Minimum reservoir water to rock ratio at optimal dam wall height
double min_max_dam_height;			// Minimum maximum dam height (m) (Before overlapping filters) to be considered a potential reservoir

vector<string> filter_filenames;
vector<double> dam_wall_heights; 	//  Wall heights to test and export

// Pairing
int min_head;						// Minimum head (m) to be considered a potential pair
int max_head;						// Maximum head (m) to be considered a potential pair
double min_pair_water_rock;			// Minimum pair water to rock ratio based on interpolated values
double min_slope;					// Minimum slope based on interpolated nearest point seperation between two reservoirs
double min_pp_slope;				// Minimum slope based on pourpoint seperation between two reservoirs
int max_lowers_per_upper;			// Maximum number of lower reservoirs to keep per upper reservoir
double tolerance_on_FOM;
double max_head_variability;		// Maximum amount the head can vary during water transfer (Default 0.35)
int num_altitude_volume_pairs;		// Number of altitude-volume pairs provided with an existing pit
int pit_height_resolution;			// Height resolution of top and bottom of pit in metres

// Common
double gravity;						// Acceleration due to gravity (m/s/s)
double generation_efficiency;		// Efficiency of generation
double usable_volume;				// Usable volume of reservoir
double water_density;				// Density of water (kg/m^3)
int max_wall_height;

// Output
bool output_FOM;			// Whether to output exact FOM or category split
int good_colour[4];
int bad_colour[4];
string upper_colour;
string lower_colour;
double volume_accuracy;				// Maximum ratio error on final volume
double dam_wall_height_resolution;	// Resolution of dam wall height (m)
double minimum_dam_height;

// FOM Calculations
double powerhouse_coeff;
double power_exp;
double head_exp;
double power_slope_factor;
double slope_int;
double head_coeff;
double power_offset;
double tunnel_fixed;
double dam_cost;

// Ocean FOM Calculations
double lining_cost;
double sea_power_scaling;
double ref_marine_cost;
double ref_power;
double ref_head;

// Reservoir Sizings
vector<Test> tests;					// Test in format {Volume (GL), Storage time (h), Maximum FOM}
vector<CategoryCutoff> category_cutoffs;

double max_bluefield_surface_area_ratio;

void parse_variables(char* filename){
    if(!file_exists(filename)){
		search_config.logger.error("No file: " + string(filename));
		throw(1);
	}
	ifstream in(filename);
	string line;
	while(getline(in, line)){
		stringstream ss1(line);
		if(getline(ss1, line, ';')){
			string variable, value;
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			stringstream ss2(line);
			if((!getline(ss2, variable, '=') || !getline(ss2, value)) && line[0]!='/' && line.length()>=4){
				printf("Syntax error: %s\n", convert_string(line));
				continue;
			}
			if(variable=="filter")
				filter_filenames.push_back(value);
			if(variable=="filter_to_tile")
				filter_filenames_to_tile.push_back(value);
			if(variable=="dem_type")
				dem_type=value;
			if(variable=="min_watershed_area"){
				min_watershed_area = stod(value);
				stream_threshold = (int)(11.1*min_watershed_area);
			}
			if(variable=="border")
				border = stoi(value);
			if(variable=="contour_height")
				contour_height = stoi(value);
			if(variable=="min_head")
				min_head = stoi(value);
			if(variable=="min_reservoir_volume")
				min_reservoir_volume = stod(value);
			if(variable=="min_reservoir_water_rock")
				min_reservoir_water_rock = stod(value);
			if(variable=="min_max_dam_height")
				min_max_dam_height = stod(value);
			if(variable=="dambatter")
				dambatter = stod(value);
			if(variable=="cwidth")
				cwidth = stod(value);
			if(variable=="freeboard")
				freeboard = stod(value);
			if(variable=="gravity")
				gravity = stod(value);
			if(variable=="generation_efficiency")
				generation_efficiency = stod(value);
			if(variable=="usable_volume")
				usable_volume = stod(value);
			if(variable=="water_density")
				water_density = stod(value);
			if(variable=="dam_wall_heights"){
				vector<string> heights = read_from_csv_file(value);
				for(string height : heights){
					dam_wall_heights.push_back(stod(height));
				}
				sort(dam_wall_heights.begin(), dam_wall_heights.end());
				max_wall_height = dam_wall_heights[dam_wall_heights.size()-1];
			}
			if(variable=="test"){
				vector<string> t = read_from_csv_file(value);
				Test test = {stod(t[0]), stoi(t[1])};
				tests.push_back(test);
				sort(tests.begin(), tests.end());
			}
			if(variable=="powerhouse_coeff")
				powerhouse_coeff = stod(value);
			if(variable=="power_exp")
				power_exp = stod(value);
			if(variable=="head_exp")
				head_exp = stod(value);
			if(variable=="power_slope_factor")
				power_slope_factor = stod(value);
			if(variable=="slope_int")
				slope_int = stod(value);
			if(variable=="head_coeff")
				head_coeff = stod(value);
			if(variable=="power_offset")
				power_offset = stod(value);
			if(variable=="tunnel_fixed")
				tunnel_fixed = stod(value);
			if(variable=="dam_cost")
				dam_cost = stod(value);
			if(variable=="min_pair_water_rock")
				min_pair_water_rock = stod(value);
			if(variable=="min_slope")
				min_slope = stod(value);
			if(variable=="min_pp_slope")
				min_pp_slope = stod(value);
			if(variable=="tasks_file")
				tasks_file = value;
			if(variable=="processes_file")
				processes_file = value;
			if(variable=="good_colour"){
				vector<string> t = read_from_csv_file(value);
				good_colour[0] = stoi(t[0]);good_colour[1] = stoi(t[1]);good_colour[2] = stoi(t[2]);good_colour[3] = stoi(t[3]);
			}
			if(variable=="bad_colour"){
				vector<string> t = read_from_csv_file(value);
				bad_colour[0] = stoi(t[0]);bad_colour[1] = stoi(t[1]);bad_colour[2] = stoi(t[2]);bad_colour[3] = stoi(t[3]);
			}
			if(variable=="volume_accuracy")
				volume_accuracy = stod(value);
			if(variable=="dam_wall_height_resolution")
				dam_wall_height_resolution = stod(value);
			if(variable.length()==1){ // Category cutoffs
				vector<string> t = read_from_csv_file(value);
				CategoryCutoff cutoff = {variable[0], stod(t[0]), stod(t[1])};
				category_cutoffs.push_back(cutoff);
				sort(category_cutoffs.begin(), category_cutoffs.end());
			}
			if(variable=="upper_colour")
				upper_colour = value;
			if(variable=="lower_colour")
				lower_colour = value;
			if(variable=="file_storage_location")
				file_storage_location = value;
			if(variable=="max_head")
				max_head = stod(value);
			if(variable=="max_lowers_per_upper")
				max_lowers_per_upper = stoi(value);
			if(variable=="minimum_dam_height")
				minimum_dam_height = stod(value);
			if(variable=="existing_reservoirs_csv")
				existing_reservoirs_csv = value;
			if(variable=="existing_reservoirs_shp")
				existing_reservoirs_shp = value;
			if(variable=="existing_reservoirs_shp_names")
				existing_reservoirs_shp_names = value;
			if(variable=="tolerance_on_FOM")
				tolerance_on_FOM = stod(value);

			if(variable=="lining_cost")
				lining_cost = stod(value);
			if(variable=="sea_power_scaling")
				sea_power_scaling = stod(value);
			if(variable=="ref_marine_cost")
				ref_marine_cost = stod(value);
			if(variable=="ref_power")
				ref_power = stod(value);
			if(variable=="ref_head")
				ref_head = stod(value);
			if(variable=="max_head_variability")
				max_head_variability = stod(value);
			if(variable=="num_altitude_volume_pairs")
				num_altitude_volume_pairs = stoi(value);
			if(variable=="pit_height_resolution")
				pit_height_resolution = stoi(value);
			if(variable=="max_bluefield_surface_area_ratio")
				max_bluefield_surface_area_ratio = stod(value);
			if(variable=="use_tiled_bluefield")
				use_tiled_bluefield = stoi(value);
		}
	}
}

