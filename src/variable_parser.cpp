#include "phes_base.h"

double resolution;         // Approx. 30 m for 1 arc-second DEM
double min_watershed_area;      // Minimum watershed area in hectares
int stream_threshold;
double contour_height;      // Contour interval for finding dam sites to test
double freeboard;            // Freeboard on dam
int MIN_HEAD;
double min_reservoir_volume;
double min_reservoir_water_rock;
double min_max_dam_height;
double dambatter;                      // Slope on sides of dam
double cwidth;                        // Width of top of dam
int border;                       // Number of cells to add as border around square
double gravity; 				  
double generation_efficiency;
double usable_volume;
double J_GWh_conversion;
double water_density;
double cubic_metres_GL_conversion;
int MAX_WALL_HEIGHT;
vector<double> dam_wall_heights; //  Wall heights to test and export
vector<string> filter_filenames;
double powerhouse_coeff;
double power_exp;
double head_exp;
double power_slope_factor;
double slope_int;
double head_coeff;
double power_offset;
double tunnel_fixed;
double dam_cost;

vector<Test> tests;

void parse_variables(char* filename){
	ifstream in(filename);
	string line;
	while(getline(in, line)){
		stringstream ss1(line);
		if(getline(ss1, line, ';')){
			string variable, value;
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			stringstream ss2(line);
			if(!getline(ss2, variable, '=') || !getline(ss2, value)){
				printf("Syntax error: %s\n", convert_string(line));
				continue;
			}
			if(variable=="filter")
				filter_filenames.push_back(value);
			if(variable=="resolution")
				resolution = stod(value);
			if(variable=="min_watershed_area"){
				min_watershed_area = stod(value);
				stream_threshold = (int)(11.1*min_watershed_area);
			}
			if(variable=="border")
				border = stoi(value);
			if(variable=="contour_height")
				contour_height = stod(value);
			if(variable=="min_head")
				MIN_HEAD = stoi(value);
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
			if(variable=="J_GWh_conversion")
				J_GWh_conversion = stod(value);
			if(variable=="water_density")
				water_density = stod(value);
			if(variable=="cubic_metres_GL_conversion")
				cubic_metres_GL_conversion = stod(value);
			if(variable=="dam_wall_heights"){
				vector<string> heights = read_from_csv_file(value);
				for(string height : heights){
					dam_wall_heights.push_back(stod(height));
				}
				sort(dam_wall_heights.begin(), dam_wall_heights.end());
				MAX_WALL_HEIGHT = dam_wall_heights[dam_wall_heights.size()-1];
			}
			if(variable=="test"){
				vector<string> t = read_from_csv_file(value);
				Test test = {stoi(t[0]), stoi(t[1]), stoi(t[2])};
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
		}
	}
}
