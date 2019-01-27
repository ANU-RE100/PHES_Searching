

//Being risky
#include <bits/stdc++.h>
using namespace std;

#include "phes_base.h"
#include "kml.h"

void update_reservoir_boundary(ArrayCoordinate (*dam_shape_bounds)[ndirections], ArrayCoordinate point, int elevation_above_pp)
{
	for (int ih=0; ih<NWALL_HEIGHTS; ih++) {
		int dam_height =  dam_wall_heights[ih];
		if (dam_height >= elevation_above_pp)
			for (int i = 0; i<ndirections; i++) {
				if ( (directions[i].row*point.row + directions[i].col * point.col) >
				     (directions[i].row*dam_shape_bounds[ih][i].row + directions[i].col*dam_shape_bounds[ih][i].col) ){
					dam_shape_bounds[ih][i] = point;
				}
			}
	}
}

RoughReservoir RoughReservoir_init(ArrayCoordinate pour_point, int elevation)
{
	RoughReservoir reservoir;
	reservoir.elevation = elevation;
	GeographicCoordinate geo_coordinate = convert_coordinates(pour_point);
	reservoir.latitude = geo_coordinate.lat;
	reservoir.longitude = geo_coordinate.lon;
	reservoir.pour_point = pour_point;
	reservoir.max_dam_height = MAX_WALL_HEIGHT;

	// initialize bounds
	for (int ih=0; ih < NWALL_HEIGHTS; ih++) {
		for (int idir=0; idir < ndirections; idir++) {
			reservoir.shape_bound[ih][idir].row = -100000 * directions[idir].row;
			reservoir.shape_bound[ih][idir].col = -100000 * directions[idir].col;
		}
	}

	return reservoir;
}

Reservoir Reservoir_init(ArrayCoordinate pour_point, int elevation)
{
	Reservoir reservoir;
	reservoir.elevation = elevation;
	GeographicCoordinate geo_coordinate = convert_coordinates(pour_point);
	reservoir.latitude = geo_coordinate.lat;
	reservoir.longitude = geo_coordinate.lon;
	reservoir.pour_point = pour_point;
	for (int idir=0; idir < ndirections; idir++) {
		reservoir.shape_bound[idir] = pour_point;
	}
	return reservoir;
}


void write_rough_reservoir_csv_header(FILE *csv_file)
{
	vector<string> header = { "Identifier", "Latitude", "Longitude", "Elevation (m)", "Maximum dam height", "Watershed area (ha)"};
	for (int i=0; i<NWALL_HEIGHTS;i++) header.push_back(dtos(dam_wall_heights[i],0)+"m reservoir volume (GL)");
	for (int i=0; i<NWALL_HEIGHTS;i++) header.push_back(dtos(dam_wall_heights[i],0)+"m reservoir area (ha)");
	for (int i=0; i<NWALL_HEIGHTS;i++) header.push_back(dtos(dam_wall_heights[i],0)+"m dam volume (GL)");
	for (int i=0; i<NWALL_HEIGHTS;i++) header.push_back(dtos(dam_wall_heights[i],0)+"m water to rock");
	write_to_csv_file(csv_file, header);
}

void write_rough_reservoir_data_header(FILE *csv_file)
{
	vector<string> header = { "Identifier", "Latitude", "Longitude", "Elevation (m)", "Maximum dam height", "Watershed area (ha)" };
	for (int i=0; i<NWALL_HEIGHTS;i++) header.push_back(dtos(dam_wall_heights[i],0)+"m reservoir volume (GL)");
	for (int i=0; i<NWALL_HEIGHTS;i++) header.push_back(dtos(dam_wall_heights[i],0)+"m dam volume (GL)");
	header.push_back("Boundary");
	write_to_csv_file(csv_file, header);
}

void write_rough_reservoir_csv(FILE *csv_file, RoughReservoir reservoir)
{
	vector<string> line = {reservoir.identifier, dtos(reservoir.latitude, 4), dtos(reservoir.longitude,4), to_string(reservoir.elevation), dtos(reservoir.max_dam_height,0), dtos(reservoir.watershed_area,1)};
	for (int i=0; i<NWALL_HEIGHTS;i++) line.push_back(dtos(reservoir.volumes[i],2));
	for (int i=0; i<NWALL_HEIGHTS;i++) line.push_back(dtos(reservoir.areas[i],2));
	for (int i=0; i<NWALL_HEIGHTS;i++) line.push_back(dtos(reservoir.dam_volumes[i],2));
	for (int i=0; i<NWALL_HEIGHTS;i++) line.push_back(dtos(reservoir.water_rocks[i],1));
	write_to_csv_file(csv_file, line);
}

void write_rough_reservoir_data(FILE *csv_file, RoughReservoir reservoir)
{
	vector<string> line = {reservoir.identifier, dtos(reservoir.latitude, 5), dtos(reservoir.longitude,5), to_string(reservoir.elevation), dtos(reservoir.max_dam_height,1), dtos(reservoir.watershed_area,2)};
	for (int i=0; i<NWALL_HEIGHTS;i++) line.push_back(dtos(reservoir.volumes[i],5));
	for (int i=0; i<NWALL_HEIGHTS;i++) line.push_back(dtos(reservoir.dam_volumes[i],5));
	for (int ih=0; ih< NWALL_HEIGHTS; ih++) {
		for (int idir=0;idir < ndirections; idir++){
			line.push_back(to_string(reservoir.shape_bound[ih][idir].row));
			line.push_back(to_string(reservoir.shape_bound[ih][idir].col));
		}
	}
	write_to_csv_file(csv_file, line);
}

vector<RoughReservoir> read_rough_reservoir_data(char* filename)
{
	vector<RoughReservoir> reservoirs;
	ifstream inputFile(filename);
	string s;
	bool header = true;
    while (getline(inputFile, s)) {

    	if(header){
    		header = false;
    		continue;
    	}
        vector<string> line = read_from_csv_file(s);
        GeographicCoordinate gc = GeographicCoordinate_init(stod(line[1]),stod(line[2]));
        GeographicCoordinate origin = get_origin(GridSquare_init((int)FLOOR(gc.lat)-EPS,(int)FLOOR(gc.lon)+EPS), border);
        RoughReservoir reservoir = RoughReservoir_init(convert_coordinates(gc, origin),stoi(line[3]));
        for (int i=0; i<NWALL_HEIGHTS;i++) reservoir.volumes[i] = stod(line[6+i]);
		for (int i=0; i<NWALL_HEIGHTS;i++) reservoir.dam_volumes[i] = stod(line[6+NWALL_HEIGHTS+i]);
		reservoir.max_dam_height = stod(line[4]);
		reservoir.watershed_area = stod(line[5]);
		reservoir.identifier = line[0];
		for (int ih=0; ih< NWALL_HEIGHTS; ih++) {
			for (int idir=0;idir < ndirections; idir++){
				//printf("%d %d:\n", ih, idir);
				//printf("%d\n", stoi(line[6+2*NWALL_HEIGHTS+(ih*NWALL_HEIGHTS+idir)*2]));
				//printf("%d\n", 6+2*NWALL_HEIGHTS+(ih*NWALL_HEIGHTS+idir)*2);
				reservoir.shape_bound[ih][idir].row = stoi(line[6+2*NWALL_HEIGHTS+(ih*ndirections+idir)*2]);
				//printf("3\n");
				reservoir.shape_bound[ih][idir].col = stoi(line[6+2*NWALL_HEIGHTS+1+(ih*ndirections+idir)*2]);
				reservoir.shape_bound[ih][idir].origin = origin;
			}
		}
        reservoirs.push_back(reservoir);
    }
    if(header){
    	throw 1;
    }

	return reservoirs;
}

void write_rough_pair_csv_header(FILE *csv_file)
{
	vector<string> header = {"Pair Identifier",
	"Upper Identifier","Upper latitude","Upper longitude","Upper elevation (m)","Upper dam height (m)","Upper max dam height (m)","Upper water to rock estimate",
	"Lower Identifier","Lower latitude","Lower longitude","Lower elevation (m)","Lower dam height (m)","Lower max dam height (m)","Lower water to rock estimate",
	"Head (m)","Pourpoint distance (km)","Distance (km)","Slope","Volume (GL)","Energy (GWh)","Storage time (h)","Figure of merit"};
	write_to_csv_file(csv_file, header);
}

void write_rough_pair_data_header(FILE *csv_file)
{
	vector<string> header = {"Pair Identifier",
	"Upper Identifier","Upper latitude","Upper longitude","Upper elevation (m)","Upper dam height (m)","Upper max dam height (m)","Upper water to rock estimate",
	"Lower Identifier","Lower latitude","Lower longitude","Lower elevation (m)","Lower dam height (m)","Lower max dam height (m)","Lower water to rock estimate",
	"Head (m)","Pourpoint distance (km)","Distance (km)","Slope","Volume (GL)","Energy (GWh)","Storage time (h)","Figure of merit"};
	write_to_csv_file(csv_file, header);
}

void write_rough_pair_csv(FILE *csv_file, Pair *pair)
{
	vector<string> line = {pair->identifier, 
	pair->upper.identifier, dtos(pair->upper.latitude,4), dtos(pair->upper.longitude,4), to_string(pair->upper.elevation), dtos(pair->upper.dam_height,1), dtos(pair->upper.max_dam_height,0), dtos(pair->upper.water_rock,1),
	pair->lower.identifier, dtos(pair->lower.latitude,4), dtos(pair->lower.longitude,4), to_string(pair->lower.elevation), dtos(pair->lower.dam_height,1), dtos(pair->lower.max_dam_height,0), dtos(pair->lower.water_rock,1),
	to_string(pair->head), dtos(pair->pp_distance, 2), dtos(pair->distance, 2), dtos(pair->slope, 2), dtos(pair->required_volume, 2), to_string(pair->energy_capacity), to_string(pair->storage_time), dtos(pair->FOM,1)};
	write_to_csv_file(csv_file, line);
}

void write_rough_pair_data(FILE *csv_file, Pair *pair)
{
	vector<string> line = {pair->identifier, 
	pair->upper.identifier, dtos(pair->upper.latitude,6), dtos(pair->upper.longitude,6), to_string(pair->upper.elevation), dtos(pair->upper.dam_height,3), dtos(pair->upper.max_dam_height,1), dtos(pair->upper.water_rock,5),
	pair->lower.identifier, dtos(pair->lower.latitude,6), dtos(pair->lower.longitude,6), to_string(pair->lower.elevation), dtos(pair->lower.dam_height,3), dtos(pair->lower.max_dam_height,1), dtos(pair->lower.water_rock,5),
	to_string(pair->head), dtos(pair->pp_distance, 5), dtos(pair->distance, 5), dtos(pair->slope, 6), dtos(pair->required_volume, 5), to_string(pair->energy_capacity), to_string(pair->storage_time), dtos(pair->FOM,3)};
	write_to_csv_file(csv_file, line);
}

array<vector<Pair>,tests.size()> read_rough_pair_data(char* filename)
{
	array<vector<Pair>,tests.size()> pairs;

	ifstream inputFile(filename);
	string s;
	bool header = true;
    while (getline(inputFile, s)) {
    	if(header){
    		header = false;
    		continue;
    	}
    	vector<string> line = read_from_csv_file(s);

    	Pair pair;

    	GeographicCoordinate gc = GeographicCoordinate_init(stod(line[2]), stod(line[3]));
    	GeographicCoordinate origin = get_origin(GridSquare_init((int)(FLOOR(gc.lat)-EPS),(int)(FLOOR(gc.lon)+EPS)), border);
    	pair.upper = Reservoir_init(convert_coordinates(gc, origin), stoi(line[4]));
    	gc = GeographicCoordinate_init(stod(line[9]), stod(line[10]));
    	origin = get_origin(GridSquare_init((int)FLOOR(gc.lat)-EPS,(int)FLOOR(gc.lon)+EPS), border);
    	pair.lower = Reservoir_init(convert_coordinates(gc, origin), stoi(line[11]));
    	//printf("%f %f %d %f\n", gc.lat, origin.lat, convert_coordinates(gc, origin).row, pair.lower.latitude);

    	pair.identifier = line[0]; 

		pair.upper.identifier = line[1];
		pair.upper.dam_height = stod(line[5]);
		pair.upper.max_dam_height = stod(line[6]);
		pair.upper.water_rock = stod(line[7]);

		pair.lower.identifier = line[8];
		pair.lower.dam_height = stod(line[12]);
		pair.lower.max_dam_height = stod(line[13]);
		pair.lower.water_rock = stod(line[14]);

		pair.head = stoi(line[15]);
		pair.pp_distance = stod(line[16]);
		pair.distance = stod(line[17]);
		pair.slope = stod(line[18]);
		pair.required_volume = stod(line[19]);
		pair.upper.volume = stod(line[19]);
		pair.lower.volume = stod(line[19]);
		pair.energy_capacity = stoi(line[20]);
		pair.storage_time = stoi(line[21]);
		pair.FOM = stod(line[22]);

		for(uint i = 0; i<tests.size(); i++)
			if(pair.energy_capacity == tests[i].energy_capacity && pair.storage_time == tests[i].storage_time)
				pairs[i].push_back(pair);
    }

    if(header)
    	throw 1;

    return pairs;
}

void write_pair_csv_header(FILE *csv_file)
{
	vector<string> header = {"Pair Identifier","Figure of Merit","Head (m)","Distance (km)","Slope (%)","Volume (GL)","Energy (GWh)","Storage time (h)","Combined water to rock ratio",
	"Upper Identifier","Upper elevation (m)","Upper latitude","Upper longitude","Upper reservoir area (ha)","Upper reservoir volume (GL)","Upper dam height (m)","Upper dam length (m)","Upper dam volume (GL)","Upper water to rock ratio",
	"Lower Identifier","Lower elevation (m)","Lower latitude","Lower longitude","Lower reservoir area (ha)","Lower reservoir volume (GL)","Lower dam height (m)","Lower dam length (m)","Lower dam volume (GL)","Lower water to rock ratio"};
	write_to_csv_file(csv_file, header);
}

void write_pair_csv(FILE *csv_file, Pair *pair)
{
	vector<string> line = {pair->identifier, dtos(pair->FOM,0),to_string(pair->head), dtos(pair->distance, 2), dtos(pair->slope*100, 0),  dtos(pair->volume, 1), to_string(pair->energy_capacity), to_string(pair->storage_time), dtos(pair->water_rock, 1),
	pair->upper.identifier, to_string(pair->upper.elevation), dtos(pair->upper.latitude,4), dtos(pair->upper.longitude,4), dtos(pair->upper.area,0), dtos(pair->upper.volume,1),dtos(pair->upper.dam_height,1),dtos(pair->upper.dam_length,0), dtos(pair->upper.dam_volume,2), dtos(pair->upper.water_rock,1),
	pair->lower.identifier, to_string(pair->lower.elevation), dtos(pair->lower.latitude,4), dtos(pair->lower.longitude,4), dtos(pair->lower.area,0), dtos(pair->lower.volume,1),dtos(pair->lower.dam_height,1),dtos(pair->lower.dam_length,0), dtos(pair->lower.dam_volume,2), dtos(pair->lower.water_rock,1)};
	write_to_csv_file(csv_file, line);
}

void write_total_csv_header(FILE *csv_file)
{
	vector<string> header = {"Grid Identifier","Number of paired sites","Total potential capacity (GWh)"};
	write_to_csv_file(csv_file, header);
}

void write_total_csv(FILE *csv_file, string square_name, int num_sites, int energy_capacity)
{
	vector<string> line = {square_name, to_string(num_sites), to_string(energy_capacity)};
	write_to_csv_file(csv_file, line);
}
