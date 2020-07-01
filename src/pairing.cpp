#include "phes_base.h"

int display = false;
bool ocean = false;
bool pit = false;
ExistingPit pit_details;

vector<int> pairs;

vector<GeographicCoordinate> find_points_to_test(vector<array<ArrayCoordinate, directions.size()>>& boundary, double& wall_height, ArrayCoordinate& pour_point)
{
	array<ArrayCoordinate, directions.size()> one_point = { pour_point, pour_point, pour_point, pour_point,
					  pour_point, pour_point, pour_point, pour_point };

	int i = 0;
	while (dam_wall_heights[i] < wall_height) {
		i+=1;
	}

	int lower_wall_height = (i) ? dam_wall_heights[i-1] : 0;
	array<ArrayCoordinate, directions.size()> lower_shape = (i) ? boundary[i-1] : one_point;
	vector<GeographicCoordinate> bound;

	double inv_wall_height_interval = 0.1;
	
	for (uint j=0; j<directions.size(); j++) {
		GeographicCoordinate point1 = convert_coordinates(lower_shape[j]);
		GeographicCoordinate point2 = convert_coordinates(boundary[i][j]);
		bound.push_back((GeographicCoordinate) { point1.lat+(point2.lat-point1.lat)*(wall_height-lower_wall_height)*inv_wall_height_interval,
				     point1.lon+(point2.lon-point1.lon)*(wall_height-lower_wall_height)*inv_wall_height_interval });
	}

	bound.push_back(convert_coordinates(pour_point));

	return bound;
}


double find_least_distance_sqd(vector<array<ArrayCoordinate, directions.size()>>& upper_boundary, vector<array<ArrayCoordinate, directions.size()>>& lower_boundary,
			      double upper_wall_height, double lower_wall_height,
			      ArrayCoordinate upper_pour_point, ArrayCoordinate lower_pour_point)
{
	double mindist2 = INF;
	vector<GeographicCoordinate> upper_points = find_points_to_test(upper_boundary, upper_wall_height, upper_pour_point);
	vector<GeographicCoordinate> lower_points = find_points_to_test(lower_boundary, lower_wall_height, lower_pour_point);

	for (uint iu = 0; iu<upper_points.size(); iu++) {
		GeographicCoordinate p1 = upper_points[iu];
		for (uint il=0; il<lower_points.size(); il++) {
			GeographicCoordinate p2 = lower_points[il];
			mindist2 = MIN(mindist2, find_distance_sqd(p1,p2));
		}
	}

	return mindist2;
}

int max_altitude(vector<AltitudeVolumePair> pairs){
	return pairs[pairs.size()-1].altitude;
}

vector<int> get_altitudes(ExistingPit& pit){
	vector<int> to_return;
	for(AltitudeVolumePair pair : pit.volumes)
		to_return.push_back(pair.altitude);
	return to_return;
}

vector<double> get_volumes(ExistingPit& pit){
	vector<double> to_return;
	for(AltitudeVolumePair pair : pit.volumes)
		to_return.push_back(pair.volume);
	return to_return;
}

vector<double> int_to_double_vector(vector<int> int_vector){
	vector<double> to_return(int_vector.begin(), int_vector.end());
	return to_return;
}

double pit_volume(ExistingPit& pit, int bottom_elevation, int top_elevation){
	return linear_interpolate(top_elevation, int_to_double_vector(get_altitudes(pit)), get_volumes(pit))-linear_interpolate(bottom_elevation, int_to_double_vector(get_altitudes(pit)), get_volumes(pit));
}

bool determine_pit_elevation_and_volume(RoughReservoir& upper, RoughReservoir& lower, double energy_capacity, ExistingPit& pit_details, double& required_volume, int& head){
	RoughReservoir& greenfield = upper;
	RoughReservoir& pit = lower;
	if(upper.brownfield){
		greenfield = lower;
		pit = upper;
	}
	while(pit.elevation < max_altitude(pit_details.volumes)){
		pit.max_dam_height=max_altitude(pit_details.volumes)-pit.elevation;
		int pit_depth = 0;
		while(pit_depth < pit.max_dam_height){
			pit_depth+=1;
			double volume = pit_volume(pit_details, pit.elevation, pit.elevation+pit_depth);
			double greenfield_wall_height = linear_interpolate(volume, greenfield.volumes, dam_wall_heights);
			head = (int) ABS(((0.5*(double) greenfield_wall_height+(double) greenfield.elevation)-(0.5*(double) pit_depth+(double) pit.elevation)));
			if ( head < min_head || head>max_head)
				continue;
			double head_ratio = (head+0.5*(greenfield_wall_height+(double) pit_depth))/(head-0.5*(greenfield_wall_height+(double) pit_depth));
			// cout << volume << " " << greenfield_wall_height << " " << greenfield.elevation << " " << pit_depth << " " << pit.elevation << " " << head << " " << head_ratio << "\n";
			if(head_ratio>(1+max_head_variability)){
				break;
			}
			if(volume<find_required_volume(energy_capacity, head)){
				continue;
			}
			required_volume = volume;
			return true;
		}
		
		pit.elevation += pit_height_resolution;
	}
	return false;
}

Pair *check_good_pair(RoughReservoir upper, RoughReservoir lower, double energy_capacity, int storage_time, Pair *pair, int max_FOM)
{
	int head = upper.elevation - lower.elevation;
	double required_volume = find_required_volume(energy_capacity, head);
	if ( (max(upper.volumes) < required_volume) ||
	     (max(lower.volumes) < required_volume) )
		return NULL;

	if(pit && !determine_pit_elevation_and_volume(upper, lower, energy_capacity, pit_details, required_volume, head))
		return NULL;

	double upper_dam_wall_height = 0;
	double lower_dam_wall_height = 0;
	double upper_water_rock_estimate = INF;
	double lower_water_rock_estimate = INF;

	if(!upper.brownfield){
		upper_dam_wall_height = linear_interpolate(required_volume, upper.volumes, dam_wall_heights);
		upper_water_rock_estimate = required_volume/linear_interpolate(upper_dam_wall_height, dam_wall_heights, upper.dam_volumes);
	}else{
		if(pit)
			upper_dam_wall_height = linear_interpolate(required_volume+pit_volume(pit_details, pit_details.reservoir.elevation, upper.elevation), get_volumes(pit_details), int_to_double_vector(get_altitudes(pit_details)))-upper.elevation;
		else
			upper_dam_wall_height = dam_wall_heights[0];
		upper_water_rock_estimate = INF;
	}
	if(!lower.brownfield && !lower.ocean){
		lower_dam_wall_height = linear_interpolate(required_volume, lower.volumes, dam_wall_heights);
		lower_water_rock_estimate = required_volume/linear_interpolate(lower_dam_wall_height, dam_wall_heights, lower.dam_volumes);
	}else{
		if(pit)
			lower_dam_wall_height = linear_interpolate(required_volume+pit_volume(pit_details, pit_details.reservoir.elevation, lower.elevation), get_volumes(pit_details), int_to_double_vector(get_altitudes(pit_details)))-lower.elevation;
		else
			lower_dam_wall_height = dam_wall_heights[0];
		lower_water_rock_estimate = INF;
	}

	if ( (!upper.brownfield && upper_dam_wall_height > upper.max_dam_height) || (!lower.brownfield && !lower.ocean && lower_dam_wall_height > lower.max_dam_height) )
		return NULL;

	if ( (upper_water_rock_estimate*lower_water_rock_estimate) < min_pair_water_rock*(upper_water_rock_estimate+lower_water_rock_estimate) )
		return NULL;

	ArrayCoordinate upper_coordinates = upper.pour_point;
	ArrayCoordinate lower_coordinates = lower.pour_point;

	double least_distance = find_least_distance_sqd(upper.shape_bound, lower.shape_bound,
							   upper_dam_wall_height, lower_dam_wall_height,
							   upper_coordinates, lower_coordinates);

	if ( SQ(head*0.001) < least_distance*SQ(min_slope))
		return NULL;

	Reservoir upper_reservoir = Reservoir_init(upper.pour_point, upper.elevation);
	Reservoir lower_reservoir = Reservoir_init(lower.pour_point, lower.elevation);

	upper_reservoir.identifier = upper.identifier;
	upper_reservoir.volume = required_volume;
	// if(upper.brownfield)
	// 	upper_reservoir.volume = upper.volumes[0];
	if(!upper.brownfield){
		upper_reservoir.dam_volume = linear_interpolate(upper_dam_wall_height, dam_wall_heights, upper.dam_volumes);
		upper_reservoir.area = linear_interpolate(upper_dam_wall_height, dam_wall_heights, upper.areas);
	}
	upper_reservoir.water_rock = upper_water_rock_estimate;
	upper_reservoir.dam_height = upper_dam_wall_height;
	upper_reservoir.max_dam_height = upper.max_dam_height;
	upper_reservoir.brownfield = upper.brownfield;
	upper_reservoir.pit = upper.pit;

	lower_reservoir.identifier = lower.identifier;
	lower_reservoir.volume = required_volume;
	// if(lower.brownfield)
	// 	lower_reservoir.volume = lower.volumes[0];
	if(!lower.brownfield){
		lower_reservoir.dam_volume = linear_interpolate(lower_dam_wall_height, dam_wall_heights, lower.dam_volumes);
		lower_reservoir.area = linear_interpolate(lower_dam_wall_height, dam_wall_heights, lower.areas);
	}
	lower_reservoir.water_rock = lower_water_rock_estimate;
	lower_reservoir.dam_height = lower_dam_wall_height;
	lower_reservoir.max_dam_height = lower.max_dam_height;
	lower_reservoir.brownfield = lower.brownfield;
	lower_reservoir.pit = lower.pit;
	lower_reservoir.ocean = lower.ocean;

	pair->identifier = upper.identifier+" & "+lower.identifier;
	pair->upper = upper_reservoir;
	pair->lower = lower_reservoir;
	pair->head = head;
	pair->distance = SQRT(least_distance);
	pair->pp_distance = find_distance(pair->upper.pour_point, pair->lower.pour_point);
	pair->energy_capacity = energy_capacity;
	pair->storage_time = storage_time;
	pair->required_volume = required_volume;
    pair->slope = pair->head/(pair->distance)*0.001;
    pair->water_rock = 1/(1/pair->upper.water_rock+1/pair->lower.water_rock);

	set_FOM(pair);

	if(pair->FOM>max_FOM)
		return NULL;
	return pair;
}

void pairing(vector<RoughReservoir>& upper_reservoirs, vector<RoughReservoir>& lower_reservoirs, FILE *csv_file, FILE *csv_data_file)
{
	vector<set<Pair>> temp_pairs;
	for (uint itest=0; itest<tests.size(); itest++) {
		pairs.push_back(0);
		set<Pair> a;
		temp_pairs.push_back(a);
	}

	RoughReservoir* upper_reservoir;
	RoughReservoir* lower_reservoir;
	for (uint iupper=0;iupper < upper_reservoirs.size(); iupper++) {
		upper_reservoir = &upper_reservoirs[iupper];
		double coslat = COS(upper_reservoir->latitude);
		for (uint ilower = 0; ilower<lower_reservoirs.size(); ilower++){
		    lower_reservoir = &lower_reservoirs[ilower];

		    int head = upper_reservoir->elevation - lower_reservoir->elevation;
			if ( head < min_head || head>max_head)
				continue;

			// Pour point separation
			double distance_sqd = find_distance_sqd(upper_reservoir->pour_point,lower_reservoir->pour_point, coslat);

			if (SQ(head*0.001) <= distance_sqd*SQ(min_pp_slope))
				continue;

			for (uint itest=0; itest<tests.size(); itest++) {
				Pair temp_pair;
				int max_FOM = (category_cutoffs[0].storage_cost*tests[itest].storage_time+category_cutoffs[0].power_cost)*(1+tolerance_on_FOM);
				if (check_good_pair(*upper_reservoir, *lower_reservoir, tests[itest].energy_capacity, tests[itest].storage_time, &temp_pair, max_FOM)) {
			 		temp_pairs[itest].insert(temp_pair);
			 		if((int)temp_pairs[itest].size()>max_lowers_per_upper || (ocean && temp_pairs[itest].size()>1))
			 			temp_pairs[itest].erase(prev(temp_pairs[itest].end())); 
				}
			}
		}
		for (uint itest=0; itest<tests.size(); itest++) {
			for(Pair pair : temp_pairs[itest]){
				write_rough_pair_csv(csv_file, &pair);
			 	write_rough_pair_data(csv_data_file, &pair);
			 	pairs[itest]++;
			}
			temp_pairs[itest].clear();
		}
	}
}


int main(int nargs, char **argv)
{

	GridSquare square_coordinate;
	bool brownfield = false;
	string fname;
	string prefix = "";
	string ocean_prefix = "";
	string arg1(argv[1]);

	int adj = 0;
	if(arg1.compare("ocean")==0){
		ocean = true;
		prefix = "ocean_";
		ocean_prefix = "ocean_";
		adj = 1;
		arg1 = argv[1+adj];
	}
	if(arg1.compare("pit")==0){
		brownfield = true;
		pit = true;
		prefix = "pit_";
		adj = 1;
		arg1 = argv[1+adj];
		fname = prefix+format_for_filename(arg1);
		if(nargs>2+adj)
			display = atoi(argv[2+adj]);
		printf("Pairing started for %s\n",argv[1+adj]);
	}else{
		try{
			int lon = stoi(arg1);
			square_coordinate = GridSquare_init(atoi(argv[2+adj]), lon);
			if(nargs>3+adj)
				display = atoi(argv[3+adj]);
			fname=prefix+str(square_coordinate);
			printf("Pairing started for %s\n",convert_string(fname));
		}catch(exception e){
			brownfield = true;
			fname = prefix+format_for_filename(arg1);
			if(nargs>2+adj)
				display = atoi(argv[2+adj]);
			printf("Pairing started for %s\n",argv[1+adj]);
		}
	}

	unsigned long t_usec = walltime_usec();
	parse_variables(convert_string("storage_location"));
	parse_variables(convert_string(file_storage_location+"variables"));

	vector<RoughReservoir> upper_reservoirs;
	if(brownfield){
		square_coordinate = get_square_coordinate(get_existing_reservoir(arg1));
		upper_reservoirs = read_rough_reservoir_data(convert_string(file_storage_location+"processing_files/reservoirs/"+fname+"_reservoirs_data.csv"));
		if (pit)
			pit_details = get_pit_details(arg1);
	} else
		upper_reservoirs = read_rough_reservoir_data(convert_string(file_storage_location+"processing_files/reservoirs/"+str(square_coordinate)+"_reservoirs_data.csv"));
	
	if(display)
		printf("Read in %zu uppers\n", upper_reservoirs.size());
	vector<RoughReservoir> lower_reservoirs;

	GridSquare neighbors[9] = {
		(GridSquare){square_coordinate.lat  ,square_coordinate.lon  },
		(GridSquare){square_coordinate.lat+1,square_coordinate.lon-1},
		(GridSquare){square_coordinate.lat+1,square_coordinate.lon  },
		(GridSquare){square_coordinate.lat+1,square_coordinate.lon+1},
		(GridSquare){square_coordinate.lat  ,square_coordinate.lon+1},
		(GridSquare){square_coordinate.lat-1,square_coordinate.lon+1},
		(GridSquare){square_coordinate.lat-1,square_coordinate.lon  },
		(GridSquare){square_coordinate.lat-1,square_coordinate.lon-1},
		(GridSquare){square_coordinate.lat  ,square_coordinate.lon-1}};

	for(int i = 0; i<9; i++){
		try{
			vector<RoughReservoir> temp = read_rough_reservoir_data(convert_string(file_storage_location+"processing_files/reservoirs/"+ocean_prefix+str(neighbors[i])+"_reservoirs_data.csv"));
			for(uint j = 0; j<temp.size(); j++)
				lower_reservoirs.push_back(temp[j]);
		}catch(int e){
			if(display)
				printf("Could not import reservoirs from %s\n", convert_string(file_storage_location+"processing_files/reservoirs/"+ocean_prefix+str(neighbors[i])+"_reservoirs_data.csv"));
		}
	}
	if(display)
		printf("Read in %zu lowers\n", lower_reservoirs.size());

	mkdir(convert_string(file_storage_location+"output/pairs"),0777);
	FILE *csv_file = fopen(convert_string(file_storage_location+"output/pairs/"+fname+"_rough_pairs.csv"), "w");
	if (!csv_file) {
	 	fprintf(stderr, "Failed to open reservoir pair CSV file\n");
		exit(1);
    }
	write_rough_pair_csv_header(csv_file);

	mkdir(convert_string(file_storage_location+"processing_files/pairs"),0777);
	FILE *csv_data_file = fopen(convert_string(file_storage_location+"processing_files/pairs/"+fname+"_rough_pairs_data.csv"), "w");
	if (!csv_data_file) {
	 	fprintf(stderr, "Failed to open reservoir pair CSV data file\n");
		exit(1);
    }
	write_rough_pair_data_header(csv_data_file);

	pairing(upper_reservoirs, lower_reservoirs, csv_file, csv_data_file);
	if(brownfield)
		pairing(lower_reservoirs, upper_reservoirs, csv_file, csv_data_file);

	int total=0;
	for (uint itest=0; itest<tests.size(); itest++) {
		if(display)
			printf("%d %.1fGWh %dh pairs\n", pairs[itest], tests[itest].energy_capacity, tests[itest].storage_time);
		total+=pairs[itest];
	}

	fclose(csv_file);
	fclose(csv_data_file);

	printf("Pairing finished for %s. Found %d pairs. Runtime: %.2f sec\n", convert_string(fname), total, 1.0e-6*(walltime_usec() - t_usec) );
}
