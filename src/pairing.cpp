#include <stdlib.h>

#include "phes_base.h"


array<vector<Pair>, tests.size()> pairs;

vector<GeographicCoordinate> find_points_to_test(ArrayCoordinate boundary[][ndirections], double wall_height, ArrayCoordinate pour_point)
{
	ArrayCoordinate one_point[ndirections] = { pour_point, pour_point, pour_point, pour_point,
					  pour_point, pour_point, pour_point, pour_point };

	int i = 0;
	while (dam_wall_heights[i] < wall_height) {
		i+=1;
	}

	int lower_wall_height = (i) ? dam_wall_heights[i-1] : 0;
	ArrayCoordinate *lower_shape = (i) ? boundary[i-1] : one_point;
	vector<GeographicCoordinate> bound;

	double inv_wall_height_interval = 0.1;
	
	for (int j=0; j<ndirections; j++) {
		GeographicCoordinate point1 = convert_coordinates(lower_shape[j]);
		GeographicCoordinate point2 = convert_coordinates(boundary[i][j]);
		bound.push_back((GeographicCoordinate) { point1.lat+(point2.lat-point1.lat)*(wall_height-lower_wall_height)*inv_wall_height_interval,
				     point1.lon+(point2.lon-point1.lon)*(wall_height-lower_wall_height)*inv_wall_height_interval });
	}

	bound.push_back(convert_coordinates(pour_point));

	return bound;
}


double find_least_distance_sqd(ArrayCoordinate upper_boundary[][ndirections], ArrayCoordinate (*lower_boundary)[ndirections],
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

Pair *check_good_pair(RoughReservoir upper, RoughReservoir lower, int energy_capacity, int storage_time, Pair *pair)
{
	int head = upper.elevation - lower.elevation;
	double required_volume = find_required_volume(energy_capacity, head);
	if ( (max_over_wall_heights(upper.volumes) < required_volume) ||
	     (max_over_wall_heights(lower.volumes) < required_volume) )
		return NULL;
	
	double upper_dam_wall_height = linear_interpolate(required_volume, upper.volumes, dam_wall_heights);
	double lower_dam_wall_height = linear_interpolate(required_volume, lower.volumes, dam_wall_heights);

	if ( (upper_dam_wall_height > upper.max_dam_height) || (lower_dam_wall_height > lower.max_dam_height) )
		return NULL;

	double upper_water_rock_estimate = required_volume/linear_interpolate(upper_dam_wall_height, dam_wall_heights, upper.dam_volumes);
	double lower_water_rock_estimate = required_volume/linear_interpolate(lower_dam_wall_height, dam_wall_heights, lower.dam_volumes);

	if ( (upper_water_rock_estimate*lower_water_rock_estimate) < 1.5*(upper_water_rock_estimate+lower_water_rock_estimate) )
		return NULL;

	ArrayCoordinate upper_coordinates = upper.pour_point;
	ArrayCoordinate lower_coordinates = lower.pour_point;

	double least_distance = find_least_distance_sqd(upper.shape_bound, lower.shape_bound,
							   upper_dam_wall_height, lower_dam_wall_height,
							   upper_coordinates, lower_coordinates);

	if ( SQ(10*head*0.001) < least_distance )
		return NULL;

	Reservoir upper_reservoir = Reservoir_init(upper.pour_point, upper.elevation);
	Reservoir lower_reservoir = Reservoir_init(lower.pour_point, lower.elevation);

	upper_reservoir.identifier = upper.identifier;
	upper_reservoir.volume = required_volume;
	upper_reservoir.dam_volume = linear_interpolate(upper_dam_wall_height, dam_wall_heights, upper.dam_volumes);
	upper_reservoir.water_rock = upper_water_rock_estimate;
	upper_reservoir.dam_height = upper_dam_wall_height;
	upper_reservoir.max_dam_height = upper.max_dam_height;

	lower_reservoir.identifier = lower.identifier;
	lower_reservoir.volume = required_volume;
	lower_reservoir.dam_volume = linear_interpolate(lower_dam_wall_height, dam_wall_heights, lower.dam_volumes);
	lower_reservoir.water_rock = lower_water_rock_estimate;
	lower_reservoir.dam_height = lower_dam_wall_height;
	lower_reservoir.max_dam_height = lower.max_dam_height;

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
	set_FOM(pair);
	return pair;
}

void pairing(vector<RoughReservoir> upper_reservoirs, vector<RoughReservoir> lower_reservoirs)
{
	RoughReservoir* upper_reservoir;
	RoughReservoir* lower_reservoir;
	for (uint iupper=0;iupper < upper_reservoirs.size(); iupper++) {
		upper_reservoir = &upper_reservoirs[iupper];
		double coslat = COS(upper_reservoir->latitude);
		for (uint ilower = 0; ilower<lower_reservoirs.size(); ilower++){
		    lower_reservoir = &lower_reservoirs[ilower];

		    int head = upper_reservoir->elevation - lower_reservoir->elevation;
			if ( head < MIN_HEAD)
				continue;

			// Pour point separation
			double distance_sqd = find_distance_sqd(upper_reservoir->pour_point,lower_reservoir->pour_point, coslat);

			if (distance_sqd >= SQ(20) || ( SQ(20.0*head*0.001) <= distance_sqd ))
				continue;

			for (uint itest=0; itest<tests.size(); itest++) {
				Pair temp_pair;
				if (check_good_pair(*upper_reservoir, *lower_reservoir, tests[itest].energy_capacity, tests[itest].storage_time, &temp_pair)) {
					pairs[itest].push_back(temp_pair);
				}
			}
		}
	}
}


int main(int nargs, char **argv)
{
	GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));

	unsigned long t_usec = walltime_usec();
	vector<RoughReservoir> upper_reservoirs = read_rough_reservoir_data(convert_string("processing_files/reservoirs/"+str(square_coordinate)+"_reservoirs_data.csv"));
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
			vector<RoughReservoir> temp = read_rough_reservoir_data(convert_string("processing_files/reservoirs/"+str(neighbors[i])+"_reservoirs_data.csv"));
			for(uint j = 0; j<temp.size(); j++)
				lower_reservoirs.push_back(temp[j]);
		}catch(int e){
			printf("Could not import reservoirs from %s\n", convert_string("processing_files/reservoirs/"+str(neighbors[i])+"_reservoirs_data.csv"));
		}
	}

	printf("Read in %zu lowers\n", lower_reservoirs.size());

	pairing(upper_reservoirs, lower_reservoirs);

	FILE *csv_file = fopen(convert_string("Output/Pairs/"+str(square_coordinate)+"_rough_pairs.csv"), "w");
	if (!csv_file) {
	 	fprintf(stderr, "failed to open reservoir pair CSV file\n");
		exit(1);
    }
	write_rough_pair_csv_header(csv_file);

	FILE *csv_data_file = fopen(convert_string("processing_files/pairs/"+str(square_coordinate)+"_rough_pairs_data.csv"), "w");
	if (!csv_data_file) {
	 	fprintf(stderr, "failed to open reservoir pair CSV data file\n");
		exit(1);
    }
	write_rough_pair_data_header(csv_data_file);

	int total=0;
	for (uint itest=0; itest<tests.size(); itest++) {

		printf("%zu reservoir pairs at energy %dGWh and storage time %dh\n", pairs[itest].size(), tests[itest].energy_capacity, tests[itest].storage_time);

		total+=pairs[itest].size();
		Pair pair;
		for (uint i=0 ; i<pairs[itest].size(); i++) {
			pair = pairs[itest][i];
			write_rough_pair_csv(csv_file, &pair);
			write_rough_pair_data(csv_data_file, &pair);
		}
	}

	fclose(csv_file);

	printf("Reservoir pairs: found %d pairs,  Runtime: %.2f sec\n", total, 1.0e-6*(walltime_usec() - t_usec) );
}
