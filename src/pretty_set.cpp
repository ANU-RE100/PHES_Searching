#include <stdlib.h>
#include <sys/time.h>
#include "shapefil.h"
#include <sys/stat.h> 
#include <sys/types.h> 

#include "model2D.h"
#include "TIFF_IO.h"
#include "phes_base.h"

//Being risky (for devel only)
#include <bits/stdc++.h>
using namespace std;

int display = false;

vector<vector<Pair> > pairs;

bool check_reservoir(Reservoir reservoir, Models models, Model_int8 *seen, vector<ArrayCoordinate> *used_points){
	Model_int16 *DEM = models.DEMs[0];
	Model_int16 *flow_directions = models.flow_directions[0];

	for(int i = 0; i<9; i++){
		if(models.neighbors[i].lat == (int)(FLOOR(reservoir.latitude)-EPS) && models.neighbors[i].lon == (int)(FLOOR(reservoir.longitude)+EPS)){
			DEM = models.DEMs[i];
			flow_directions = models.flow_directions[i];
		}
	}

	double volume = 0;
	double req_volume = reservoir.volume;
	double wall_height = reservoir.dam_height;
	vector<ArrayCoordinate> temp_used_points;

	char last = 'd';

	while(volume*(1+0.5/reservoir.water_rock)<0.99*req_volume || volume*(1+0.5/reservoir.water_rock)>1.01*req_volume){
		temp_used_points.clear();
		volume = 0;
		queue<ArrayCoordinate> q;
		q.push(reservoir.pour_point);
		while (!q.empty()) {
			ArrayCoordinate p = q.front();
			q.pop();
			ArrayCoordinate big_ac = convert_coordinates(convert_coordinates(p), models.origin);

			temp_used_points.push_back(big_ac);

			if (seen->d[big_ac.row][big_ac.col]){
				return false;
			}
			volume+=(wall_height-(DEM->d[p.row][p.col]-DEM->d[reservoir.pour_point.row][reservoir.pour_point.col]))*find_area(p)/100;
			for (uint d=0; d<directions.size(); d++) {
				ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
				if (check_within(neighbor, flow_directions->shape) &&
				    flows_to(neighbor, p, flow_directions) &&
				    ((DEM->d[neighbor.row][neighbor.col]-DEM->d[reservoir.pour_point.row][reservoir.pour_point.col]) < wall_height) ) {
					q.push(neighbor);
				}
			}

		}
		
		if(volume*(1+0.5/reservoir.water_rock)<0.99*req_volume){
			wall_height+=0.1;
            if(wall_height>reservoir.max_dam_height)
                return false;
            last = 'u';
		}
                
        if(volume*(1+0.5/reservoir.water_rock)>1.01*req_volume){
        	if (last == 'u')
                break;
            wall_height-=0.1;
            last = 'd';
        }            
	}
	for(uint i = 0; i<temp_used_points.size(); i++){
		used_points->push_back(temp_used_points[i]);
	}
	return true;
}

bool check_pair(Pair pair, Models models, Model_int8 *seen){
	vector<ArrayCoordinate> used_points;
	if(!check_reservoir(pair.upper, models, seen, &used_points))
		return false;
	if(!check_reservoir(pair.lower, models, seen, &used_points))
		return false;
	for(uint i = 0; i<used_points.size();i++){
		seen->d[used_points[i].row][used_points[i].col] = true;
	}
	return true;
}

int main(int nargs, char **argv)
{

	GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));
	if(nargs>3)
		display = atoi(argv[3]);

	printf("Pretty set started for %s\n",convert_string(str(square_coordinate)));

	TIFF_IO_init();
	parse_variables(convert_string("variables"));
	unsigned long t_usec = walltime_usec();
	
	Models models = Models_init(square_coordinate);

	for(int i = 0; i<9; i++){
		GridSquare gs = models.neighbors[i];
		try{
			models.DEMs[i] = read_DEM_with_borders(gs);
			models.flow_directions[i] = TIFF_Read_int16(convert_string("processing_files/flow_directions/"+str(gs)+"_flow_directions.tif"), NULL, NULL);
		}catch(int e){
			if(display)
				printf("Could not find %s\n", convert_string(str(gs)));
		}
		
	}

	pairs = read_rough_pair_data(convert_string("processing_files/pairs/"+str(square_coordinate)+"_rough_pairs_data.csv"));

	Model_int8 *seen;

	mkdir("processing_files/pretty_set_pairs",0777);
	FILE *csv_data_file = fopen(convert_string("processing_files/pretty_set_pairs/"+str(square_coordinate)+"_rough_pretty_set_pairs_data.csv"), "w");
	write_rough_pair_data_header(csv_data_file);

	for(uint i = 0; i<tests.size(); i++){
		sort(pairs[i].begin(), pairs[i].end());
		int big_shape[2] = {models.DEMs[0]->shape[0]*3-border*4, models.DEMs[0]->shape[1]*3-border*4};
		seen = Model_int8_create(big_shape, MODEL_SET_ZERO);
		int count = 0;
		for(uint j=0; j<pairs[i].size(); j++){
			if(check_pair(pairs[i][j], models, seen)){
				write_rough_pair_data(csv_data_file, &pairs[i][j]);
				count++;
			}
		}
		if(display)
			printf("%d %dGWh Pairs with storage time %dh\n", count, tests[i].energy_capacity, tests[i].storage_time);
	}
	printf("Pretty set finished for %s. Runtime: %.2f sec\n", convert_string(str(square_coordinate)), 1.0e-6*(walltime_usec() - t_usec) );
}
