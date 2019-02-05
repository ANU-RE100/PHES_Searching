#include "phes_base.h"
#include "kml.h"

int display = false;

vector<vector<Pair>> pairs;

// Returns a tuple with the two cells adjacent to an edge defined by two points
ArrayCoordinate* get_adjacent_cells(ArrayCoordinate point1, ArrayCoordinate point2){
    double average_row = (point1.row+point2.row)/2.0;
    double average_col = (point1.col+point2.col)/2.0;
    
    ArrayCoordinate* points = (ArrayCoordinate*)malloc(2*sizeof(ArrayCoordinate));

    if(average_row-(int)average_row>0.9 || average_row-(int)average_row<0.1){
        points[0] = ArrayCoordinate_init((int)(average_row+EPS),(int)(average_col-0.5+EPS), point1.origin);
    	points[1] = ArrayCoordinate_init((int)(average_row-1+EPS),(int)(average_col-0.5+EPS), point1.origin);
    }else{
    	points[0] = ArrayCoordinate_init((int)(average_row-0.5+EPS),(int)(average_col+EPS), point1.origin);
    	points[1] = ArrayCoordinate_init((int)(average_row-0.5+EPS),(int)(average_col-1+EPS), point1.origin);
    }
    return points;
}

// Determines if two points give the edge of a reservoir given its raster model
bool is_edge(ArrayCoordinate point1, ArrayCoordinate point2, Model_int16* model, int threshold){
    
    if (point1.row<0 || point1.col<0 || point1.row>model->shape[0] || point1.col>model->shape[1])
        return false;
    if (point2.row<0 || point2.col<0 || point2.row>model->shape[0] || point2.col>model->shape[1])
        return false;
    
    ArrayCoordinate* to_check = get_adjacent_cells(point1, point2);

    if((!check_within(to_check[0], model->shape) && model->d[to_check[1].row][to_check[1].col]>=threshold)||(!check_within(to_check[1], model->shape) && model->d[to_check[0].row][to_check[0].col]>=threshold))
        return true;

    if (model->d[to_check[0].row][to_check[0].col]<threshold && model->d[to_check[1].row][to_check[1].col]>=threshold)
        return true;
    if (model->d[to_check[1].row][to_check[1].col]<threshold && model->d[to_check[0].row][to_check[0].col]>=threshold)
        return true;
    
    return false;
}

// Determines if an edge of a reservoir between two points requires a dam wall
bool is_dam_wall(ArrayCoordinate point1, ArrayCoordinate point2, Model_int16* DEM, double wall_elevation){
    if (point1.row<0 || point1.col<0 || point1.row>DEM->shape[0] || point1.col>DEM->shape[1])
        return false;
    if (point2.row<0 || point2.col<0 || point2.row>DEM->shape[0] || point2.col>DEM->shape[1])
        return false;
    
    double average_row = (point1.row+point2.row)/2.0;
    double average_col = (point1.col+point2.col)/2.0;
    
    ArrayCoordinate* to_check = get_adjacent_cells(point1, point2);

    if(average_row-(int)average_row>0.9 || average_row-(int)average_row<0.1){
        to_check[0] = ArrayCoordinate_init((int)(average_row+EPS),(int)(average_col-0.5+EPS), point1.origin);
    	to_check[1] = ArrayCoordinate_init((int)(average_row-1+EPS),(int)(average_col-0.5+EPS), point1.origin);
    }else{
    	to_check[0] = ArrayCoordinate_init((int)(average_row-0.5+EPS),(int)(average_col+EPS), point1.origin);
    	to_check[1] = ArrayCoordinate_init((int)(average_row-0.5+EPS),(int)(average_col-1+EPS), point1.origin);
    }

    if(!check_within(to_check[0], DEM->shape) || !check_within(to_check[0], DEM->shape))
        return false;

    if(DEM->d[to_check[0].row][to_check[0].col]<wall_elevation and DEM->d[to_check[1].row][to_check[1].col]<wall_elevation)
        return true;
    return false;
}

// Converts a raster model to a polygon given a raster model and a point on the interior edge of the polygon
vector<ArrayCoordinate> convert_to_polygon(Model_int16* model, ArrayCoordinate pour_point, int threshold){
    
    vector<ArrayCoordinate> to_return;
    int dir_def[4][2] = {{-1,0},{1,0},{0,-1},{0,1}};
    int dir_to_do[4][3] = {{2,0,3},{3,1,2},{1,2,0},{0,3,1}};

    int tests[] = {1, 2, 0, 3};
    ArrayCoordinate test_coordinates[] = {
    	ArrayCoordinate_init(pour_point.row+1, pour_point.col+1, pour_point.origin),
    	ArrayCoordinate_init(pour_point.row+1, pour_point.col, pour_point.origin),
    	ArrayCoordinate_init(pour_point.row, pour_point.col, pour_point.origin),
    	ArrayCoordinate_init(pour_point.row, pour_point.col+1, pour_point.origin)};

    for(int i = 0; i<4; i++){
    	vector<ArrayCoordinate> temp_to_return;
        int last_dir = tests[i];
        ArrayCoordinate initial = test_coordinates[i];
        ArrayCoordinate last = initial;
        while(true){
        	for(int id = 0; id<3; id++){
        		int d = dir_to_do[last_dir][id];
        		ArrayCoordinate next = ArrayCoordinate_init(last.row+dir_def[d][0],last.col+dir_def[d][1], pour_point.origin);
        		if(is_edge(last, next, model, threshold)){
        			temp_to_return.push_back(next);
        			last = next;
        			last_dir = d;
        			break;
        		}
        	}
        	if(last.row==initial.row && last.col == initial.col)
        		break;
        }
        if(temp_to_return.size()>to_return.size()){
        	to_return.clear();
        	for(uint j = 0; j<temp_to_return.size(); j++){
        		to_return.push_back(temp_to_return[j]);
        	}
        }
    }
    return to_return;
}


vector<GeographicCoordinate> convert_poly(vector<ArrayCoordinate> polygon){
    vector<GeographicCoordinate> to_return;
    for(uint i = 0; i<polygon.size(); i++){
    	to_return.push_back(convert_coordinates(polygon[i]));
    }
    return to_return;
}

// "Smooths" a polygon by cutting the corners
vector<GeographicCoordinate> corner_cut_poly(vector<GeographicCoordinate> polygon){
    vector<GeographicCoordinate> to_return;
    for(uint i = 0; i<polygon.size(); i++){
    	to_return.push_back(GeographicCoordinate_init((polygon[i].lat*2+polygon[(i+1)%polygon.size()].lat*2)/4.0, (polygon[i].lon*2+polygon[(i+1)%polygon.size()].lon*2)/4.0));
    }
    to_return.push_back(to_return[0]);
    return to_return;
}

// Compresses a polygon by removing collinear points
vector<GeographicCoordinate> compress_poly(vector<GeographicCoordinate> polygon){
    vector<GeographicCoordinate> to_return;
    to_return.push_back(polygon[0]);
    for(uint i = 1; i<polygon.size()-1; i++){
        if ((polygon[i+1].lon-polygon[i].lon)/(polygon[i+1].lat-polygon[i].lat) < (polygon[i].lon-polygon[i-1].lon)/(polygon[i].lat-polygon[i-1].lat)-EPS 
        	|| (polygon[i+1].lon-polygon[i].lon)/(polygon[i+1].lat-polygon[i].lat) > (polygon[i].lon-polygon[i-1].lon)/(polygon[i].lat-polygon[i-1].lat)+EPS ){
            to_return.push_back(polygon[i]);
        }
    }
    to_return.push_back(polygon[polygon.size()-1]);                           
    return to_return;
}

string str(vector<GeographicCoordinate> polygon, double elevation){
	string to_return = "";
	for(uint i = 0; i<polygon.size(); i++){
		to_return += dtos(polygon[i].lon, 5)+","+dtos(polygon[i].lat, 5)+","+dtos(elevation, 1)+" ";
	}
	return to_return;
}

bool model_reservoir(Reservoir* reservoir, Models models, Reservoir_KML_Coordinates* coordinates, Model<bool>* seen, bool* non_overlap, vector<ArrayCoordinate> *used_points, BigModel big_model){
	
	Model_int16 *DEM = models.DEMs[0];
	Model<char>* flow_directions = big_model.flow_directions[0];
	Model_int16 *cur_model = Model_int16_create(models.DEMs[0]->shape, MODEL_SET_ZERO);

	for(int i = 0; i<9; i++)
		if(models.neighbors[i].lat == (int)(FLOOR(reservoir->latitude)-EPS) && models.neighbors[i].lon == (int)(FLOOR(reservoir->longitude)+EPS)){
			DEM = models.DEMs[i];
			flow_directions = big_model.flow_directions[i];
		}

	double req_volume = reservoir->volume;
	reservoir->volume = 0;
	reservoir->area = 0;
	vector<ArrayCoordinate> temp_used_points;
	
	// RESERVOIR
	char last_dir = 'd';
	while(reservoir->volume*(1+0.5/reservoir->water_rock)<(1-volume_accuracy)*req_volume || reservoir->volume*(1+0.5/reservoir->water_rock)>(1+volume_accuracy)*req_volume){
		temp_used_points.clear();
		reservoir->volume = 0;
		reservoir->area = 0;

		queue<ArrayCoordinate> q;
		q.push(reservoir->pour_point);
		while (!q.empty()) {
			ArrayCoordinate p = q.front();
			q.pop();

			cur_model->d[p.row][p.col]=1;

			ArrayCoordinate big_ac = convert_coordinates(convert_coordinates(p), models.origin);

			temp_used_points.push_back(big_ac);
			if (seen->get(big_ac.row,big_ac.col)){
				*non_overlap = false;
			}

			reservoir->volume+=(reservoir->dam_height-(DEM->d[p.row][p.col]-DEM->d[reservoir->pour_point.row][reservoir->pour_point.col]))*find_area(p)/100;
			reservoir->area+=find_area(p);
			update_reservoir_boundary(reservoir->shape_bound, p);

			for (uint d=0; d<directions.size(); d++) {
				ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
				if (flow_directions->check_within(neighbor.row, neighbor.col) &&
				    flows_to(neighbor, p, flow_directions) &&
				    ((DEM->d[neighbor.row][neighbor.col]-DEM->d[reservoir->pour_point.row][reservoir->pour_point.col]) < reservoir->dam_height) ) {
					q.push(neighbor);
				}
			}
		}
		
		if(reservoir->volume*(1+0.5/reservoir->water_rock)<(1-volume_accuracy)*req_volume){
			reservoir->dam_height+=dam_wall_height_resolution;
            if(reservoir->dam_height>reservoir->max_dam_height)
                return false;
            last_dir = 'u';
		}
                
        if(reservoir->volume*(1+0.5/reservoir->water_rock)>(1+volume_accuracy)*req_volume){
        	if (last_dir == 'u')
                break;
            reservoir->dam_height-=dam_wall_height_resolution;
            last_dir = 'd';
        }            
	}
	vector<ArrayCoordinate> reservoir_polygon = convert_to_polygon(cur_model, reservoir->pour_point, 1);
	for(uint i = 0; i<temp_used_points.size(); i++){
		used_points->push_back(temp_used_points[i]);
	}
	
	// DAM
	reservoir->dam_volume = 0;
	reservoir->dam_length = 0;
	bool is_turkeys_nest = true;
	vector<vector<ArrayCoordinate> > dam_polygon;
	bool last = false;
	vector<bool> polygon_bool;
	for(uint i = 0; i<reservoir_polygon.size(); i++){
		ArrayCoordinate point1 = reservoir_polygon[i];
		ArrayCoordinate point2 = reservoir_polygon[(i+1)%reservoir_polygon.size()];
		if(is_dam_wall(point1, point2, DEM, reservoir->elevation+reservoir->dam_height)){
            polygon_bool.push_back(true);
            
            ArrayCoordinate* adjacent = get_adjacent_cells(point1, point2);
            double average_height = (DEM->d[adjacent[0].row][adjacent[0].col]+DEM->d[adjacent[1].row][adjacent[1].col])/2.0;
            if(cur_model->d[adjacent[0].row][adjacent[0].col] == 0)
                cur_model->d[adjacent[0].row][adjacent[0].col] = 2;
            else
            	cur_model->d[adjacent[1].row][adjacent[1].col] = 2;
            if(!last){
            	vector<ArrayCoordinate> temp;
            	temp.push_back(point1);
            	dam_polygon.push_back(temp);
            }
            double height = reservoir->elevation+reservoir->dam_height-average_height;
            double length = find_distance(point1, point2)*1000;
            reservoir->dam_volume += convert_to_dam_volume(height, length);
            reservoir->dam_length += length;
            dam_polygon[dam_polygon.size()-1].push_back(point2);
            last = true;
        }else{
        	is_turkeys_nest = false;
        	polygon_bool.push_back(false);
        	last = false;
        }
	}
	if(polygon_bool[0] && polygon_bool[dam_polygon.size()-1] && !is_turkeys_nest){
		for(uint i = 0; i<dam_polygon[dam_polygon.size()-1].size();i++){
			dam_polygon[0].push_back(dam_polygon[dam_polygon.size()-1][i]);
		}
        dam_polygon.pop_back();
    }

    if(is_turkeys_nest){
    	ArrayCoordinate* adjacent = get_adjacent_cells(dam_polygon[0][0], dam_polygon[0][1]);
    	ArrayCoordinate to_check = adjacent[1];
    	if(cur_model->d[adjacent[0].row][adjacent[0].col]==2)
    		to_check = adjacent[0];
    	vector<GeographicCoordinate> polygon = compress_poly(corner_cut_poly(convert_poly(convert_to_polygon(cur_model, to_check, 2))));
    	string polygon_string = str(polygon, reservoir->elevation+reservoir->dam_height+freeboard);
    	coordinates->dam.push_back(polygon_string);
    }else{
    	for(uint i = 0; i<dam_polygon.size();i++){
    		ArrayCoordinate* adjacent = get_adjacent_cells(dam_polygon[i][0], dam_polygon[i][1]);
	    	ArrayCoordinate to_check = adjacent[1];
	    	if(cur_model->d[adjacent[0].row][adjacent[0].col]==2)
	    		to_check = adjacent[0];
	    	vector<GeographicCoordinate> polygon = compress_poly(corner_cut_poly(convert_poly(convert_to_polygon(cur_model, to_check, 2))));
	    	string polygon_string = str(polygon, reservoir->elevation+reservoir->dam_height+freeboard);
	    	coordinates->dam.push_back(polygon_string);
    	}
    }

    reservoir->volume+=(reservoir->dam_volume)/2;
    reservoir->water_rock = reservoir->volume/reservoir->dam_volume;
    reservoir->average_water_depth = reservoir->volume/reservoir->area;

    string polygon_string = str(compress_poly(corner_cut_poly(convert_poly(reservoir_polygon))), reservoir->elevation+reservoir->dam_height);
    coordinates->reservoir = polygon_string;

    coordinates->is_turkeys_nest = is_turkeys_nest;
    Model_int16_free(cur_model);
	return true;
}

bool model_pair(Pair* pair, Models models, Pair_KML* pair_kml, Model<bool>* seen, bool* non_overlap, int max_FOM, BigModel big_model){

	vector<ArrayCoordinate> used_points;
	*non_overlap = true;

	if(!model_reservoir(&pair->upper, models, &pair_kml->upper, seen, non_overlap, &used_points, big_model))
		return false;

	if(!model_reservoir(&pair->lower, models, &pair_kml->lower, seen, non_overlap, &used_points, big_model))
		return false;

	if(*non_overlap){
		for(uint i = 0; i<used_points.size();i++){
			seen->set(used_points[i].row,used_points[i].col,true);
		}
	}

	ArrayCoordinate upper_closest_point = pair->upper.pour_point;
    ArrayCoordinate lower_closest_point = pair->lower.pour_point;
    double mindist = find_distance(upper_closest_point, lower_closest_point);
    for(uint iupper = 0; iupper<directions.size(); iupper++)
    	for(uint ilower = 0; ilower<directions.size(); ilower++)
    		if (find_distance(pair->upper.shape_bound[iupper], pair->lower.shape_bound[ilower])<mindist){
                upper_closest_point = pair->upper.shape_bound[iupper];
                lower_closest_point = pair->lower.shape_bound[ilower];
                mindist = find_distance(pair->upper.shape_bound[iupper], pair->lower.shape_bound[ilower]);
            }

    pair->distance = mindist;
    pair->slope = pair->head/(pair->distance*1000);
    pair->volume = min(pair->upper.volume, pair->lower.volume);
    pair->water_rock = 1/((1/pair->upper.water_rock)+(1/pair->lower.water_rock));
    set_FOM(pair);
    if(pair->FOM>max_FOM || pair->category=='Z'){
        return false;
    }
    GeographicCoordinate average = GeographicCoordinate_init((convert_coordinates(upper_closest_point).lat+convert_coordinates(lower_closest_point).lat)/2,
    	((convert_coordinates(upper_closest_point).lon+convert_coordinates(lower_closest_point).lon)/2));
    pair_kml->point = dtos(average.lon,5)+","+dtos(average.lat,5)+",0";
    pair_kml->line = dtos(convert_coordinates(upper_closest_point).lon,5)+","+dtos(convert_coordinates(upper_closest_point).lat,5)+",0 "+dtos(convert_coordinates(lower_closest_point).lon,5)+","+dtos(convert_coordinates(lower_closest_point).lat,5)+",0";
	return true;
}

int main(int nargs, char **argv)
{

	GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));
    if(nargs>3)
        display = atoi(argv[3]);

    printf("Constructor started for %s\n",convert_string(str(square_coordinate)));

    TIFF_IO_init();
    parse_variables(convert_string(file_storage_location+"variables"));
    unsigned long t_usec = walltime_usec();

    BigModel big_model = BigModel_init(square_coordinate);
	Models models = Models_init(square_coordinate);

	for(int i = 0; i<9; i++){
		GridSquare gs = models.neighbors[i];
		try{
			models.DEMs[i] = read_DEM_with_borders(gs);
			models.flow_directions[i] = TIFF_Read_int16(convert_string(file_storage_location+"processing_files/flow_directions/"+str(gs)+"_flow_directions.tif"), NULL, NULL);
		}catch(int e){
            if(display)
			 printf("Could not find %s\n", convert_string(str(gs)));
		}
		
	}

	Model<bool>* seen = new Model<bool>(big_model.DEM->nrows(), big_model.DEM->nrows(), MODEL_SET_ZERO);
    seen->set_geodata(big_model.DEM->get_geodata());

	pairs = read_rough_pair_data(convert_string(file_storage_location+"processing_files/pretty_set_pairs/"+str(square_coordinate)+"_rough_pretty_set_pairs_data.csv"));
	mkdir(convert_string(file_storage_location+"output/final_output"), 0777);
	mkdir(convert_string(file_storage_location+"output/final_output/"+str(square_coordinate)),0777);

	FILE *total_csv_file = fopen(convert_string(file_storage_location+"output/final_output/"+str(square_coordinate)+"/"+str(square_coordinate)+"_total.csv"), "w");
	write_total_csv_header(total_csv_file);

	int total_count = 0;
	int total_capacity = 0;
	for(uint i = 0; i<tests.size(); i++){
        sort(pairs[i].begin(), pairs[i].end());
        
		FILE *csv_file = fopen(convert_string(file_storage_location+"output/final_output/"+str(square_coordinate)+"/"+str(square_coordinate)+"_"+str(tests[i])+".csv"), "w");
		write_pair_csv_header(csv_file);

		ofstream kml_file(convert_string(file_storage_location+"output/final_output/"+str(square_coordinate)+"/"+str(square_coordinate)+"_"+str(tests[i])+".kml"), ios::out);
		KML_Holder kml_holder;

		//FILE *fusion_csv_file = fopen(convert_string(file_storage_location+"Output/Final Output/"+str(square_coordinate)+"/"+str(square_coordinate)+"_"+str(tests[i])+"_Fusion_Table.csv"), "w");
		//write_fusion_csv_header(fusion_csv_file);

		sort(pairs[i].begin(), pairs[i].end());
		int count = 0;
		for(uint j=0; j<pairs[i].size(); j++){
			Pair_KML pair_kml;
			bool non_overlap;
			if(model_pair(&pairs[i][j], models, &pair_kml, seen, &non_overlap, tests[i].max_FOM, big_model)){
				write_pair_csv(csv_file, &pairs[i][j]);
				//write_fusion_csv(fusion_csv_file, &pairs[i][j], &pair_kml);
				update_kml_holder(&kml_holder, &pairs[i][j], &pair_kml);
				count++;
				if(non_overlap){
					total_count++;
					total_capacity+=tests[i].energy_capacity;
				}
			}
		}

		kml_file << output_kml(&kml_holder, str(square_coordinate), tests[i]);
        if(display)
            printf("%d %dGWh %dh Pairs\n", count, tests[i].energy_capacity, tests[i].storage_time);
		kml_file.close();
        fclose(csv_file);
	}
	write_total_csv(total_csv_file, str(square_coordinate), total_count, total_capacity);
    printf("Constructor finished for %s. Found %d non-overlapping pairs with a total of %dGWh. Runtime: %.2f sec\n", convert_string(str(square_coordinate)), total_count, total_capacity, 1.0e-6*(walltime_usec() - t_usec) );
}
