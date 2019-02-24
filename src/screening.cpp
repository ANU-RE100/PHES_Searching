#include "phes_base.h"

bool debug_output = false;
bool debug = false;
int display = false;

// find_polygon_intersections returns an array containing the longitude of all line. Assumes last coordinate is same as first
vector<double> find_polygon_intersections(int row, vector<GeographicCoordinate> &polygon, Model<bool>* filter){
    vector<double> to_return;
    double lat = filter->get_coordinate(row, 0).lat;
    for(uint i = 0; i<polygon.size()-1; i++){
        GeographicCoordinate line[2] = {polygon[i], polygon[(i+1)]};
        if((line[0].lat < lat && line[1].lat>=lat) || (line[0].lat >= lat && line[1].lat < lat)){
            to_return.push_back(line[0].lon+(lat-line[0].lat)/(line[1].lat-line[0].lat)*(line[1].lon-line[0].lon));
        }
    }
    sort(to_return.begin(), to_return.end());
    return to_return;
}

void read_shp_filter(string filename, Model<bool>* filter){
	char *shp_filename = new char[filename.length() + 1];
	strcpy(shp_filename, filename.c_str());
    if(!file_exists(shp_filename)){
		if(display)
			cout << "No file: "+filename;
		throw(1);
	}
	SHPHandle SHP = SHPOpen(convert_string(filename), "rb" );
	if(SHP != NULL ){
    	int	nEntities;
    	vector<vector<GeographicCoordinate>> relevant_polygons;
    	SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL );
	    for( int i = 0; i < nEntities; i++ )
	    {
	        SHPObject	*shape;
	        shape = SHPReadObject( SHP, i );
	        if( shape == NULL ){
	            fprintf( stderr,"Unable to read shape %d, terminating object reading.\n",i);
	            break;
	        }
	        vector<GeographicCoordinate> temp_poly;
	        bool to_keep = false;
	        for(int j = 0, iPart = 1; j < shape->nVertices; j++ )
	        {
	            if( iPart < shape->nParts && shape->panPartStart[iPart] == j )
	            {
	            	if(to_keep)
	        			relevant_polygons.push_back(temp_poly);
	            	to_keep = false;
	            	temp_poly.clear();
	                iPart++;
	            }
	            GeographicCoordinate temp_point = GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
	            to_keep = (to_keep || filter->check_within(temp_point));
	            temp_poly.push_back(temp_point);
	        }
	        if(to_keep)
	        	relevant_polygons.push_back(temp_poly);
	        SHPDestroyObject( shape );
	    }
	    if(display)
	    	printf("%d Polygons imported from %s\n", (int)relevant_polygons.size(), convert_string(filename));
	    for(uint i = 0; i<relevant_polygons.size(); i++){
	    	for(int row =0; row<filter->nrows(); row++){
                vector<double> polygon_intersections = find_polygon_intersections(row, relevant_polygons[i], filter);
                for(uint j = 0; j<polygon_intersections.size();j++)
                	polygon_intersections[j] = (convert_coordinates(GeographicCoordinate_init(0, polygon_intersections[j]),filter->get_origin()).col);
                for(uint j = 0; j<polygon_intersections.size()/2;j++)
                    for(int col=polygon_intersections[2*j];col<polygon_intersections[2*j+1];col++)
                        if(filter->check_within(row, col))
                        	filter->set(row,col,true);
            }
	    }
    }else{
    	throw(1);
    }
    SHPClose(SHP);
}

void read_tif_filter(string filename, Model<bool>* filter, unsigned char value_to_filter){
	try{
		Model<unsigned char>* tif_filter = new Model<unsigned char>(filename, GDT_Byte);
		GeographicCoordinate point;
		for(int row = 0; row<filter->nrows(); row++){
			for(int col = 0; col<filter->ncols(); col++){
				point = filter->get_coordinate(row, col);
				if(tif_filter->check_within(point) && tif_filter->get(point)==value_to_filter)
					filter->set(row,col,true);
			}
		}
		delete tif_filter;
	}catch(exception& e){
		if(display)
			printf("Problem with %s\n", convert_string(filename));
	}catch(int e){
		if(display)
			printf("Problem with %s\n", convert_string(filename));
	}
}

string find_world_urban_filename(GeographicCoordinate point){
	char clat = 'A'+floor((point.lat+96)/8);
	if(clat>=73)clat++;
	if(clat>=79)clat++;
	int nlon = floor((point.lon+180)/6+1);
	char strlon[3];
	snprintf(strlon, 3, "%02d", nlon);
	return "input/WORLD_URBAN/"+string(strlon)+string(1, clat)+"_hbase_human_built_up_and_settlement_extent_geographic_30m.tif";
}

Model<bool>* read_filter(Model<short>* DEM, vector<string> filenames)
{
	Model<bool>* filter = new Model<bool>(DEM->nrows(), DEM->ncols(), MODEL_SET_ZERO);
	filter->set_geodata(DEM->get_geodata());
	for(string filename:filenames){
		if(filename=="use_world_urban"){
			if(display)
				printf("Using world urban data as filter\n");
			vector<string> done;
			for(GeographicCoordinate corner: DEM->get_corners()){
				string urban_filename = find_world_urban_filename(corner);
				if(find(done.begin(), done.end(), urban_filename)==done.end()){
					read_tif_filter(file_storage_location+urban_filename, filter, 201);
					done.push_back(urban_filename);
				}
			}
		}else if(filename=="use_tiled_filter"){
			if(display)
				printf("Using tiled filter\n");
			GridSquare sc = {(int)((filter->get_coordinate(filter->nrows(), filter->ncols()).lat+filter->get_origin().lat)/2.0)-1,(int)((filter->get_coordinate(filter->nrows(), filter->ncols()).lon+filter->get_origin().lon)/2.0)};
			GridSquare neighbors[9] = {
				(GridSquare){sc.lat  ,sc.lon  },
				(GridSquare){sc.lat+1,sc.lon-1},
				(GridSquare){sc.lat+1,sc.lon  },
				(GridSquare){sc.lat+1,sc.lon+1},
				(GridSquare){sc.lat  ,sc.lon+1},
				(GridSquare){sc.lat-1,sc.lon+1},
				(GridSquare){sc.lat-1,sc.lon  },
				(GridSquare){sc.lat-1,sc.lon-1},
				(GridSquare){sc.lat  ,sc.lon-1}};
			for(int i = 0; i<9; i++){
				try{
					read_shp_filter(file_storage_location+"input/shapefile_tiles/"+str(neighbors[i])+"_shapefile_tile.shp", filter);
				}catch(int e){
					if(i==0)
						exit(1);
				}
			}
		}else{
			read_shp_filter(file_storage_location+filename, filter);
		}		
	}
	for(int row = 0; row<DEM->nrows(); row++)
		for(int col = 0; col<DEM->ncols(); col++)
			if(DEM->get(row,col)<-1000)
				filter->set(row,col,true);
	return filter;
}


Model<double>* fill(Model<short>* DEM)
{
	Model<double>* DEM_filled_no_flat = new Model<double>(DEM->nrows(), DEM->ncols(), MODEL_UNSET);
	DEM_filled_no_flat->set_geodata(DEM->get_geodata());
	Model<bool>* seen = new Model<bool>(DEM->nrows(), DEM->ncols(), MODEL_SET_ZERO);

	queue<ArrayCoordinateWithHeight> q;
	priority_queue<ArrayCoordinateWithHeight> pq;
	for (int row=0; row<DEM->nrows(); row++)
		for (int col=0; col<DEM->ncols();col++)
			DEM_filled_no_flat->set(row,col,(double)DEM->get(row,col));

	for (int row=0; row<DEM->nrows()-1;row++) {
		pq.push(ArrayCoordinateWithHeight_init(row+1, DEM->ncols()-1, (double)DEM->get(row+1,DEM->ncols()-1)));
		seen->set(row+1,DEM->ncols()-1,true);
		pq.push(ArrayCoordinateWithHeight_init(row, 0, (double)DEM->get(row,0)));
		seen->set(row,0,true);
	}

	for (int col=0; col<DEM->ncols()-1;col++) {
		pq.push(ArrayCoordinateWithHeight_init(DEM->nrows()-1, col, (double)DEM->get(DEM->ncols()-1,col)));
		seen->set(DEM->nrows()-1,col,true);
		pq.push(ArrayCoordinateWithHeight_init(0, col+1, (double)DEM->get(0,col+1)));
		seen->set(0,col+1,true);
	}

	ArrayCoordinateWithHeight c;
	ArrayCoordinateWithHeight neighbor;
	while ( !q.empty() || !pq.empty()) {
		
		if (q.empty()){
			c = pq.top();
			pq.pop();
		}else{
			c = q.front();
			q.pop();
		}

		for (uint d=0; d<directions.size(); d++) {
			neighbor = ArrayCoordinateWithHeight_init(c.row+directions[d].row,c.col+directions[d].col,0);		
			if (!DEM->check_within(c.row+directions[d].row, c.col+directions[d].col) || seen->get(neighbor.row,neighbor.col))
				continue;
			neighbor.h = DEM_filled_no_flat->get(neighbor.row,neighbor.col);

			seen->set(neighbor.row,neighbor.col,true);

			if (neighbor.h<=c.h) {
				DEM_filled_no_flat->set(neighbor.row,neighbor.col,DEM_filled_no_flat->get(c.row,c.col) + EPS);
				neighbor.h = DEM_filled_no_flat->get(neighbor.row,neighbor.col);
				q.push(neighbor);
			}
			else {
				pq.push(neighbor);
			}
		}
	}
	delete seen;
	return DEM_filled_no_flat;
}

// Find the lowest neighbor of a point given the cos of the latitude (for speed optimization)
int find_lowest_neighbor(int row, int col, Model<double>* DEM_filled_no_flat, double coslat)
{
	int result = 0;
	double min_drop = 0;
	double min_dist = 100000;
	for (uint d=0; d<directions.size(); d++) {
		int row_neighbor = row+directions[d].row;
		int col_neighbor = col+directions[d].col;
		if (DEM_filled_no_flat->check_within(row_neighbor, col_neighbor)) {
			double drop = DEM_filled_no_flat->get(row_neighbor,col_neighbor)-DEM_filled_no_flat->get(row,col);
			if(drop>=0)
				continue;
			double dist = find_distance(DEM_filled_no_flat->get_coordinate(row, col), DEM_filled_no_flat->get_coordinate(row_neighbor, col_neighbor), coslat);
			if (drop*min_dist < min_drop*dist) {
				min_drop = drop;
				min_dist = dist;
				result = d;
			}
		}
	}
	if(min_drop==0)
		if(display)
			printf("Alert: Minimum drop of 0 at %d %d\n", row, col);
	return result;
}

// Find the direction of flow for each square in a filled DEM
static Model<char>* flow_direction(Model<double>* DEM_filled_no_flat)
{
	Model<char>* flow_dirn = new Model<char>(DEM_filled_no_flat->nrows(), DEM_filled_no_flat->ncols(), MODEL_UNSET);
	flow_dirn->set_geodata(DEM_filled_no_flat->get_geodata());
	double coslat = COS(RADIANS(flow_dirn->get_origin().lat-(0.5+border/(double)(flow_dirn->nrows()-2*border))));
	for (int row=1; row<flow_dirn->nrows()-1;row++)
		for (int col=1; col<flow_dirn->ncols()-1; col++) 
			flow_dirn->set(row,col,find_lowest_neighbor(row, col, DEM_filled_no_flat, coslat));

	for (int row=0; row<flow_dirn->nrows()-1;row++) {
		flow_dirn->set(row,0,4);
		flow_dirn->set(flow_dirn->nrows()-row-1,flow_dirn->ncols()-1,0);
	}
	for (int col=0; col<flow_dirn->ncols()-1; col++) {
 		flow_dirn->set(0,col+1,6);
		flow_dirn->set(flow_dirn->nrows()-1,flow_dirn->ncols()-col-2,2);
	}
	flow_dirn->set(0,0,5);
	flow_dirn->set(0,flow_dirn->ncols()-1,7);
	flow_dirn->set(flow_dirn->nrows()-1,flow_dirn->ncols()-1,1);
	flow_dirn->set(flow_dirn->nrows()-1,0,3);
	return flow_dirn;
}

// Calculate flow accumulation given a DEM and the flow directions
static Model<int>* find_flow_accumulation(Model<char>* flow_directions, Model<double>* DEM_filled_no_flat)
{
	Model<int>* flow_accumulation = new Model<int>(DEM_filled_no_flat->nrows(), DEM_filled_no_flat->ncols(), MODEL_SET_ZERO);
	flow_accumulation->set_geodata(DEM_filled_no_flat->get_geodata());

	ArrayCoordinateWithHeight* to_check = new ArrayCoordinateWithHeight[flow_accumulation->nrows()*flow_accumulation->ncols()];

	for (int row=0; row<flow_accumulation->nrows(); row++)
		for (int col=0; col<flow_accumulation->ncols();col++) {
			to_check[DEM_filled_no_flat->ncols()*row+col].col = col;
			to_check[DEM_filled_no_flat->ncols()*row+col].row = row;
			to_check[DEM_filled_no_flat->ncols()*row+col].h = DEM_filled_no_flat->get(row,col);
		}

	// start at the highest point and distribute flow downstream
	sort(to_check, to_check+flow_accumulation->nrows()*flow_accumulation->ncols());

	for (int i=0; i<DEM_filled_no_flat->nrows()*DEM_filled_no_flat->ncols();i++) {
		ArrayCoordinateWithHeight p = to_check[i];
		ArrayCoordinateWithHeight downstream = ArrayCoordinateWithHeight_init(p.row+directions[flow_directions->get(p.row,p.col)].row, p.col+directions[flow_directions->get(p.row,p.col)].col,0);
		if (flow_accumulation->check_within( downstream.row, downstream.col) ){
			//printf("%4d %4d %3d %3d %2d %2d %1d\n", downstream.row,downstream.col,flow_accumulation->get(p.row,p.col), flow_accumulation->get(p.row,p.col)+1, directions[flow_directions->get(p.row,p.col)].row, directions[flow_directions->get(p.row,p.col)].col, flow_directions->get(p.row,p.col));
			flow_accumulation->set(downstream.row,downstream.col, flow_accumulation->get(p.row,p.col)+flow_accumulation->get(downstream.row,downstream.col)+1);
		}
	}
	delete to_check;
	return flow_accumulation;
}

// Find streams given the flow accumulation
static Model<bool>* find_streams(Model<int>* flow_accumulation)
{
	Model<bool>* streams = new Model<bool>(flow_accumulation->nrows(), flow_accumulation->ncols(), MODEL_SET_ZERO);
	streams->set_geodata(flow_accumulation->get_geodata());
	int stream_site_count=0;
	for (int row=0; row<flow_accumulation->nrows(); row++)
		for (int col=0; col<flow_accumulation->ncols();col++) 
			if(flow_accumulation->get(row,col) >= stream_threshold){
				streams->set(row,col,true);
				stream_site_count++;
			}
	if(display)
		printf("Number of stream sites = %d\n",  stream_site_count);
	return streams;
}


// Find dam sites to check given the streams, flow diresctions and DEM
static Model<bool>* find_pour_points(Model<bool>* streams, Model<char>* flow_directions, Model<short>* DEM_filled)
{
	Model<bool>* pour_points = new Model<bool>(streams->nrows(), streams->ncols(), MODEL_SET_ZERO);
	pour_points->set_geodata(streams->get_geodata());
	int pour_point_count=0;
	for (int row = border; row < border+pour_points->nrows()-2*border; row++)
		for (int col = border; col <  border+pour_points->ncols()-2*border; col++)
			if (streams->get(row,col)) {
				ArrayCoordinate downstream = ArrayCoordinate_init(row+directions[flow_directions->get(row,col)].row,col+directions[flow_directions->get(row,col)].col, GeographicCoordinate_init(0,0));
				if ( flow_directions->check_within(downstream.row, downstream.col) &&
					DEM_filled->get(row,col)-DEM_filled->get(row,col)%contour_height>DEM_filled->get(downstream.row,downstream.col)) {
						pour_points->set(row,col,true);
						pour_point_count++;
					}
			}
	if(display)
		printf("Number of dam sites = %d\n",  pour_point_count);
	return pour_points;
}

// Find details of possible reservoirs at pour_point
static RoughReservoir model_reservoir(ArrayCoordinate pour_point, Model<char>* flow_directions, Model<short>* DEM_filled, Model<bool>* filter,
				  Model<int>* modelling_array, int iterator)
{

	RoughReservoir reservoir = RoughReservoir_init(pour_point, (int)(DEM_filled->get(pour_point.row,pour_point.col)));

	double area_at_elevation[max_wall_height+1] = {0};
	double cumulative_area_at_elevation[max_wall_height+1] = {0};
	double volume_at_elevation[max_wall_height+1] = {0};
	double dam_length_at_elevation[max_wall_height+1] = {0};

	queue<ArrayCoordinate> q;
	q.push(pour_point);
	while (!q.empty()) {
		ArrayCoordinate p = q.front();
		q.pop();

		int elevation = (int)(DEM_filled->get(p.row,p.col));
		int elevation_above_pp = MAX(elevation - reservoir.elevation, 0);

		update_reservoir_boundary(reservoir.shape_bound, p, elevation_above_pp);

		if (filter->get(p.row,p.col))
			reservoir.max_dam_height = MIN(reservoir.max_dam_height,elevation_above_pp);

		area_at_elevation[elevation_above_pp+1] += find_area(p); 
		modelling_array->set(p.row,p.col,iterator);

		for (uint d=0; d<directions.size(); d++) {
			ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
			if (flow_directions->check_within(neighbor.row, neighbor.col) &&
			    flows_to(neighbor, p, flow_directions) &&
			    ((int)(DEM_filled->get(neighbor.row,neighbor.col)-DEM_filled->get(pour_point.row,pour_point.col)) < max_wall_height) ) {
				q.push(neighbor);
			}
		}

	}

	for (int ih=1; ih<max_wall_height+1;ih++) {
		cumulative_area_at_elevation[ih] = cumulative_area_at_elevation[ih-1] + area_at_elevation[ih];
		volume_at_elevation[ih] = volume_at_elevation[ih-1] + 0.01*cumulative_area_at_elevation[ih]; // area in ha, vol in GL
	}

	q.push(pour_point);
	while (!q.empty()) {
		ArrayCoordinate p = q.front();
		q.pop();
		int elevation = (int)(DEM_filled->get(p.row,p.col));
		int elevation_above_pp = MAX(elevation - reservoir.elevation,0);
		for (uint d=0; d<directions.size(); d++) {
			ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
			if (flow_directions->check_within(neighbor.row, neighbor.col)){
				if(flows_to(neighbor, p, flow_directions) &&
			    ((int)(DEM_filled->get(neighbor.row,neighbor.col)-DEM_filled->get(pour_point.row,pour_point.col)) < max_wall_height) ) {
					q.push(neighbor);
				}
				if ((directions[d].row * directions[d].col == 0) // coordinate orthogonal directions
				    && (modelling_array->get(neighbor.row,neighbor.col) < iterator ) ){
					dam_length_at_elevation[MIN(MAX(elevation_above_pp, (int)(DEM_filled->get(neighbor.row,neighbor.col)-reservoir.elevation)),max_wall_height)] +=find_orthogonal_nn_distance(p, neighbor);	//WE HAVE PROBLEM IF VALUE IS NEGATIVE???	
				}
			}
		}
	}

	for (uint ih =0 ; ih< dam_wall_heights.size(); ih++) {
		int height = dam_wall_heights[ih];
		reservoir.areas.push_back(cumulative_area_at_elevation[height]);
		reservoir.dam_volumes.push_back(0);
		for (int jh=0; jh < height; jh++)
			reservoir.dam_volumes[ih] += convert_to_dam_volume(height-jh, dam_length_at_elevation[jh]);
		reservoir.volumes.push_back(volume_at_elevation[height] + 0.5*reservoir.dam_volumes[ih]);
		reservoir.water_rocks.push_back(reservoir.volumes[ih]/reservoir.dam_volumes[ih]);
	}

	return reservoir;
}


static int model_reservoirs(GridSquare square_coordinate, Model<bool>* pour_points, Model<char>* flow_directions, Model<short>* DEM_filled, Model<int>* flow_accumulation, Model<bool>* filter)
{
	mkdir(convert_string(file_storage_location+"output/reservoirs"),0777);
	FILE *csv_file = fopen(convert_string(file_storage_location+"output/reservoirs/"+str(square_coordinate)+"_reservoirs.csv"), "w");
	if (!csv_file) {
	 	fprintf(stderr, "failed to open reservoir CSV file\n");
		exit(1);
    }
	write_rough_reservoir_csv_header(csv_file);

	mkdir(convert_string(file_storage_location+"processing_files/reservoirs"),0777);
	FILE *csv_data_file = fopen(convert_string(file_storage_location+"processing_files/reservoirs/"+str(square_coordinate)+"_reservoirs_data.csv"), "w");
	if (!csv_file) {
	 	fprintf(stderr, "failed to open reservoir CSV data file\n");
		exit(1);
    }
	write_rough_reservoir_data_header(csv_data_file);

	int i = 0;
	int count = 0;
	Model<int>* model = new Model<int>(pour_points->nrows(), pour_points->ncols(), MODEL_SET_ZERO);

	for (int row=border;row <border+DEM_filled->nrows()-2*border; row++)
		for (int col=border;col <border+DEM_filled->ncols()-2*border; col++) {
			if (!pour_points->get(row,col) || filter->get(row,col)) continue;
			ArrayCoordinate pour_point = {row, col , get_origin(square_coordinate, border)};
			i++;
			RoughReservoir reservoir = model_reservoir(pour_point, flow_directions, DEM_filled, filter, model, i);

			if ( max(reservoir.volumes)>=min_reservoir_volume &&
			     max(reservoir.water_rocks)>min_reservoir_water_rock &&
			     reservoir.max_dam_height>=min_max_dam_height) {
				reservoir.watershed_area = find_area(pour_point)*flow_accumulation->get(row,col);

				reservoir.identifier = str(square_coordinate)+"_RES"+str(i);

				write_rough_reservoir_csv(csv_file, reservoir);
				write_rough_reservoir_data(csv_data_file, reservoir);
				count++;
			}
		}
	fclose(csv_file);
	fclose(csv_data_file);
	return count;
}

int main(int nargs, char **argv)
{
	GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));
	if(nargs>3)
		display = atoi(argv[3]);

	printf("Screening started for %s\n",convert_string(str(square_coordinate)));

	GDALAllRegister();
	parse_variables(convert_string("storage_location"));
	parse_variables(convert_string(file_storage_location+"variables"));
	unsigned long start_usec = walltime_usec();
	unsigned long t_usec = start_usec;

	Model<short>* DEM = read_DEM_with_borders(square_coordinate, border);

	if (display) {
		printf("\nAfter border added:\n");
		DEM->print();
	}
	if(debug_output){
		mkdir(convert_string("debug"),0777);
		mkdir(convert_string("debug/input"),0777);
		DEM->write("debug/input/"+str(square_coordinate)+"_input.tif", GDT_Int16);
	}

	Model<char>* flow_directions;
	Model<bool>* pour_points;
	Model<int>* flow_accumulation;
	Model<short>* DEM_filled;
	Model<bool>* filter;

	if (debug) {
		filter = new Model<bool>(file_storage_location+"debug/filter/"+str(square_coordinate)+"_filter.tif", GDT_Byte);
		DEM_filled = new Model<short>(file_storage_location+"debug/DEM_filled/"+str(square_coordinate)+"_DEM_filled.tif", GDT_Int16);
		flow_directions = new Model<char>(file_storage_location+"debug/flow_directions/"+str(square_coordinate)+"_flow_directions.tif", GDT_Byte);
		flow_accumulation =  new Model<int>(file_storage_location+"debug/flow_accumulation/"+str(square_coordinate)+"_flow_accumulation.tif", GDT_Int32);
		pour_points = new Model<bool>(file_storage_location+"debug/pour_points/"+str(square_coordinate)+"_pour_points.tif", GDT_Byte);
	}else {
		t_usec = walltime_usec();
		filter = read_filter(DEM, filter_filenames);
		if (display) {
			printf("\nFilter:\n");
			filter->print();
			printf("Filter Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - t_usec) );
		}
		if(debug_output){
			mkdir(convert_string(file_storage_location+"debug/filter"),0777);
			filter->write(file_storage_location+"debug/filter/"+str(square_coordinate)+"_filter.tif", GDT_Byte);
		}

		t_usec = walltime_usec();
		Model<double>* DEM_filled_no_flat = fill(DEM);
		DEM_filled = new Model<short>(DEM->nrows(), DEM->ncols(), MODEL_SET_ZERO);
		DEM_filled->set_geodata(DEM->get_geodata());
		for(int row = 0; row<DEM->nrows();row++)
			for(int col = 0; col<DEM->ncols();col++)
				DEM_filled->set(row, col, convert_to_int(DEM_filled_no_flat->get(row, col)));
		if (display) {
			printf("\nFilled No Flats:\n");
			DEM_filled_no_flat->print();
			printf("Fill Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - t_usec) );
		}
		if(debug_output){
			mkdir(convert_string(file_storage_location+"debug/DEM_filled"),0777);
			DEM_filled->write(file_storage_location+"debug/DEM_filled/"+str(square_coordinate)+"_DEM_filled.tif", GDT_Int16);
			DEM_filled_no_flat->write(file_storage_location+"debug/DEM_filled/"+str(square_coordinate)+"_DEM_filled_no_flat.tif",GDT_Float64);
		}

		t_usec = walltime_usec();
		flow_directions = flow_direction(DEM_filled_no_flat);
		if (display) {
			printf("\nFlow Directions:\n");
			flow_directions->print();
			printf("Flow directions Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - t_usec) );
		}
		if(debug_output){
			mkdir(convert_string(file_storage_location+"debug/flow_directions"),0777);
			flow_directions->write(file_storage_location+"debug/flow_directions/"+str(square_coordinate)+"_flow_directions.tif",GDT_Byte);
		}
		mkdir(convert_string(file_storage_location+"processing_files/flow_directions"),0777);
		flow_directions->write(file_storage_location+"processing_files/flow_directions/"+str(square_coordinate)+"_flow_directions.tif",GDT_Byte);

		t_usec = walltime_usec();
		flow_accumulation = find_flow_accumulation(flow_directions, DEM_filled_no_flat);
		if (display) {
			printf("\nFlow Accumulation:\n");
			flow_accumulation->print();
			printf("Flow accumulation Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - t_usec) );
		}
		if(debug_output){
			mkdir(convert_string(file_storage_location+"debug/flow_accumulation"),0777);
			flow_accumulation->write(file_storage_location+"debug/flow_accumulation/"+str(square_coordinate)+"_flow_accumulation.tif", GDT_Int32);
		}
		delete DEM_filled_no_flat;
		
		Model<bool>* streams = find_streams(flow_accumulation);
		if (display) {
			printf("\nStreams (Greater than %d accumulation):\n", stream_threshold);
			streams->print();
		}
		if(debug_output){
			mkdir(convert_string(file_storage_location+"debug/streams"),0777);
			streams->write(file_storage_location+"debug/streams/"+str(square_coordinate)+"_streams.tif",GDT_Byte);
		}

		pour_points = find_pour_points(streams, flow_directions, DEM_filled);
		if (display) {
			printf("\nPour points (Streams every %dm):\n", contour_height);
			pour_points->print();
		}
		if(debug_output){
			mkdir(convert_string(file_storage_location+"debug/pour_points"),0777);
			pour_points->write(file_storage_location+"debug/pour_points/"+str(square_coordinate)+"_pour_points.tif",GDT_Byte);
		}
		delete streams;
	}

	t_usec = walltime_usec();
	int count = model_reservoirs(square_coordinate, pour_points, flow_directions, DEM_filled, flow_accumulation, filter);
	if(display)
		printf("Found %d reservoirs. Runtime: %.2f sec\n", count, 1.0e-6*(walltime_usec() - t_usec));
	printf(convert_string("Screening finished for "+str(square_coordinate)+". Runtime: %.2f sec\n"), 1.0e-6*(walltime_usec() - start_usec) );
}
