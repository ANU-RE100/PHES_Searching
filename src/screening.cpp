#include "phes_base.h"

bool debug_output = false;
bool debug = false;
int display = false;

bool check_intersection(double lat, GeographicCoordinate line[2]){
    if ((line[0].lat < lat && line[1].lat>=lat) || (line[0].lat >= lat && line[1].lat < lat))
        return true;
    return false;
}

// find_intersection finds the longitude of the intersection of a latitude (lat) and a line, assuming they intersect
double find_intersection(double lat, GeographicCoordinate line[2]){
    return (line[0].lon+(lat-line[0].lat)/(line[1].lat-line[0].lat)*(line[1].lon-line[0].lon));
}

// find_polygon_intersections returns an array containing the longitude of all line. Assumes last coordinate is same as first
vector<double> find_polygon_intersections(int row, vector<GeographicCoordinate> &polygon, GeographicCoordinate origin){
    vector<double> to_return;
    double lat = convert_coordinates(ArrayCoordinate_init(row, 0, origin)).lat;
    for(uint i = 0; i<polygon.size()-1; i++){
        GeographicCoordinate line[2] = {polygon[i], polygon[(i+1)]};
        if(check_intersection(lat, line)){
            to_return.push_back(find_intersection(lat, line));
        }
    }
    sort(to_return.begin(), to_return.end());
    return to_return;
}

static void read_shp_filter(string filename, GeographicCoordinate origin, int DEM_shape[2], Model_int16* filter){
	SHPHandle SHP = SHPOpen(convert_string(filename), "rb" );
	if(SHP != NULL ){
    	int	nEntities, nShapeType;
    	double 	adfMinBound[4], adfMaxBound[4];
    	vector<vector<GeographicCoordinate> > relevant_polygons;
    	SHPGetInfo(SHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );
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
	            to_keep = (to_keep || check_within(convert_coordinates(temp_point, origin),DEM_shape));
	            temp_poly.push_back(temp_point);
	        }
	        if(to_keep)
	        	relevant_polygons.push_back(temp_poly);
	        SHPDestroyObject( shape );
	    }
	    if(display)
	    	printf("%d Polygons imported from %s\n", (int)relevant_polygons.size(), convert_string(filename));
	    for(uint i = 0; i<relevant_polygons.size(); i++){
	    	for(int row =0; row<DEM_shape[0]; row++){
                vector<double> polygon_intersections = find_polygon_intersections(row, relevant_polygons[i], origin);
                for(uint j = 0; j<polygon_intersections.size();j++){
                	polygon_intersections[j] = (convert_coordinates(GeographicCoordinate_init(0, polygon_intersections[j]),origin).col);
                }
                for(uint j = 0; j<polygon_intersections.size()/2;j++){
                    for(int col=polygon_intersections[2*j];col<polygon_intersections[2*j+1];col++){
                    	ArrayCoordinate coordinate = ArrayCoordinate_init(row, col, origin);
                        if(check_within(coordinate,DEM_shape)){
                        	filter->d[row][col]=true;
                        }
                    }
                }
            }
	    }
    }
    SHPClose( SHP );
}

void read_tif_filter(string filename, GeographicCoordinate origin, Model_int16* filter, int value_to_filter){
	try{
		double geotransform[6];
		Model_int8* tif_filter = TIFF_Read_int8(convert_string(filename), geotransform, NULL);
		GeographicCoordinate tif_origin = {geotransform[3], geotransform[0]};
		for(int row = 0; row<filter->shape[0]; row++){
			for(int col = 0; col<filter->shape[1]; col++){
				GeographicCoordinate point = convert_coordinates(ArrayCoordinate_init(row, col, origin));
				ArrayCoordinate tif_point = convert_coordinates(point, tif_origin, geotransform[5], geotransform[1]);
				if(check_within(tif_point, tif_filter->shape) && tif_filter->d[tif_point.row][tif_point.col]==value_to_filter)
					filter->d[row][col] = true;
			}
		}
	}catch(int e){
		if(display)
			printf("Problem with %s\n", convert_string(filename));
	}
}

string find_world_urban_filename(GeographicCoordinate point){
	char clat = 65+floor((point.lat+90)/8+1);
	if(clat>=73)
		clat++;
	if(clat>=79)
		clat++;
	int nlon = floor((point.lon+180)/6+1);

	char strlon[3];
	snprintf(strlon, 3, "%02d", nlon);
	return "input/WORLD_URBAN/"+string(strlon)+string(1, clat)+"_hbase_human_built_up_and_settlement_extent_geographic_30m.tif";
}

vector<GeographicCoordinate> get_corners(GeographicCoordinate origin, int* DEM_shape){
	vector<GeographicCoordinate> to_return = {
		origin, 
		convert_coordinates(ArrayCoordinate_init(0, DEM_shape[1], origin)),
		convert_coordinates(ArrayCoordinate_init(DEM_shape[0], 0, origin)),
		convert_coordinates(ArrayCoordinate_init(DEM_shape[0], DEM_shape[1], origin))
	};
	return to_return;
}

static Model_int16 *read_filter(vector<string> filenames, GeographicCoordinate origin, int DEM_shape[2])
{
	vector<string> filter_files;
	Model_int16 *filter = Model_int16_create(DEM_shape, MODEL_SET_ZERO);
	for(string filename:filenames){
		if(filename=="use_world_urban"){
			if(display)
				printf("Using world urban data\n");
			vector<string> required;
			for(GeographicCoordinate corner: get_corners(origin, DEM_shape)){
				string urban_filename = find_world_urban_filename(corner);
				if(find(required.begin(), required.end(), urban_filename)==required.end()){
					required.push_back(urban_filename);
					read_tif_filter(urban_filename, origin, filter, -55);
				}
			}
		}else{
			read_shp_filter(filename, origin, DEM_shape, filter);
		}		
	}
	return filter;
}


static Model_double *fill(Model_int16 *DEM)
{
	Model_double *DEM_filled = Model_double_create(DEM->shape, MODEL_UNSET);

	queue<ArrayCoordinateWithHeight> q;
	priority_queue<ArrayCoordinateWithHeight> pq;

	Model_int8 *mseen = Model_int8_create(DEM->shape, MODEL_SET_ZERO);
	char **seen = mseen->d;

	for (int row=0; row<DEM->shape[0]; row++){
		for (int col=0; col<DEM->shape[1];col++) {
			DEM_filled->d[row][col] = (double)DEM->d[row][col];
		}
	}

	
	for (int row=0; row<DEM->shape[0]-1;row++) {
		pq.push(ArrayCoordinateWithHeight_init(row+1, DEM->shape[1]-1, (double)DEM->d[row+1][DEM->shape[1]-1]));
		seen[row+1][DEM->shape[1]-1]=true;
		pq.push(ArrayCoordinateWithHeight_init(row, 0, (double)DEM->d[row][0]));
		seen[row][0]=true;
	}

	for (int col=0; col<DEM->shape[1]-1;col++) {
		pq.push(ArrayCoordinateWithHeight_init(DEM->shape[0]-1, col, (double)DEM->d[DEM->shape[0]-1][col]));
		seen[DEM->shape[0]-1][col] = true;
		pq.push(ArrayCoordinateWithHeight_init(0, col+1, (double)DEM->d[0][col+1]));
		seen[0][col+1] = true;
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
			if (!check_within(neighbor, DEM->shape) || seen[neighbor.row][neighbor.col]){
				continue;
			}
			neighbor.h = DEM->d[neighbor.row][neighbor.col];

			seen[neighbor.row][neighbor.col] = true;

			if (neighbor.h<=c.h) {
				DEM_filled->d[neighbor.row][neighbor.col] = DEM_filled->d[c.row][c.col] + EPS; //EPS_lift;
				neighbor.h = DEM_filled->d[c.row][c.col] + EPS;
				q.push(neighbor);
			}
			else {
				pq.push(neighbor);
			}
		}
	}
	Model_int8_free(mseen);
	return DEM_filled;
}

// Find the direction of flow for each square in a filled DEM
static Model_int16 *flow_direction(Model_double *DEM_filled, GeographicCoordinate origin)
{
	Model_int16 *flow_dirn = Model_int16_create(DEM_filled->shape, MODEL_UNSET);
	ArrayCoordinate c;
	c.origin = origin;
	double coslat = COS(RADIANS(origin.lat-(0.5+border/(double)(DEM_filled->shape[0]-2*border))));
	for (int row=0; row<DEM_filled->shape[0];row++)
		for (int col=0; col<DEM_filled->shape[1]; col++) {
			c.row = row;
			c.col = col;
			flow_dirn->d[row][col] = find_lowest_neighbor(c, DEM_filled, coslat);
		}

	for (int row=0; row<DEM_filled->shape[0]-1;row++) {
		flow_dirn->d[row][0]=16; // -dx
		flow_dirn->d[DEM_filled->shape[0]-row-1][DEM_filled->shape[1]-1]=1; // +dx
	}
	for (int col=0; col<DEM_filled->shape[1]-1; col++) {
 		flow_dirn->d[0][col+1]=64; // +dy
		flow_dirn->d[DEM_filled->shape[0]-1][DEM_filled->shape[1]-col-2]=4; // -dy
	}
	flow_dirn->d[0][0]=32;
	flow_dirn->d[0][DEM_filled->shape[1]-1]=128;
	flow_dirn->d[DEM_filled->shape[0]-1][DEM_filled->shape[1]-1]=2;
	flow_dirn->d[DEM_filled->shape[0]-1][0]=8;

	return flow_dirn;
}

// Calculate flow accumulation given a DEM and the flow directions
static Model_int32 *find_flow_accumulation(Model_int16 *flow_directions, Model_double *DEM_filled)
{
	Model_int32 *flow_accumulation = Model_int32_create(DEM_filled->shape, MODEL_UNSET);

	ArrayCoordinateWithHeight *to_check = (ArrayCoordinateWithHeight*)malloc(DEM_filled->shape[0]*DEM_filled->shape[1]*sizeof(ArrayCoordinateWithHeight));

	for (int row=0; row<DEM_filled->shape[0]; row++){
		for (int col=0; col<DEM_filled->shape[1];col++) {
			flow_accumulation->d[row][col] = 0;
			to_check[DEM_filled->shape[1]*row+col].col = col;
			to_check[DEM_filled->shape[1]*row+col].row = row;
			to_check[DEM_filled->shape[1]*row+col].h = DEM_filled->d[row][col];
		}
	}

	// start at the highest point and distribute flow downstream
	sort(to_check, to_check+DEM_filled->shape[0]*DEM_filled->shape[1]);

	for (int i=0; i<DEM_filled->shape[0]*DEM_filled->shape[1];i++) {
		ArrayCoordinateWithHeight p = to_check[i];
		ArrayCoordinateWithHeight downstream = ArrayCoordinateWithHeight_init(p.row+directions[flow_directions->d[p.row][p.col]].row, p.col+directions[flow_directions->d[p.row][p.col]].col,0);
		if (check_within( downstream, DEM_filled->shape) )
			flow_accumulation->d[downstream.row][downstream.col] += flow_accumulation->d[p.row][p.col]+1;
	}

	free(to_check);

	return flow_accumulation;
}


// Find streams given the flow accumulation
static Model_int8 *find_streams(Model_int32 *flow_accumulation)
{
	Model_int8 *streams = Model_int8_create(flow_accumulation->shape, MODEL_UNSET);

	int stream_site_count=0;
	for (int row=0; row<flow_accumulation->shape[1]; row++){
		for (int col=0; col<flow_accumulation->shape[0];col++) {
			streams->d[row][col] = (flow_accumulation->d[row][col] >= stream_threshold);
			if (streams->d[row][col]) {
				stream_site_count++;
			}
		}
	}

	if(display)
		printf("Number of stream sites = %d\n",  stream_site_count);

	return streams;
}


// Find dam sites to check given the streams, flow diresctions and DEM
static Model_int8 *find_pour_points(Model_int8 *streams, Model_int16 *flow_directions, Model_double *DEM_filled)
{
	Model_int8 *pour_points = Model_int8_create(flow_directions->shape, MODEL_SET_ZERO);

	int pour_point_count=0;
	for (int row = border; row < border+streams->shape[0]-2*border; row++)
		for (int col = border; col <  border+streams->shape[1]-2*border; col++) {
			if (streams->d[row][col]) {
				ArrayCoordinate downstream = ArrayCoordinate_init(row+directions[flow_directions->d[row][col]].row,col+directions[flow_directions->d[row][col]].col, GeographicCoordinate_init(0,0));
				if ( check_within(downstream, flow_directions->shape)) {
					double new_height = DEM_filled->d[downstream.row][downstream.col];
					double cur_height = DEM_filled->d[row][col];
					if ( contour_height*floorf(cur_height/contour_height) > new_height) {
						pour_points->d[row][col] = true;
						pour_point_count++;
					}
				}
			}
		}
	if(display)
		printf("Number of dam sites = %d\n",  pour_point_count);

	return pour_points;
}

// Find details of possible reservoirs at pour_point
static RoughReservoir model_reservoir(ArrayCoordinate pour_point, Model_int16 *flow_directions, Model_int16 *DEM, Model_int16 *filter,
				  Model_int32 *modelling_array, int iterator)
{

	RoughReservoir reservoir = RoughReservoir_init(pour_point, (int)(DEM->d[pour_point.row][pour_point.col]));

	double area_at_elevation[max_wall_height+1] = {0};
	double cumulative_area_at_elevation[max_wall_height+1] = {0};
	double volume_at_elevation[max_wall_height+1] = {0};
	double dam_length_at_elevation[max_wall_height+1] = {0};

	queue<ArrayCoordinate> q;
	q.push(pour_point);
	while (!q.empty()) {
		ArrayCoordinate p = q.front();
		q.pop();

		int elevation = (int)(DEM->d[p.row][p.col]);
		int elevation_above_pp = MAX(elevation - reservoir.elevation, 0);

		update_reservoir_boundary(reservoir.shape_bound, p, elevation_above_pp);

		if (filter->d[p.row][p.col])
			reservoir.max_dam_height = MIN(reservoir.max_dam_height,elevation_above_pp);

		area_at_elevation[elevation_above_pp+1] += find_area(p); 
		modelling_array->d[p.row][p.col] = iterator;

		for (uint d=0; d<directions.size(); d++) {
			ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
			if (check_within(neighbor, flow_directions->shape) &&
			    flows_to(neighbor, p, flow_directions) &&
			    ((int)(DEM->d[neighbor.row][neighbor.col]-DEM->d[pour_point.row][pour_point.col]) < max_wall_height) ) {
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
		int elevation = (int)(DEM->d[p.row][p.col]);
		int elevation_above_pp = MAX(elevation - reservoir.elevation,0);
		for (uint d=0; d<directions.size(); d++) {
			ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
			if (check_within(neighbor, flow_directions->shape)){
				if(flows_to(neighbor, p, flow_directions) &&
			    ((int)(DEM->d[neighbor.row][neighbor.col]-DEM->d[pour_point.row][pour_point.col]) < max_wall_height) ) {
					q.push(neighbor);
				}
				if ((directions[d].row * directions[d].col == 0) // coordinate orthogonal directions
				    && (modelling_array->d[neighbor.row][neighbor.col] < iterator ) ){
					dam_length_at_elevation[MIN(MAX(elevation_above_pp, (int)(DEM->d[neighbor.row][neighbor.col]-reservoir.elevation)),max_wall_height)] +=find_orthogonal_nn_distance(p, neighbor);	//WE HAVE PROBLEM IF VALUE IS NEGATIVE???	
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


static int model_reservoirs(GridSquare square_coordinate, Model_int8 *pour_points, Model_int16 *flow_directions, Model_int16 *DEM, Model_int32 *flow_accumulation, Model_int16 *filter)
{
	mkdir("output/reservoirs",0777);
	FILE *csv_file = fopen(convert_string("output/reservoirs/"+str(square_coordinate)+"_reservoirs.csv"), "w");
	if (!csv_file) {
	 	fprintf(stderr, "failed to open reservoir CSV file\n");
		exit(1);
    }
	write_rough_reservoir_csv_header(csv_file);

	mkdir("processing_files/reservoirs",0777);
	FILE *csv_data_file = fopen(convert_string("processing_files/reservoirs/"+str(square_coordinate)+"_reservoirs_data.csv"), "w");
	if (!csv_file) {
	 	fprintf(stderr, "failed to open reservoir CSV data file\n");
		exit(1);
    }
	write_rough_reservoir_data_header(csv_data_file);

	int i = 0;
	int count = 0;
	Model_int32 *model = Model_int32_create(flow_directions->shape, MODEL_SET_ZERO);

	for (int row=border;row <border+pour_points->shape[0]-2*border; row++)
		for (int col=border;col <border+pour_points->shape[1]-2*border; col++) {

			if (!pour_points->d[row][col]) continue;
			ArrayCoordinate pour_point = {row, col , get_origin(square_coordinate, border)};
			i++;
			RoughReservoir reservoir = model_reservoir(pour_point, flow_directions, DEM, filter, model, i);

			if ( max(reservoir.volumes)>=min_reservoir_volume &&
			     max(reservoir.water_rocks)>min_reservoir_water_rock &&
			     reservoir.max_dam_height>=min_max_dam_height) {
				reservoir.watershed_area = find_area(pour_point)*flow_accumulation->d[row][col];

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

	// Point origin;
	GeographicCoordinate origin = get_origin(square_coordinate, 600);

	TIFF_IO_init();
	parse_variables(convert_string("variables"));
	unsigned long start_usec = walltime_usec();
	unsigned long t_usec = start_usec;

	char *geoprojection;
	double geotransform[6];
	Model_int16 *DEM = read_DEM_with_borders(square_coordinate);
	Model_int16 *DEM_temp = TIFF_Read_int16(convert_string("input/"+str(square_coordinate)+"_1arc_v3.tif"), geotransform, &geoprojection);
	Model_int16_free(DEM_temp);

	geotransform[0] = origin.lon;
	geotransform[3] = origin.lat;

	if (display) {
		printf("\nAfter border added:\n");
		Model_int16_print(DEM);
	}
	if(debug_output){
		mkdir("debug",0777);
		mkdir("debug/input",0777);
		TIFF_Write_int16(convert_string("debug/input/"+str(square_coordinate)+"_input.tif"), geotransform, geoprojection, DEM);
	}

	

	Model_int16 *flow_directions;
	Model_int16 *flow_directions_esri;
	Model_int8 *pour_points;
	Model_int32 *flow_accumulation;
	Model_double *DEM_filled;
	Model_int16 *DEM_filled_int;
	Model_int16 *filter;


	if (debug) {
		// Read in preprocessed files to save time during debug
		filter = TIFF_Read_int16(convert_string("debug/filter/"+str(square_coordinate)+"_filter.tif"), NULL, NULL);
		DEM_filled = TIFF_Read_double(convert_string("debug/DEM_filled/"+str(square_coordinate)+"_DEM_filled.tif"), NULL, NULL);
		DEM_filled_int = TIFF_Read_int16(convert_string("debug/DEM_filled/"+str(square_coordinate)+"_DEM_filled_int.tif"), NULL, NULL);
		flow_directions = TIFF_Read_int16(convert_string("debug/flow_directions/"+str(square_coordinate)+"_flow_directions.tif"), NULL, NULL);
		flow_accumulation =  TIFF_Read_int32(convert_string("debug/flow_accumulation/"+str(square_coordinate)+"_flow_accumulation.tif"), NULL, NULL);
		pour_points = TIFF_Read_int8(convert_string("debug/pour_points/"+str(square_coordinate)+"_pour_points.tif"), NULL, NULL);
	}else {
		// If not using preprocessed files, process files
		t_usec = walltime_usec();
		filter = read_filter(filter_filenames, origin, DEM->shape);
		if (display) {
			printf("\nFilter:\n");
			Model_int16_print(filter);
			printf("Filter Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - t_usec) );
		}
		if(debug_output){
			mkdir("debug/filter",0777);
			TIFF_Write_int16(convert_string("debug/filter/"+str(square_coordinate)+"_filter.tif"), geotransform, geoprojection, filter);
		}

		t_usec = walltime_usec();
		DEM_filled = fill(DEM);
		if (display) {
			printf("\nFilled No Flats:\n");
			Model_double_print(DEM_filled);
			printf("Fill Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - t_usec) );
		}
		DEM_filled_int = Model_double_to_int16(DEM_filled);
		if(debug_output){
			mkdir("debug/DEM_filled",0777);
			TIFF_Write_double(convert_string("debug/DEM_filled/"+str(square_coordinate)+"_DEM_filled.tif"), geotransform, geoprojection, DEM_filled);
			TIFF_Write_int16(convert_string("debug/DEM_filled/"+str(square_coordinate)+"_DEM_filled_int.tif"), geotransform, geoprojection, DEM_filled_int);
		}

		t_usec = walltime_usec();
		flow_directions = flow_direction(DEM_filled, origin);
		if (display) {
			printf("\nFlow Directions:\n");
			Model_int16_print(flow_directions);
			flow_directions_esri = Model_int16_create(flow_directions->shape, MODEL_UNSET);
			for(int row=0 ; row < flow_directions->shape[0] ; row++)
				for(int col=0 ; col < flow_directions->shape[1] ; col++)
					flow_directions_esri->d[row][col] = directions[flow_directions->d[row][col]].val;
			printf("\nFlow Directions ESRI:\n");
			Model_int16_print(flow_directions_esri);
			printf("Flow directions Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - t_usec) );
		}
		if(debug_output){
			mkdir("debug/flow_directions",0777);
			TIFF_Write_int16(convert_string("debug/flow_directions/"+str(square_coordinate)+"_flow_directions.tif"), geotransform, geoprojection, flow_directions);
		}
		mkdir("processing_files/flow_directions",0777);
		TIFF_Write_int16(convert_string("processing_files/flow_directions/"+str(square_coordinate)+"_flow_directions.tif"), geotransform, geoprojection, flow_directions);

		t_usec = walltime_usec();
		flow_accumulation = find_flow_accumulation(flow_directions, DEM_filled);
		if (display) {
			printf("\nFlow Accumulation:\n");
			Model_int32_print(flow_accumulation);
			printf("Flow accumulation Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - t_usec) );
		}
		if(debug_output){
			mkdir("debug/flow_accumulation",0777);
			TIFF_Write_int32(convert_string("debug/flow_accumulation/"+str(square_coordinate)+"_flow_accumulation.tif"), geotransform, geoprojection, flow_accumulation);
		}
		

		Model_int8 *streams = find_streams(flow_accumulation);
		if (display) {
			printf("\nStreams (Greater than %d accumulation):\n", stream_threshold);
			Model_int8_print(streams);
		}
		if(debug_output){
			mkdir("debug/streams",0777);
			TIFF_Write_int8(convert_string("debug/streams/"+str(square_coordinate)+"_streams.tif"), geotransform, geoprojection, streams);
		}

		pour_points = find_pour_points(streams, flow_directions, DEM_filled);
		Model_int8_free(streams);
		if (display) {
			printf("\nPour points (Streams every %.1fm):\n", contour_height);
			Model_int8_print(pour_points);
		}
		if(debug_output){
			mkdir("debug/pour_points",0777);
			TIFF_Write_int8(convert_string("debug/pour_points/"+str(square_coordinate)+"_pour_points.tif"), geotransform, geoprojection, pour_points);
		}
	}

	t_usec = walltime_usec();
	int count = model_reservoirs(square_coordinate, pour_points, flow_directions, DEM_filled_int, flow_accumulation, filter);
	if(display)
		printf("Found %d reservoirs. Runtime: %.2f sec\n", count, 1.0e-6*(walltime_usec() - t_usec));
	printf(convert_string("Screening finished for "+str(square_coordinate)+". Runtime: %.2f sec\n"), 1.0e-6*(walltime_usec() - start_usec) );
}
