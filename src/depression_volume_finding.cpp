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



int main(int nargs, char **argv)
{
	GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));
	display = true;

	printf("Volume finding started for %s\n",convert_string(str(square_coordinate)));

	GDALAllRegister();
	parse_variables(convert_string(file_storage_location+"variables"));
	unsigned long start_usec = walltime_usec();

	Model<short>* DEM = read_DEM_with_borders(square_coordinate, border);
	Model<bool>* extent = new Model<bool>(DEM->nrows(), DEM->ncols(), MODEL_SET_ZERO);
	extent->set_geodata(DEM->get_geodata());
	string rs(argv[3]);
	read_shp_filter(rs, extent);

	int min_elevation = 100000000;

	for(int row = 0; row<extent->nrows(); row++)
		for(int col = 0; col<extent->ncols(); col++){
			if(extent->get(row, col))
				min_elevation = MIN(DEM->get(row, col), min_elevation);
		}

	double area_at_elevation[1001] = {0};
	double volume_at_elevation[1001] = {0};
	double cumulative_area_at_elevation[1001] = {0};
	
	for(int row = 0; row<extent->nrows(); row++)
		for(int col = 0; col<extent->ncols(); col++)
			if(extent->get(row, col)){
				int elevation_above_pp = MAX(DEM->get(row,col) - min_elevation, 0);
				area_at_elevation[elevation_above_pp+1] += find_area(ArrayCoordinate_init(row, col, DEM->get_origin())); 
			}

	for (int ih=1; ih<200;ih++) {
 		cumulative_area_at_elevation[ih] = cumulative_area_at_elevation[ih-1] + area_at_elevation[ih];
 		volume_at_elevation[ih] = volume_at_elevation[ih-1] + 0.01*cumulative_area_at_elevation[ih]; // area in ha, vol in GL
 		printf("%d %d %f   ", ih, min_elevation+ih, volume_at_elevation[ih]);
	}

	
	printf(convert_string("Volume finding finished for "+str(square_coordinate)+". Runtime: %.2f sec\n"), 1.0e-6*(walltime_usec() - start_usec) );
}
