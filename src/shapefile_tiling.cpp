#include "phes_base.h"

int display = false;


// Reads a list of cells to process from the tasks_file (Eg. 148 -36)
vector<GridSquare> read_tasklist(char *tasks_file)
{
	FILE *fd = fopen(tasks_file, "r");
	if (!fd)  {
		fprintf(stderr, "failed to open task file %s: %s\n", tasks_file, strerror(errno));
		exit(1);
	}
	vector<GridSquare> tasklist;
	int lon, lat, rc=0;
	while (rc != EOF) {	
		rc = fscanf(fd, "%d %d", &lon, &lat);
		tasklist.push_back(GridSquare_init(lat, lon));
	}
	printf("read %zu tasks\n", tasklist.size());
	fclose(fd);
	return tasklist;
}

int main()
{
	printf("Tiling started\n");

	parse_variables(convert_string(file_storage_location+"variables"));
	unsigned long start_usec = walltime_usec();

	vector<GridSquare> tasklist = read_tasklist(convert_string(tasks_file));

	Model<vector<int>>* relevant_polygons = new Model<vector<int>>(180, 360);
	Model<bool>* to_keep = new Model<bool>(180, 360, MODEL_SET_ZERO);
	vector<vector<GeographicCoordinate>> polygons;
	int iterator = 0;

	for(string filename:filter_filenames_to_tile){
		SHPHandle SHP = SHPOpen(convert_string(file_storage_location+filename), "rb" );
		if(SHP != NULL ){
	    	int	nEntities;
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

		        
		        for(int j = 0, iPart = 1; j < shape->nVertices; j++ )
		        {
		            if( iPart < shape->nParts && shape->panPartStart[iPart] == j )
		            {
		            	
		            	for(int lat = -90; lat<90; lat++)
		            		for(int lon = -180; lon<180; lon++)
		            			if(to_keep->get(lat+90, lon+180)){
		            				relevant_polygons->get_pointer(lat+90, lon+180)->push_back(iterator);
		            				to_keep->set(lat+90, lon+180,false);
		            			}
		            	polygons.push_back(temp_poly);
		            	iterator++;
		            	temp_poly.clear();
		                iPart++;
		            }
		            GeographicCoordinate temp_point = GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
		            to_keep->set((int)(temp_point.lat+90), (int)(temp_point.lon+180), true);
		            temp_poly.push_back(temp_point);
		        }
		        for(int lat = -90; lat<90; lat++)
            		for(int lon = -180; lon<180; lon++)
            			if(to_keep->get(lat+90, lon+180)){
             				relevant_polygons->get_pointer(lat+90, lon+180)->push_back(iterator);
            				to_keep->set(lat+90, lon+180,false);
            			}
            	polygons.push_back(temp_poly);
		        iterator++;
		        SHPDestroyObject( shape );
		    }
		    printf("%d Polygons imported from %s\n", iterator, convert_string(filename));
	    }
	    SHPClose(SHP);
	}

	mkdir(convert_string(file_storage_location+"input/shapefile_tiles"),0777);
	int z = 0;
	for(GridSquare gs:tasklist){
		z++;
		// if(z<44557)
		// 	continue;
		int lat = gs.lat;
		int lon = gs.lon;
		SHPHandle SHP = SHPCreate(convert_string(file_storage_location+"input/shapefile_tiles/"+str(gs)+"_shapefile_tile.shp"), SHPT_POLYGON);
	    if( SHP == NULL ){
			printf( "Unable to create:%s\n", convert_string(file_storage_location+"input/shapefile_tiles/"+str(gs)+"_shapefile_tile.shp"));
			exit(1);
	    }

	    SHPObject* psObject;
	    for(uint ipolygon = 0; ipolygon<relevant_polygons->get(lat+90,lon+180).size(); ipolygon++){
	    	int nVertices = polygons[relevant_polygons->get(lat+90,lon+180)[ipolygon]].size();
	    	int panParts[1] = {0};

	    	double*	padfX = new double [nVertices];
	    	double* padfY = new double [nVertices];
	    	double *padfZ = NULL, *padfM = NULL;

	    	for(int i = 0; i<nVertices; i++){
	    		GeographicCoordinate temp = polygons[relevant_polygons->get(lat+90,lon+180)[ipolygon]][i];
	    		padfX[i] = temp.lon;
	    		padfY[i] = temp.lat;
	    	}
	    	psObject = SHPCreateObject( SHPT_POLYGON, -1, 0, panParts, NULL,
                            nVertices, padfX, padfY, padfZ, padfM );
	    	
		    SHPWriteObject( SHP, -1, psObject );
		    SHPDestroyObject( psObject );
		    delete padfY;
		    delete padfX;
	    }
	    SHPClose(SHP);
	}
			
	printf("Tiling finished. Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - start_usec) );
}

