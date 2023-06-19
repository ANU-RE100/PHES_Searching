/*
Data for mining tenaments is obtained from the following 2 sources:
    1. Global-scale mining polygons v2 (https://doi.pangaea.de/10.1594/PANGAEA.942325)
    2. Open Street Maps landuse=quarry, industrial=mine, historical=mine (./tools/osm_mine_downloader.py)

The OSM shapefiles can be merged into a single GPKG using GDAL. An example shell script is given in ./tools/merge_shp_to_gpkg.sh

The gpkg_tiling program will take merged_gpkg and split it into Shapefiles containing polygons that overlap with each
1 degree latitude by 1 degree longitude gridsquare. This will allow the brownfield screening to only access polygons
that are relevant to each specific gridsquare, dramatically reducing PHES Searching time.
*/

#include "phes_base.h"
#include "gdal/ogrsf_frmts.h"
#include <fstream>

bool canReadFile(const std::string& filePath) {
    std::ifstream fileStream(filePath);
    return fileStream.good();
}

// Reads a list of cells to process from the tasks_file (Eg. 148 -36)
vector<GridSquare> read_tasklist(char *tasks_file)
{
	ifstream fd(tasks_file);
	if (!fd)  {
		fprintf(stderr, "failed to open task file %s: %s\n", tasks_file, strerror(errno));
		exit(1);
	}
	vector<GridSquare> tasklist;
	string line;
	while(getline(fd, line)){
		line.erase(remove(line.begin(), line.end(), '\n'), line.end());
		line.erase(remove(line.begin(), line.end(), '\r'), line.end());
		istringstream is(line);
	    int lon, lat;
	    if(is>>lon)
	    	if(is>>lat)
	    		tasklist.push_back(GridSquare_init(lat, lon));
	}
	printf("read %zu tasks\n", tasklist.size());
	fd.close();
	return tasklist;
}

int main()
{
	printf("Tiling started\n");

    GDALAllRegister();

	parse_variables(convert_string("storage_location"));
	parse_variables(convert_string(file_storage_location+"variables"));
	unsigned long start_usec = walltime_usec();

	vector<GridSquare> tasklist = read_tasklist(convert_string(file_storage_location + tasks_file));

    GDALDataset *poDS_g;

    poDS_g = (GDALDataset*) GDALOpenEx( gpkg_path.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );
    
    if( poDS_g == NULL ){
        printf("Attempted to open gpkg_path file at: %s\n", gpkg_path.c_str());
        if (canReadFile(gpkg_path))
            printf("gpkg_file exists and user has permission to read the file.\n");
        else
            printf("Either the gpkg_file does not exist or user does not have permission to read the file.\n");
        printf( "Could not open gpkg_path file. GDAL Error: %s\n", CPLGetLastErrorMsg() );
        exit( 1 );
    }

    OGRLayer  *poLayer_g;
    poLayer_g = poDS_g->GetLayerByName(gpkg_layer.c_str());

    if (poLayer_g == nullptr) {
        printf("gpkg_layer not found in dataset.\n");
        exit(1);
    }
    
    std::string lat_prefix;
    std::string lon_prefix;

    for(GridSquare gs:tasklist){
		int lat = gs.lat;
		int lon = gs.lon;

        OGRPolygon square;
        OGRLinearRing ring;
        ring.addPoint(lon, lat);
        ring.addPoint(lon, lat + 1);
        ring.addPoint(lon + 1, lat + 1);
        ring.addPoint(lon + 1, lat);
        ring.addPoint(lon, lat); // Repeat the first point to close the ring
        square.addRing(&ring);

        poLayer_g->SetSpatialFilter(&square);

        const char *pszDriverName = "ESRI Shapefile";
        GDALDriver *poDriver;
        poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);

        GDALDataset *poDS_s;
        string filename = mining_tenament_shp + str(search_config.grid_square) + ".shp";        

        printf("Started tiling grid square %i %i with feature count %i...\n", search_config.grid_square.lat, search_config.grid_square.lon, (int) poLayer_g->GetFeatureCount());
        
        poDS_s = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
        if( poDS_s == NULL ){
            printf( "Creation of output file failed.\n" );
            exit( 1 );
        }

        OGRLayer *poLayer_s;
        poLayer_s = poDS_s->CreateLayer("result", poLayer_g->GetSpatialRef(), wkbPolygon, NULL );
        if( poLayer_s == NULL ){
            printf( "Layer creation failed.\n" );
            exit( 1 );
        }

        OGRFeature *poFeature;
        while ((poFeature = poLayer_g->GetNextFeature()) != NULL){
            OGRGeometry *poGeometry = poFeature->GetGeometryRef();
        
            // Check if the geometry is either Polygon or MultiPolygon
            if (poGeometry != NULL && 
                (wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon || 
                wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon)) {

                // Clone the feature 
                OGRFeature *poDstFeature = poFeature->Clone();
                
                if( poLayer_s->CreateFeature( poDstFeature ) != OGRERR_NONE ){
                    printf( "Failed to create feature in shapefile.\n" );
                    exit( 1 );
                }
                
                OGRFeature::DestroyFeature( poDstFeature );
            }

            OGRFeature::DestroyFeature( poFeature );
        }        

        GDALClose(poDS_s);
    }

    GDALClose(poDS_g);

    printf("GPKG Tiling finished. Runtime: %.2f sec\n", 1.0e-6*(walltime_usec() - start_usec) );
}
