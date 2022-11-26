#include "phes_base.h"

SearchConfig search_config;

int main(int nargs, char **argv)
{
  if(nargs < 3){
    cout << "Not enough arguements" << endl;
    return -1;
  }
	GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));
  search_config.logger = Logger::DEBUG;

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
