#include "phes_base.h"
#include "coordinates.h"
#include "model2D.h"
#include "search_config.hpp"

int convert_to_int(double f)
{
	if(f>=0)
		return (int) (f+0.5);
	else
		return (int) (f-0.5);
}

double max(vector<double> a)
{
	double amax = -1.0e20;
	for (uint ih=0; ih<a.size(); ih++)
		amax = MAX(amax, a[ih]);

	return amax;
}

double convert_to_dam_volume(int height, double length)
{
	return (((height+freeboard)*(cwidth+dambatter*(height+freeboard)))/1000000)*length;
}

double convert_to_dam_volume(double height, double length)
{
	return (((height+freeboard)*(cwidth+dambatter*(height+freeboard)))/1000000)*length;
}

double linear_interpolate(double value, vector<double> x_values, vector<double> y_values)
{
	uint i = 0;
	while (x_values[i]<value-EPS) {
		if (i==x_values.size()-1)
			return INF;
		else
			i++;
	}

	double xlower = (i) ? x_values[i-1] : 0;
	double ylower = (i) ? y_values[i-1] : 0;
	double r = x_values[i]-xlower;

	return (ylower+(y_values[i]-ylower)*(value-xlower)/r);
}

string str(int i)
{
	char buf[32];
	sprintf(buf, "%d", i);
	string to_return(buf);
	return to_return;
}

unsigned long walltime_usec()
{
	struct timeval now;
	gettimeofday(&now,(struct timezone*)0);
	return (1000000*now.tv_sec + now.tv_usec);
}

double find_required_volume(int energy, int head)
{
	return (((double)(energy)*J_GWh_conversion)/((double)(head)*water_density*gravity*generation_efficiency*usable_volume*cubic_metres_GL_conversion));
}

char* convert_string(string str){
  // NUKE THIS, mem leak galore
	char *c = new char[str.length() + 1];
	strcpy(c, str.c_str());
	return c;
}

string dtos(double f, int nd) {
	stringstream ss;
	ss << fixed << std::setprecision(nd) << f;
	return ss.str();
}

std::string get_dem_filename(GridSquare gs){
	std::string to_return;
	if (dem_type == "SRTM") {
		to_return = file_storage_location+"/input/DEMs/"+str(gs)+"_1arc_v3.tif";
	}
	else if (dem_type == "FABDEM"){
		to_return = file_storage_location+"/input/FABDEMs/"+str_fabdem(gs)+"_FABDEM_V1-2.tif";
	}
	else {
		printf("Invalid dem_type specified.");
		exit(1);
	}
	return to_return;
}

Model<short>* read_DEM_with_borders(GridSquare sc, int border){
	Model<short>* DEM = new Model<short>(0, 0, MODEL_UNSET);
	const int neighbors[9][4][2] = {
		//[(Tile coordinates) , 				(Tile base)		 		  									, (Tile limit)				  																				, (Tile offset)	 	       ]
		{ {sc.lat  ,sc.lon  } , {border,      						border	 						}, {border+model_size-tile_overlap, 	model_size+border-tile_overlap		}, {border-tile_overlap,    			border  		  				 	} },
		{ {sc.lat+1,sc.lon-1} , {0,			  						0		 						}, {border, 	    					border	 	 						}, {border-model_size, 	   				border-(model_size-tile_overlap)	} },
		{ {sc.lat+1,sc.lon  } , {0,	      	  						border	 						}, {border,	    						model_size+border-tile_overlap		}, {border-model_size, 					border     							} },
		{ {sc.lat+1,sc.lon+1} , {0,	      	  						model_size+border-tile_overlap	}, {border,        						model_size+2*border-tile_overlap	}, {border-model_size, 					border+(model_size-tile_overlap)	} },
		{ {sc.lat  ,sc.lon+1} , {border-tile_overlap,    			model_size+border-tile_overlap	}, {model_size+border-tile_overlap,   	model_size+2*border-tile_overlap	}, {border-tile_overlap,   				border+(model_size-tile_overlap)	} },
		{ {sc.lat-1,sc.lon+1} , {model_size+border-tile_overlap,	model_size+border-tile_overlap	}, {model_size+2*border-tile_overlap, 	model_size+2*border-tile_overlap	}, {border+(model_size-2*tile_overlap),	border+(model_size-tile_overlap)	} },
		{ {sc.lat-1,sc.lon  } , {model_size+border-tile_overlap,	border							}, {model_size+2*border-tile_overlap, 	model_size+border					}, {border+(model_size-2*tile_overlap), border     							} },
		{ {sc.lat-1,sc.lon-1} , {model_size+border-tile_overlap,	0		 	 					}, {model_size+2*border-tile_overlap, 	border	 							}, {border+(model_size-2*tile_overlap), border-(model_size-tile_overlap)	} },
		{ {sc.lat  ,sc.lon-1} , {border-tile_overlap,    			0		 	 					}, {model_size+border-tile_overlap,   	border	 	 						}, {border-tile_overlap,    			border-(model_size-tile_overlap)	} }
	};
	for (int i=0; i<9; i++) {
		GridSquare gs = GridSquare_init(neighbors[i][0][0], neighbors[i][0][1]);
		ArrayCoordinate tile_start = ArrayCoordinate_init(neighbors[i][1][0], neighbors[i][1][1], get_origin(gs, border));
		ArrayCoordinate tile_end = ArrayCoordinate_init(neighbors[i][2][0], neighbors[i][2][1], get_origin(gs, border));
		ArrayCoordinate tile_offset = ArrayCoordinate_init(neighbors[i][3][0], neighbors[i][3][1], get_origin(gs, border));
		try{
			Model<short>* DEM_temp = new Model<short>(get_dem_filename(gs), GDT_Int16);
			if (i==0) {
				DEM = new Model<short>(DEM_temp->nrows()+2*border-tile_overlap,DEM_temp->ncols()+2*border-tile_overlap, MODEL_SET_ZERO);
				DEM->set_geodata(DEM_temp->get_geodata());
				GeographicCoordinate origin = get_origin(gs, border);
				DEM->set_origin(origin.lat, origin.lon);
			}
			for(int row = tile_start.row ; row < tile_end.row ; row++)
				for(int col = tile_start.col ; col < tile_end.col; col++){
					DEM->set(row, col, DEM_temp->get(row-tile_offset.row,col-tile_offset.col));
				}
			delete DEM_temp;
		}catch (int e){
			search_config.logger.debug("Could not find file "+get_dem_filename(gs)+" " + strerror(errno));
			if (i==0)
				throw(1);
		}
	}
	return DEM;
}


BigModel BigModel_init(GridSquare sc){
	BigModel big_model;
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
		big_model.neighbors[i] = neighbors[i];
	}
	big_model.DEM = read_DEM_with_borders(sc, (model_size-tile_overlap));
	for(int i = 0; i<9; i++){
		GridSquare gs = big_model.neighbors[i];
		try{
			big_model.flow_directions[i] = new Model<char>(file_storage_location+"processing_files/flow_directions/"+str(gs)+"_flow_directions.tif",GDT_Byte);
		}catch(int e){
			search_config.logger.debug("Could not find " + str(gs));
		}
	}
	return big_model;
}

double calculate_power_house_cost(double power, double head){
	return powerhouse_coeff*pow(MIN(power,800),(power_exp))/pow(head,head_exp);
}

double calculate_tunnel_cost(double power, double head, double seperation){
	return ((power_slope_factor*MIN(power,800)+slope_int)*pow(head,head_coeff)*seperation*1000)+(power_offset*MIN(power,800)+tunnel_fixed);
}

void set_FOM(Pair* pair){
	double seperation = pair->distance;
	double head = (double)pair->head;
	double power = 1000*pair->energy_capacity/pair->storage_time;
	double energy_cost = dam_cost*1/(pair->water_rock*generation_efficiency * usable_volume*water_density*gravity*head)*J_GWh_conversion/cubic_metres_GL_conversion;
	double power_cost;
	double tunnel_cost;
	double power_house_cost;
	if (head > 800) {
		power_house_cost = 2*calculate_power_house_cost(power/2, head/2);
		tunnel_cost = 2*calculate_tunnel_cost(power/2, head/2, seperation);
		power_cost = 0.001*(power_house_cost+tunnel_cost)/MIN(power, 800);
	}
	else {
		power_house_cost = calculate_power_house_cost(power, head);
		tunnel_cost = calculate_tunnel_cost(power, head, seperation);
		power_cost = 0.001*(power_house_cost+tunnel_cost)/MIN(power, 800);
		if(pair->lower.ocean){
			double total_lining_cost = lining_cost*pair->upper.area*meters_per_hectare;
			power_house_cost = power_house_cost*sea_power_scaling;
			double marine_outlet_cost = ref_marine_cost*power*ref_head/(ref_power*head);
			power_cost = 0.001*((power_house_cost+tunnel_cost)/MIN(power, 800) + marine_outlet_cost/power);
			energy_cost += 0.000001*total_lining_cost/pair->energy_capacity;
		}
	}

	pair->FOM = power_cost+energy_cost*pair->storage_time;
	pair->category = 'Z';
	uint i = 0;
	while(i<category_cutoffs.size() && pair->FOM<category_cutoffs[i].power_cost+pair->storage_time*category_cutoffs[i].storage_cost){
		pair->category = category_cutoffs[i].category;
		i++;
	}
}

string energy_capacity_to_string(double energy_capacity){
	return to_string(convert_to_int(energy_capacity));
}

string str(Test test){
	return energy_capacity_to_string(test.energy_capacity)+"GWh_"+to_string(test.storage_time)+"h";
}

bool file_exists (char* name) {
	ifstream infile(name);
    return infile.good();
}

bool file_exists (string name) {
	ifstream infile(name.c_str());
    return infile.good();
}

ExistingReservoir get_existing_reservoir(string name) {
  ExistingReservoir to_return;
  string filename = file_storage_location + "input/existing_reservoirs/" + existing_reservoirs_csv;
  if(!file_exists(convert_string(filename))){
    cout << "File " << filename << " does not exist." << endl;
    throw 1;
  }
  vector<ExistingReservoir> reservoirs = read_existing_reservoir_data(
      convert_string(filename));

  for (ExistingReservoir r : reservoirs)
    if (r.identifier == name)
      to_return = r;

  int i = 0;
  for (string s : read_names(convert_string(file_storage_location +
                                            "input/existing_reservoirs/" +
                                            existing_reservoirs_shp_names))) {
    if (s == name)
      break;
    else
      i++;
  }

  filename = file_storage_location + "input/existing_reservoirs/" +
                    existing_reservoirs_shp;
  char *shp_filename = new char[filename.length() + 1];
  strcpy(shp_filename, filename.c_str());
  if (!file_exists(shp_filename)) {
    search_config.logger.debug("No file: " + filename);
    throw(1);
  }
  SHPHandle SHP = SHPOpen(convert_string(filename), "rb");
  if (SHP != NULL) {
    int nEntities;
    SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL);

    SHPObject *shape;
    shape = SHPReadObject(SHP, i);
    if (shape == NULL) {
      fprintf(stderr, "Unable to read shape %d, terminating object reading.\n",
              i);
      throw(1);
    }
    for (int j = 0; j < shape->nVertices; j++) {
      // if(shape->panPartStart[iPart] == j )
      //  break;
      GeographicCoordinate temp_point =
          GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
      to_return.polygon.push_back(temp_point);
    }
    SHPDestroyObject(shape);
  } else {
    cout << "Could not read shapefile " << filename << endl;
    throw(1);
  }
  SHPClose(SHP);

  to_return.area = geographic_polygon_area(to_return.polygon);

  return to_return;
}

vector<ExistingReservoir> get_existing_reservoirs(GridSquare grid_square) {
  vector<ExistingReservoir> to_return;
  string filename = file_storage_location + "input/existing_reservoirs/" + existing_reservoirs_csv;
  if (!file_exists(convert_string(filename))) {
    cout << "File " << filename << " does not exist." << endl;
    throw 1;
  }
  vector<ExistingReservoir> reservoirs = read_existing_reservoir_data(convert_string(filename));

  vector<string> names = read_names(convert_string(
      file_storage_location + "input/existing_reservoirs/" + existing_reservoirs_shp_names));

  filename = file_storage_location + "input/existing_reservoirs/" + existing_reservoirs_shp;

  char *shp_filename = new char[filename.length() + 1];
  strcpy(shp_filename, filename.c_str());
  if (!file_exists(shp_filename)) {
    search_config.logger.debug("No file: " + filename);
    throw(1);
  }
  SHPHandle SHP = SHPOpen(convert_string(filename), "rb");
  if (SHP != NULL) {
    int nEntities;
    vector<vector<GeographicCoordinate>> relevant_polygons;
    SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL);

    SHPObject *shape;
    for(int i = 0; i<nEntities; i++){
      shape = SHPReadObject(SHP, i);
      if (shape == NULL) {
        fprintf(stderr, "Unable to read shape %d, terminating object reading.\n", i);
      }
      int idx = -1;
      for (uint r = 0; r < reservoirs.size(); r++) {
        if (reservoirs[r].identifier == names[i]) {
          idx = r;
        }
      }
      if (idx < 0) {
        search_config.logger.debug("Could not find reservoir with id " + names[i]);
      }
      ExistingReservoir reservoir = reservoirs[idx];
      //GeographicCoordinate gc = GeographicCoordinate_init(reservoir.latitude, reservoir.longitude);
      //if(!check_within(gc, grid_square))
        //continue;
      vector<GeographicCoordinate> temp_poly;
      for (int j = 0; j < shape->nVertices; j++) {
        // if(shape->panPartStart[iPart] == j )
        //  break;
        GeographicCoordinate temp_point =
            GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
        reservoir.polygon.push_back(temp_point);
      }
	  // Coordinates in existing_reservoirs_csv are based on geometric centre.
	  // Require the same calculation of coordinates here to prevent
	  // disconnect between reservoirs.csv and existing_reservoirs_csv
      bool overlaps_grid_cell = false;
	  double centre_gc_lat = 0;
	  double centre_gc_lon = 0;

      for(GeographicCoordinate gc : reservoir.polygon) {
		centre_gc_lat += gc.lat;
		centre_gc_lon += gc.lon;
	  }
	  GeographicCoordinate centre_gc = GeographicCoordinate_init(centre_gc_lat / reservoir.polygon.size(), centre_gc_lon / reservoir.polygon.size());
	  if(check_within(centre_gc, grid_square)){
        overlaps_grid_cell = true;
      }

      SHPDestroyObject(shape);
      if(overlaps_grid_cell){
        reservoir.area = geographic_polygon_area(reservoir.polygon);
        to_return.push_back(reservoir);
      }
    }
  } else {
    cout << "Could not read shapefile " << filename << endl;
    throw(1);
  }
  SHPClose(SHP);
  return to_return;
}

RoughBfieldReservoir existing_reservoir_to_rough_reservoir(ExistingReservoir r){
	RoughBfieldReservoir reservoir;
	reservoir.identifier = r.identifier;
    reservoir.brownfield = true;
    reservoir.ocean = false;
	reservoir.latitude = r.latitude;
	reservoir.longitude = r.longitude;
	reservoir.elevation = r.elevation;
	reservoir.bottom_elevation = r.elevation;
	for(uint i = 0; i<dam_wall_heights.size(); i++){
		reservoir.volumes.push_back(r.volume);
		reservoir.dam_volumes.push_back(0);
		reservoir.areas.push_back(r.area);
		reservoir.water_rocks.push_back(1000000000);
  }

	GeographicCoordinate origin = get_origin(search_config.grid_square, border);
	for(GeographicCoordinate c : r.polygon)
    reservoir.shape_bound.push_back(convert_coordinates(c, origin));
	return reservoir;
}

vector<ExistingPit> get_pit_details(GridSquare grid_square){
	vector<ExistingPit> gridsquare_pits;

	vector<ExistingPit> pits = read_existing_pit_data(convert_string(file_storage_location+"input/existing_reservoirs/"+existing_reservoirs_csv));

	for(ExistingPit p : pits){
		if (check_within(GeographicCoordinate_init(p.reservoir.latitude, p.reservoir.longitude), grid_square))
			gridsquare_pits.push_back(p);
	}
	return gridsquare_pits;
}

ExistingPit get_pit_details(string pitname){
	ExistingPit pit;
	vector<ExistingPit> pits = read_existing_pit_data(convert_string(file_storage_location+"input/existing_reservoirs/"+existing_reservoirs_csv));

	for(ExistingPit p : pits){
		if (p.reservoir.identifier==pitname)
			pit = p;
	}
	return pit;
}

RoughBfieldReservoir pit_to_rough_reservoir(BulkPit pit, GeographicCoordinate lowest_point){
	RoughBfieldReservoir reservoir;
	reservoir.identifier = pit.res_identifier;
	reservoir.pit = true;
    reservoir.brownfield = true;
    reservoir.ocean = false;
	reservoir.turkey = false;
	reservoir.latitude = lowest_point.lat;
	reservoir.longitude = lowest_point.lon;
	reservoir.elevation = pit.min_elevation;

	for(uint i = 0; i<pit.fill_elevations.size(); i++){
		reservoir.volumes.push_back(pit.volumes[i]);
		reservoir.fill_depths.push_back(pit.fill_depths[i]);
		reservoir.areas.push_back(pit.areas[i]);
		reservoir.water_rocks.push_back(1000000000);
		reservoir.dam_volumes.push_back(0);
  	}

	GeographicCoordinate origin = get_origin(search_config.grid_square, border);
	for(GeographicCoordinate c : pit.brownfield_polygon)
    	reservoir.shape_bound.push_back(convert_coordinates(c, origin));
	return reservoir;
}