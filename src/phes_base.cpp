#include "phes_base.h"

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
	char *c = new char[str.length() + 1];
	strcpy(c, str.c_str());
	return c;
}

string dtos(double f, int nd) {
	stringstream ss;
	ss << fixed << std::setprecision(nd) << f;
	return ss.str();
}

Model<short>* read_DEM_with_borders(GridSquare sc, int border){
	Model<short>* DEM = new Model<short>(0, 0, MODEL_UNSET);
	const int neighbors[9][4][2] = {
		//[(Tile coordinates) , (Tile base)		 		  , (Tile limit)				  , (Tile offset)	 	       ]
		{ {sc.lat  ,sc.lon  } , {border,      border	 }, {border+3600,  	3600+border	 }, {border-1,    border     } },
		{ {sc.lat+1,sc.lon-1} , {0,			  0		 	 }, {border, 	    border	 	 }, {border-3601, border-3600} },
		{ {sc.lat+1,sc.lon  } , {0,	      	  border	 }, {border,	    3600+border	 }, {border-3601, border     } },
		{ {sc.lat+1,sc.lon+1} , {0,	      	  3600+border}, {border,        3600+2*border}, {border-3601, border+3600} },
		{ {sc.lat  ,sc.lon+1} , {border-1,    3600+border}, {3600+border,   3600+2*border}, {border-1,    border+3600} },
		{ {sc.lat-1,sc.lon+1} , {3600+border, 3600+border}, {3600+2*border, 3600+2*border}, {border+3599, border+3600} },
		{ {sc.lat-1,sc.lon  } , {3600+border, border	 }, {3600+2*border, 3601+border	 }, {border+3599, border     } },
		{ {sc.lat-1,sc.lon-1} , {3600+border, 0		 	 }, {3600+2*border, border	 	 }, {border+3599, border-3600} },
		{ {sc.lat  ,sc.lon-1} , {border-1,    0		 	 }, {3600+border,   border	 	 }, {border-1,    border-3600} }
	};
	for (int i=0; i<9; i++) {
		GridSquare gs = GridSquare_init(neighbors[i][0][0], neighbors[i][0][1]);
		ArrayCoordinate tile_start = ArrayCoordinate_init(neighbors[i][1][0], neighbors[i][1][1], get_origin(gs, border));
		ArrayCoordinate tile_end = ArrayCoordinate_init(neighbors[i][2][0], neighbors[i][2][1], get_origin(gs, border));
		ArrayCoordinate tile_offset = ArrayCoordinate_init(neighbors[i][3][0], neighbors[i][3][1], get_origin(gs, border));
		try{
			Model<short>* DEM_temp = new Model<short>(file_storage_location+"input/DEMs/"+str(gs)+"_1arc_v3.tif", GDT_Int16);
			if (i==0) {
				DEM = new Model<short>(DEM_temp->nrows()+2*border-1,DEM_temp->ncols()+2*border-1, MODEL_SET_ZERO);
				DEM->set_geodata(DEM_temp->get_geodata());
				GeographicCoordinate origin = get_origin(gs, border);
				DEM->set_origin(origin.lat, origin.lon);
			}
			for(int row = tile_start.row ; row < tile_end.row ; row++)
				for(int col = tile_start.col ; col < tile_end.col; col++)
					DEM->set(row, col, DEM_temp->get(row-tile_offset.row,col-tile_offset.col));
			delete DEM_temp;
		}catch (int e){
			search_config.logger.debug("Could not find file "+file_storage_location+"input/DEMs/"+str(gs)+"_1arc_v3.tif " + strerror(errno));
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
	big_model.DEM = read_DEM_with_borders(sc, 3600);
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

void set_FOM(Pair* pair){
	double seperation = pair->distance;
	double head = (double)pair->head;
	double power = 1000*pair->energy_capacity/pair->storage_time;
	double power_house_cost = powerhouse_coeff*pow(MIN(power,800),(power_exp))/pow(head,head_exp);
	double tunnel_cost = ((power_slope_factor*MIN(power,800)+slope_int)*pow(head,head_coeff)*seperation*1000)+(power_offset*MIN(power,800)+tunnel_fixed);
	double power_cost = 0.001*(power_house_cost+tunnel_cost)/MIN(power, 800);
	double energy_cost = dam_cost*1/(pair->water_rock*generation_efficiency * usable_volume*water_density*gravity*pair->head)*J_GWh_conversion/cubic_metres_GL_conversion;

	if(pair->lower.ocean){
		double total_lining_cost = lining_cost*pair->upper.area*meters_per_hectare;
		power_house_cost = power_house_cost*sea_power_scaling;
		double marine_outlet_cost = ref_marine_cost*power*ref_head/(ref_power*head);
		power_cost = 0.001*((power_house_cost+tunnel_cost)/MIN(power, 800) + marine_outlet_cost/power);
		energy_cost += 0.000001*total_lining_cost/pair->energy_capacity;
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
	if(energy_capacity<10-EPS)
		return dtos(energy_capacity,1);
	else
		return to_string((int)(energy_capacity+EPS));
}

string str(Test test){
	return energy_capacity_to_string(test.energy_capacity)+"GWh_"+to_string(test.storage_time)+"h";
}

bool file_exists (char* name) {
	ifstream infile(name);
    return infile.good();
}

string format_for_filename(string s){
	replace(s.begin(), s.end(), ' ' , '_');
	s.erase(remove(s.begin(), s.end(), '"'), s.end());
	return s;
}

GeographicCoordinate get_origin(double latitude, double longitude, int border){
	return GeographicCoordinate_init(FLOOR(latitude)+1+(border/3600.0),FLOOR(longitude)-(border/3600.0));
}

ExistingReservoir get_existing_reservoir(string name) {
  ExistingReservoir to_return;
  vector<ExistingReservoir> reservoirs = read_existing_reservoir_data(
      convert_string(file_storage_location + "input/existing_reservoirs/" +
                     existing_reservoirs_csv));

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

  string filename = file_storage_location + "input/existing_reservoirs/" +
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
    vector<vector<GeographicCoordinate>> relevant_polygons;
    SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL);

    SHPObject *shape;
    shape = SHPReadObject(SHP, i);
    if (shape == NULL) {
      fprintf(stderr, "Unable to read shape %d, terminating object reading.\n",
              i);
    }
    vector<GeographicCoordinate> temp_poly;
    for (int j = 0; j < shape->nVertices; j++) {
      // if(shape->panPartStart[iPart] == j )
      //  break;
      GeographicCoordinate temp_point =
          GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
      to_return.polygon.push_back(temp_point);
    }
    SHPDestroyObject(shape);
  } else {
    throw(1);
  }
  SHPClose(SHP);

  return to_return;
}

RoughReservoir get_existing_rough_reservoir(string name){
	RoughReservoir reservoir;
	ExistingReservoir r = get_existing_reservoir(name);
	reservoir.identifier = r.identifier;
    reservoir.brownfield = true;
    reservoir.ocean = false;
	reservoir.latitude = r.latitude;
	reservoir.longitude = r.longitude;
	reservoir.elevation = r.elevation;
	for(uint i = 0; i<dam_wall_heights.size(); i++){
		reservoir.volumes.push_back(r.volume);
		reservoir.dam_volumes.push_back(0);
		reservoir.areas.push_back(0);
		reservoir.water_rocks.push_back(1000000000);
	}
	for (uint ih=0; ih < dam_wall_heights.size(); ih++) {
		array<ArrayCoordinate, directions.size()> temp_array;
		for (uint idir=0; idir < directions.size(); idir++) {
			temp_array[idir].row = -100000 * directions[idir].row;
			temp_array[idir].col = -100000 * directions[idir].col;
		}
		reservoir.shape_bound.push_back(temp_array);
	}

	GeographicCoordinate origin = get_origin(r.latitude, r.longitude, border);
	
	for(GeographicCoordinate c : r.polygon){
		ArrayCoordinate p = convert_coordinates(c, origin);
		update_reservoir_boundary(reservoir.shape_bound, p, 0);
	}
	return reservoir;
}

ExistingPit get_pit_details(string pitname){
	ExistingPit pit;
	vector<ExistingPit> pits = read_existing_pit_data(convert_string(file_storage_location+"input/existing_reservoirs/"+existing_reservoirs_csv));

	for(ExistingPit p : pits)
		if(p.reservoir.identifier==pitname)
			pit = p;
	return pit;
}
