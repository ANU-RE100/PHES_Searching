#include "phes_base.h"
#include "coordinates.h"
#include "model2D.h"
#include "search_config.hpp"

std::string get_mining_tenament_path(){
	std::string lat_prefix;
    std::string lon_prefix;

	if (search_config.grid_square.lat >= 0)
        lat_prefix = "n";
	else
        lat_prefix = "s";
    if (search_config.grid_square.lon >= 0)
        lon_prefix = "e";
    else
        lon_prefix = "w";
        
    // Define the latitude/longitude string identifiers
    std::stringstream ss_lat;
    std::stringstream ss_lon;
    ss_lat << std::setw(2) << std::setfill('0') << abs(search_config.grid_square.lat);
    ss_lon << std::setw(3) << std::setfill('0') << abs(search_config.grid_square.lon);
    std::string lat_leading = ss_lat.str();
    std::string lon_leading = ss_lon.str();

    std::string lat_str = lat_prefix + lat_leading;
    std::string lon_str = lon_prefix + lon_leading;

	string filename = mining_tenament_shp;
	filename += lat_str + "_" + lon_str + ".shp";

	return filename;
}


// Remove this deprecated functionality???
void depression_volume_finding(Model<short>* DEM) {
	vector<vector<string> > csv_modified_lines;
	vector<int> csv_modified_line_numbers;
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
		SHPGetInfo(SHP, &nEntities, NULL, NULL, NULL);

		SHPObject *shape;

		for(int i = 0; i<nEntities; i++){
			Model<bool>* extent = new Model<bool>(DEM->nrows(), DEM->ncols(), MODEL_SET_ZERO);
			extent->set_geodata(DEM->get_geodata());
			short min_elevation = 32767;
			short max_elevation = 0;
			vector<string> csv_modified_line(2*num_altitude_volume_pairs+2);

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
			GeographicCoordinate gc = GeographicCoordinate_init(reservoir.latitude, reservoir.longitude);
			if(!check_within(gc, search_config.grid_square)) {
				SHPDestroyObject(shape);
				delete extent;
				continue;
			}
			vector<GeographicCoordinate> temp_poly;
			for (int j = 0; j < shape->nVertices; j++) {
				GeographicCoordinate temp_point = GeographicCoordinate_init(shape->padfY[j], shape->padfX[j]);
				temp_poly.push_back(temp_point);
					
			}
			polygon_to_raster(temp_poly, extent);
			SHPDestroyObject(shape);

			// Find lowest elevation within mine polygon (pour point)
			for(int row = 0; row<extent->nrows(); row++)
				for(int col = 0; col<extent->ncols(); col++){
					if(extent->get(row, col)) {
						min_elevation = MIN(DEM->get(row, col), min_elevation);
						max_elevation = MAX(DEM->get(row,col), max_elevation);
					}
				}

			double area_at_elevation[max_elevation + 1] = {0};
			double volume_at_elevation[max_elevation + 1] = {0};
			double cumulative_area_at_elevation[max_elevation + 1] = {0};
			double pit_elevations[num_altitude_volume_pairs] = {0};

			// Determine the elevations for altitude-volume pairs
			for (int ih = 1; ih <= num_altitude_volume_pairs; ih++) {
				pit_elevations[ih-1] = min_elevation + std::round(ih * (max_elevation - min_elevation)/num_altitude_volume_pairs);
			}

			// Find the area of cells within mine polygon at each elevation above the pour point
			for(int row = 0; row<extent->nrows(); row++)
				for(int col = 0; col<extent->ncols(); col++)
					if(extent->get(row, col)){
						area_at_elevation[min_elevation + 1] += find_area(ArrayCoordinate_init(row, col, DEM->get_origin()));
					}

			// Find the surface area and volume of reservoir at each elevation above pour point 
			for (int ih=1; ih<max_elevation+1-min_elevation;ih++) {
				cumulative_area_at_elevation[min_elevation + ih] = cumulative_area_at_elevation[min_elevation + ih-1] + area_at_elevation[min_elevation + ih];
				volume_at_elevation[min_elevation + ih] = volume_at_elevation[min_elevation + ih-1] + 0.01*cumulative_area_at_elevation[min_elevation + ih]; // area in ha, vol in GL
			}

			// Find the altitude-volume pairs for the pit
			csv_modified_line[0] = to_string(min_elevation);
			// csv_modified_line[1] = to_string(volume_at_elevation[max_elevation]); // MAIN CODE USE THIS
			csv_modified_line[1] = to_string(3 + volume_at_elevation[max_elevation]); // ASEAN CODE USE THIS
			for (int ih =0 ; ih < num_altitude_volume_pairs; ih++) {
				int height = pit_elevations[ih];
				csv_modified_line[2+2*ih] = to_string(height);
				// csv_modified_line[2+2*ih + 1] = to_string(volume_at_elevation[height]); MAIN CODE USE THIS
				csv_modified_line[2+2*ih + 1] = to_string(3 + volume_at_elevation[height]); // ASEAN CODE USE THIS
			}

			// Add the line to the vector to be written to the pits CSV
			csv_modified_lines.push_back(csv_modified_line);
			csv_modified_line_numbers.push_back(i+1);

			//extent->write(file_storage_location+"debug/extent/"+str(search_config.grid_square)+"_extent.tif", GDT_Byte);

			delete extent;
			temp_poly.clear();
		}
	} else {
		cout << "Could not read shapefile " << filename << endl;
		throw(1);
	}
	SHPClose(SHP);			

	// Write the altitude-volume pairs to the CSV
	std::ifstream inputFile(file_storage_location + "input/existing_reservoirs/" + existing_reservoirs_csv);
	vector<string> lines;

	if (!inputFile.is_open()) {
		printf("Error opening pits CSV\n");
	}

	string line;
	int line_number = 0;
	while(std::getline(inputFile, line)) {
		istringstream lineStream(line);
		string cell;
		int column = 1;
		ostringstream modifiedLine;

		while (std::getline(lineStream, cell, ',')) {
			
			if (column >= 4 && column <= 5 + 2*num_altitude_volume_pairs && std::count(csv_modified_line_numbers.begin(), csv_modified_line_numbers.end(), line_number)) {
				std::vector<int>::iterator vector_index_itr = find(csv_modified_line_numbers.begin(), csv_modified_line_numbers.end(), line_number);
				int vector_index = std::distance(csv_modified_line_numbers.begin(), vector_index_itr);
				cell = string(csv_modified_lines[vector_index][column-4]); 
			}

			modifiedLine << cell;

			if (column < 5 + 2*num_altitude_volume_pairs) {
				modifiedLine << ",";
			}

			++column;
		}

		lines.push_back(modifiedLine.str());

		line_number++;
	}

	inputFile.close();

	std::ofstream outputFile(file_storage_location + "input/existing_reservoirs/" + existing_reservoirs_csv);
	for (const auto &line : lines) {
		outputFile << line << std::endl;
	}

	outputFile.close();
}

double pit_area_calculator(int row, int col, Model<bool> *pit_mask, Model<bool> *seen, Model<bool> *individual_pit_mask, std::vector<GeographicCoordinate> &brownfield_polygon){
	double pit_lake_area = 0;

	// Find all cells interconnected to within the pit and add them to the individual mask
	ArrayCoordinate c = ArrayCoordinate_init(row,col,pit_mask->get_origin());
	queue<ArrayCoordinate> q;
	q.push(c);

	while (!q.empty()) {
		ArrayCoordinate p = q.front();
		bool edge_check = false;
		q.pop();

		if (pit_mask->get(p.row,p.col)){
			seen->set(p.row,p.col,true);
			
			individual_pit_mask->set(p.row,p.col,true);
			pit_lake_area += 10000*find_area(p);

			// Add all perpendicular neighbors to the queue. If at least one prependicular neighbor is outside the mask, p is on the pit edge
			for (uint d=0; d<directions.size(); d++) {
				if (!pit_mask->get(row,col))
					edge_check = true;

				ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
				if (!pit_mask->check_within(neighbor.row,neighbor.col))
					continue;
				if ((directions[d].row * directions[d].col == 0) && !seen->get(neighbor.row,neighbor.col)) {
					q.push(neighbor);
				}
			}

			if (edge_check){
				GeographicCoordinate edge_point = pit_mask->get_coordinate(p.row,p.col);
				brownfield_polygon.push_back(edge_point);
			}
		}
	}

	return pit_lake_area;
}

ArrayCoordinate find_lowest_point_pit_lake(Model<bool> *individual_pit_mask) {
	// The lowest point is assumed to be at the Point of Inaccessibility (POI) for the pit polygon
	// For a perfect circle, the POI would be the centre of the circle
	
	queue<ArrayCoordinate> q;
    vector<vector<int>> dist(individual_pit_mask->nrows(), vector<int>(individual_pit_mask->ncols(), INT_MAX));
    
    // Push all boundary cells (value 0 in raster) into the queue and set their distance to 0
    for(int row = 0; row<individual_pit_mask->nrows();row++) {
		for(int col = 0; col<individual_pit_mask->ncols();col++) {	
            if (!individual_pit_mask->get(row,col)) {
				ArrayCoordinate bc = ArrayCoordinate_init(row,col,individual_pit_mask->get_origin());
                q.push(bc);
            }
        }
    }

	// Find distance between internal polygon cells and the polygon boundary
    while (!q.empty()) {
        ArrayCoordinate p = q.front();
        q.pop();

		for (uint d=0; d<directions.size(); d++) {
			ArrayCoordinate neighbor = ArrayCoordinate_init(p.row + directions[d].row, p.col + directions[d].col, individual_pit_mask->get_origin());
			if (directions[d].row * directions[d].col == 0)
				continue;
			if ((!individual_pit_mask->check_within(neighbor.row,neighbor.col)) || (!individual_pit_mask->get(neighbor.row,neighbor.col))) {
				continue;
			}

			if (dist[p.row][p.col] + 1 < dist[neighbor.row][neighbor.col]) {
				dist[neighbor.row][neighbor.col] = dist[p.row][p.col] + 1;
				q.push(neighbor);
			}
		}
    }
	
	// Determine the Point of Inaccessibility based on the maximum distance between any cell and the polygon boundary
	int max_distance = -1;
    ArrayCoordinate lowest_point = {-1, -1, individual_pit_mask->get_origin()};
    
    for(int row = 0; row<individual_pit_mask->nrows();row++) {
		for(int col = 0; col<individual_pit_mask->ncols();col++) {
            if (dist[row][col] > max_distance) {
                max_distance = dist[row][col];
                lowest_point = {row, col, individual_pit_mask->get_origin()};
            }
        }
    }

	return lowest_point;
}

double find_volume_pit_lake(double pit_area, int pit_depth) {
	// Pit lake volumes are calculated by assuming that they are conical structures with a flat bottom
	// Find height of cone with pit_lake surface as base
	double pit_lake_surface_r = sqrt(pit_area / pi);
	double pit_lake_bottom_r = sqrt(pit_lake_relative_area * pit_area / pi);
	double surface_bottom_r_diff = pit_lake_surface_r - pit_lake_bottom_r;
	double wall_vertical_angle = atan(surface_bottom_r_diff/pit_depth);
	double wall_horizontal_angle = pi/2 - wall_vertical_angle;
	double large_cone_height = pit_lake_surface_r*tan(wall_horizontal_angle);

	// Find volume of cone with pit_lake surface as base
	double large_cone_volume = pi * (pit_lake_surface_r * pit_lake_surface_r) * large_cone_height / 3;

	// Find volume of cone with pit_lake bottom as base
	double small_cone_height = large_cone_height - pit_depth;
	double small_cone_volume = pi * (pit_lake_bottom_r * pit_lake_bottom_r) * small_cone_height / 3;

	// Estimate volume of pit_lake
	double pit_volume = large_cone_volume - small_cone_volume;

	return pit_volume;
}

void find_depression_attributes(Model<bool> *individual_pit_mask, Model<short> *DEM, ArrayCoordinate &lowest_point, double &pit_volume, int &depression_elevation, vector<GeographicCoordinate> brownfield_polygon){
	int lowest_edge_elevation = INT_MAX;

	// Find lowest point on pit edge
	for (GeographicCoordinate point : brownfield_polygon) {
		if (DEM->get(point) < lowest_edge_elevation)
			lowest_edge_elevation = DEM->get(point);
	}
	
	for(int row = 0; row<individual_pit_mask->nrows();row++) {
		for(int col = 0; col<individual_pit_mask->ncols();col++) {
            if (!individual_pit_mask->get(row,col)) {
                continue;
            }
			// Find lowest point in depression
			if (DEM->get(row,col) < depression_elevation){
				depression_elevation = DEM->get(row,col) < depression_elevation;
				lowest_point = ArrayCoordinate_init(row,col,DEM->get_origin());
			}
			// Calculate volume of depression
			if (DEM->get(row,col) < lowest_edge_elevation) {
				ArrayCoordinate c = {row,col,DEM->get_origin()};
				pit_volume += find_area(c) * (lowest_edge_elevation - DEM->get(row,col));
			}
        }
    }
	return;
}

double determine_circularity(Model<bool> *individual_pit_mask, ArrayCoordinate lowest_point, double pit_area){
	double pit_radius = sqrt(pit_area / pi);
	double area_in_circle = 0;
	double pit_circularity = 0;

	for(int row = 0; row<individual_pit_mask->nrows();row++) {
		for(int col = 0; col<individual_pit_mask->ncols();col++) {
			if (!individual_pit_mask->get(row,col))
				continue;
			ArrayCoordinate c = {row,col,individual_pit_mask->get_origin()};
			double distance_to_poi = find_distance(lowest_point,c);

			if (distance_to_poi <= pit_radius)
				area_in_circle+=find_area(c);
		}
	}

	pit_circularity = area_in_circle / pit_area;

	return pit_circularity;
}