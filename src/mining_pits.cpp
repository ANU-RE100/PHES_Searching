#include "phes_base.h"
#include "coordinates.h"
#include "model2D.h"
#include "polygons.h"
#include "search_config.hpp"
#include "polygons.h"
#include "mining_pits.h"
#include "constructor_helpers.hpp"

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

double pit_area_calculator(int row, int col, Model<bool> *pit_mask, Model<bool> *overlap_mask, Model<bool> *seen, Model<bool> *individual_pit_mask, bool &pit_overlap, ArrayCoordinate &pit_new_seed){
	double pit_area = 0;

	// Find all cells interconnected to within the pit and add them to the individual mask
	ArrayCoordinate c = ArrayCoordinate_init(row,col,pit_mask->get_origin());
	queue<ArrayCoordinate> q;
	q.push(c);

	while (!q.empty()) {
		ArrayCoordinate p = q.front();
		q.pop();

		if(seen->get(p.row,p.col))
			continue;

		seen->set(p.row,p.col,true);

		if(pit_mask->get(p.row,p.col)){			
			individual_pit_mask->set(p.row,p.col,true);
			pit_area += find_area(p);

			if(overlap_mask->get(p.row,p.col)){
				pit_overlap = true;
				pit_new_seed = {p.row,p.col,overlap_mask->get_origin()};
			}

			// Add all perpendicular neighbors to the queue
			for (uint d=0; d<directions.size(); d++) {
				ArrayCoordinate neighbor = {p.row+directions[d].row, p.col+directions[d].col, p.origin};
				if (!pit_mask->check_within(neighbor.row,neighbor.col))
					continue;
				if ((directions[d].row * directions[d].col == 0) && (!seen->get(neighbor.row,neighbor.col))) {
					q.push(neighbor);
				}
			}
		}
	}

	return pit_area;
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
				dist[row][col] = 0;
            }
        }
    }

	// Find distance between internal polygon cells and the polygon boundary
    while (!q.empty()) {
        ArrayCoordinate p = q.front();
        q.pop();

		for (uint d=0; d<directions.size(); d++) {
			ArrayCoordinate neighbor = ArrayCoordinate_init(p.row + directions[d].row, p.col + directions[d].col, individual_pit_mask->get_origin());
			if (directions[d].row * directions[d].col != 0)
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

void model_pit_lakes(PitCharacteristics &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, 
						Model<bool> *seen_pl, Model<bool> *individual_pit_mask, Model<short> *DEM){
	int seed_row = pit.seed_point.row;
	int seed_col = pit.seed_point.col;	

	double pit_lake_area = pit_area_calculator(seed_row, seed_col, pit_lake_mask, depression_mask, seen_pl, individual_pit_mask, pit.pit_overlap, pit.seed_point);
	pit.pit_area = MAX(pit_lake_area,pit.pit_area);
	pit.pit_lake_area = pit_lake_area;

	// Estimate the lowest point of the pit lake
	ArrayCoordinate pit_lake_lowest_point = find_lowest_point_pit_lake(individual_pit_mask);

	// Estimate maximum depth of the pit lake
	double pit_lake_min_elevation = DEM->get(seed_row,seed_col) - pit_lake_relative_depth * 2*sqrt(pit_lake_area / pi); // Relative depth assumpion x diameter of pit, assuming it is a circle. DEM has equal elevation across entire pit lake.
	if (pit_lake_min_elevation < pit.pit_min_elevation){
		pit.lowest_point = pit_lake_lowest_point;
		pit.pit_min_elevation = pit_lake_min_elevation;
	}
	
	// Estimate the pit volume
	double pit_lake_depth = DEM->get(seed_row,seed_col) - pit_lake_min_elevation;
	double pit_lake_volume = find_volume_pit_lake(pit_lake_area, pit_lake_depth);
	pit.pit_depth = MAX(pit.pit_depth,pit_lake_depth);
	pit.pit_volume += pit_lake_volume;

	return;
}

void model_depression(PitCharacteristics &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, 
						Model<bool> *seen_d, Model<bool> *individual_pit_mask, Model<short> *DEM) {

	int seed_row = pit.seed_point.row;
	int seed_col = pit.seed_point.col;	
				
	double depression_area = pit_area_calculator(seed_row, seed_col, depression_mask, pit_lake_mask, seen_d, individual_pit_mask, pit.pit_overlap, pit.seed_point);
	pit.pit_area = MAX(depression_area, pit.pit_area);
	
	ArrayCoordinate offset = ArrayCoordinate_init(0,0,DEM->get_origin());
	ArrayCoordinate edge_point = find_edge(pit.seed_point, individual_pit_mask);

	// Check that there is at least one connected cell
	bool single_point = true;
	for (uint d=0; d<directions.size(); d++) {
		ArrayCoordinate neighbor = ArrayCoordinate_init(edge_point.row + directions[d].row, edge_point.col + directions[d].col, individual_pit_mask->get_origin());
		if (directions[d].row * directions[d].col != 0)
			continue;
		if ((!individual_pit_mask->check_within(neighbor.row,neighbor.col)) || (!individual_pit_mask->get(neighbor.row,neighbor.col))) {
			continue;
		}
		single_point = false;
	}

	std::vector<GeographicCoordinate> depression_polygon;
	if (single_point)
		depression_polygon.push_back(convert_coordinates(edge_point));
	else
	 	depression_polygon = convert_poly(convert_to_polygon(depression_mask, offset, edge_point));

	// Find lowest point on pit edge
	int lowest_edge_elevation = INT_MAX;
	ArrayCoordinate edge_lowest_point = {-1, -1, DEM->get_origin()};

	for (GeographicCoordinate point : depression_polygon) {
		if (DEM->get(point) < lowest_edge_elevation){
			lowest_edge_elevation = DEM->get(point);
			edge_lowest_point = DEM->get_array_coord(point.lat, point.lon);
		}
	}
	
	double depression_elevation_sum = 0;
	double depression_volume = 0;
	int cell_count = 0;

	for(int row = 0; row<individual_pit_mask->nrows();row++) {
		for(int col = 0; col<individual_pit_mask->ncols();col++) {
            if (!individual_pit_mask->get(row,col)) {
                continue;
            }

			// Average elevation is used for pit altitude due to sinks in unfilled DEM
			depression_elevation_sum += DEM->get(row,col);
			cell_count++;

			// Calculate volume of depression
			if (DEM->get(row,col) < lowest_edge_elevation) {
				ArrayCoordinate c = {row,col,DEM->get_origin()};
				depression_volume += 0.01*find_area(c) * (lowest_edge_elevation - DEM->get(row,col));
			}
        }
    }

	if(depression_elevation_sum / cell_count < pit.pit_min_elevation){
		pit.pit_min_elevation = depression_elevation_sum / cell_count;
		pit.lowest_point = edge_lowest_point;
	}
	
	pit.pit_volume += depression_volume;

	return;
}