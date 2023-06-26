#include "phes_base.h"
#include "coordinates.h"
#include "model2D.h"
#include "polygons.h"
#include "search_config.hpp"
#include "polygons.h"
#include "mining_pits.h"
#include "constructor_helpers.hpp"

double pit_area_calculator(int row, int col, Model<bool> *pit_mask, Model<bool> *overlap_mask, Model<bool> *seen, Model<bool> *overlap_seen, vector<ArrayCoordinate> &individual_pit_points, bool &pit_overlap, ArrayCoordinate &pit_new_seed){
	double pit_area = 0;

	// Find all cells interconnected within the pit and add them to the individual_pit_points vector
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
			individual_pit_points.push_back(p);
			pit_area += find_area(p);

			if(overlap_mask->get(p.row,p.col) && (pit_overlap == false) && !overlap_seen->get(p.row,p.col)){
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

Circle find_pole_of_inaccessibility(vector<ArrayCoordinate> polygon_points) {
	// The lowest point of a pit lake is assumed to be at the Pole of Inaccessibility (POI) for the pit polygon
	// For a perfect circle, the POI would be the centre of the circle

	// Create a vector of all points on the polygon boundary
	vector<ArrayCoordinate> polygon_boundary = find_edge(polygon_points);

	// Find the point with the largest distance between all boundary points
	double max_clearance = 0;
	ArrayCoordinate pole_point = polygon_points[0];
	for (ArrayCoordinate fill_point : polygon_points){
		
		double clearance = INT_MAX;
		for (ArrayCoordinate boundary_point : polygon_boundary) {
			double distance = find_distance(fill_point,boundary_point);
			clearance = MIN(clearance,distance);
		}
		if(clearance > max_clearance){
			max_clearance = clearance;
			pole_point = fill_point;
		}
	}

	Circle pole = {pole_point,max_clearance};

	return pole;
}

void find_pit_lake_depth_characteristics(int pit_lake_max_depth, double pit_lake_max_area, int pit_lake_test_depth, double &pit_lake_test_volume, double &pit_lake_test_area) {
	// Pit lake volumes are calculated by assuming that they are conical structures with a flat bottom
	// Find height of cone with pit_lake surface as base
	double pit_lake_surface_r = sqrt(10000*pit_lake_max_area / pi); // Area from Ha to m^2
	double pit_lake_bottom_r = sqrt(pit_lake_relative_area * 10000 * pit_lake_max_area / pi);
	double surface_bottom_r_diff = pit_lake_surface_r - pit_lake_bottom_r;
	double large_cone_height = pit_lake_surface_r*(pit_lake_max_depth/surface_bottom_r_diff);

	// Find volume of cone with pit_lake bottom as base
	double small_cone_height = large_cone_height - pit_lake_max_depth;
	double small_cone_volume = pi * (pit_lake_bottom_r * pit_lake_bottom_r) * small_cone_height / 3;

	// Find volume of cone for test depth
	double test_cone_height = small_cone_height + pit_lake_test_depth;
	double pit_lake_test_r = test_cone_height * (pit_lake_bottom_r / small_cone_height);
	pit_lake_test_area = pi * pit_lake_test_r * pit_lake_test_r / 10000; // Area from m^2 to Ha
	double test_cone_volume = pi * (pit_lake_test_r * pit_lake_test_r) * test_cone_height / 3;

	// Estimate volume of pit_lake
	pit_lake_test_volume = (test_cone_volume - small_cone_volume) / pow(10,6); // Volume from m^3 to GL

	return;
}

double determine_circularity(std::vector<ArrayCoordinate> individual_pit_points, ArrayCoordinate lowest_point, double pit_area){
	double pit_radius = sqrt(pit_area / pi);
	double area_in_circle = 0;
	double pit_circularity = 0;

	for (ArrayCoordinate point : individual_pit_points) {
		double distance_to_poi = find_distance(lowest_point, point);

		if (distance_to_poi <= pit_radius)
			area_in_circle+=find_area(point);
	}

	pit_circularity = area_in_circle / pit_area;

	return pit_circularity;
}

void model_pit_lakes(BulkPit &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, 
						Model<bool> *seen_pl, Model<bool> *seen_d, vector<ArrayCoordinate> &individual_pit_lake_points, Model<short> *DEM){
	int seed_row = pit.seed_point.row;
	int seed_col = pit.seed_point.col;	

	double pit_lake_max_area = pit_area_calculator(seed_row, seed_col, pit_lake_mask, depression_mask, seen_pl, seen_d, individual_pit_lake_points, pit.overlap, pit.seed_point);
	pit.pit_lake_area = pit_lake_max_area;

	// Estimate the lowest point of the pit lake
	ArrayCoordinate pit_lake_lowest_point = find_pole_of_inaccessibility(individual_pit_lake_points).centre_point;

	// Estimate maximum depth of the pit lake
	double pit_lake_min_elevation = DEM->get(seed_row,seed_col) - pit_lake_relative_depth * 2*sqrt(pit_lake_max_area * 10000 / pi); // Relative depth assumpion x diameter of pit, assuming it is a circle. DEM has equal elevation across entire pit lake.
	if (pit_lake_min_elevation < pit.min_elevation){
		pit.lowest_point = pit_lake_lowest_point;
		pit.min_elevation = pit_lake_min_elevation;
	}
	
	// Estimate the pit volumes and areas at different depths
	double pit_lake_max_depth = DEM->get(seed_row,seed_col) - pit_lake_min_elevation;
	uint depth_tests = pit.fill_elevations.size();
	
	if (pit.overlap == true){
		depth_tests = depth_tests/2;
	}

	for(uint i=0; i < depth_tests; i++) {
		int pit_lake_test_depth = ((i+1.0)/depth_tests) * pit_lake_max_depth;
		double pit_lake_test_volume = 0;
		double pit_lake_test_area = 0;

		find_pit_lake_depth_characteristics(pit_lake_max_depth, pit_lake_max_area, pit_lake_test_depth, pit_lake_test_volume, pit_lake_test_area);

		pit.fill_depths[i] = pit_lake_test_depth;
		pit.areas[i] = pit_lake_test_area;
		pit.volumes[i] += pit_lake_test_volume;
		pit.fill_elevations[i] = DEM->get(seed_row, seed_col) - pit_lake_test_depth;
		//printf("%i %i %i %.2f %.2f %.2f %i\n",i,pit.fill_elevations[i],pit_lake_test_depth, pit_lake_test_area, pit_lake_test_volume, pit_lake_max_depth, pit.overlap);
	}	

	// Add maximum pit lake volume and depth to all depression volumes
	if (pit.overlap == true)
		for(uint i=(pit.fill_elevations.size()-depth_tests); i < pit.fill_elevations.size(); i++){
			pit.volumes[i] += pit.volumes[depth_tests-1];
			pit.fill_depths[i] += pit_lake_max_depth;
		}

	return;
}

void model_depression(BulkPit &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, 
						Model<bool> *seen_d, Model<bool> *seen_pl, std::vector<ArrayCoordinate> &individual_depression_points, Model<short> *DEM) {

	int seed_row = pit.seed_point.row;
	int seed_col = pit.seed_point.col;	
				
	pit_area_calculator(seed_row, seed_col, depression_mask, pit_lake_mask, seen_d, seen_pl, individual_depression_points, pit.overlap, pit.seed_point);
	
	std::vector<GeographicCoordinate> depression_polygon = convert_coordinates(find_edge(individual_depression_points), 0);

	// Find lowest point on pit edge
	int lowest_edge_elevation = INT_MAX;
	ArrayCoordinate edge_lowest_point = {-1, -1, DEM->get_origin()};

	for (GeographicCoordinate point : depression_polygon) {
		for (uint d=0; d<directions.size(); d++) {
			if ((directions[d].row * directions[d].col != 0))
				continue;
			ArrayCoordinate point_ac = convert_coordinates(point, depression_mask->get_origin());
			ArrayCoordinate neighbor = {point_ac.row+directions[d].row, point_ac.col+directions[d].col, point_ac.origin};
			if (!depression_mask->check_within(neighbor.row,neighbor.col))
				continue;
			if (depression_mask->get(neighbor.row,neighbor.col))
				continue;
			if (DEM->get(neighbor.row,neighbor.col) < lowest_edge_elevation){
				lowest_edge_elevation = DEM->get(neighbor.row,neighbor.col);
				edge_lowest_point = neighbor;
			}
		}
	}

	// Define tests for volume-altitude pairs
	uint depth_tests = pit.fill_elevations.size();
	
	if (pit.overlap == true){
		depth_tests = depth_tests/2;
	}

	// Find the lowest point in the depression
	int depression_min_elevation = INT_MAX;
	for (ArrayCoordinate point : individual_depression_points) {
		depression_min_elevation = MIN(depression_min_elevation,DEM->get(point.row,point.col));
		if (depression_min_elevation < pit.min_elevation){			
			pit.lowest_point = point;
			pit.min_elevation = DEM->get(pit.lowest_point.row,pit.lowest_point.col);
		}	
	}

	// Determing each of the test elevations within the depression
	for(uint i=(pit.fill_elevations.size()-depth_tests); i < pit.fill_elevations.size(); i++) {
		pit.fill_elevations[i] = pit.min_elevation + (i / depth_tests) * (lowest_edge_elevation - pit.min_elevation);
	}
	
	// Calculate volume and area of depression at different fill elevations
	for (ArrayCoordinate point : individual_depression_points) {
		for(uint i=(pit.fill_elevations.size()-depth_tests); i < pit.fill_elevations.size(); i++) {
			if (DEM->get(point.row,point.col) <= pit.fill_elevations[i]) {
				pit.areas[i] += find_area(point);
				pit.volumes[i] += 0.01*find_area(point) * (pit.fill_elevations[i] - DEM->get(point.row,point.col));
			}
		}		
	}

	// Determine each of the test depths and fix areas for pit lakes that area almost entirely filled with water
	for(uint i=(pit.fill_elevations.size()-depth_tests); i < pit.fill_elevations.size(); i++) {
		pit.areas[i] = MAX(pit.areas[i], pit.areas[depth_tests-1]);
		pit.fill_depths[i] += pit.fill_elevations[i] - depression_min_elevation;
		printf("%i %i %i %.2f %.2f %i %i %i %i %.4f %.4f\n",i,pit.fill_elevations[i],pit.fill_depths[i], pit.areas[i], pit.volumes[i], pit.min_elevation, pit.overlap, pit.fill_depths[i], lowest_edge_elevation,convert_coordinates({seed_row, seed_col, pit.seed_point.origin}).lat,convert_coordinates({seed_row, seed_col, pit.seed_point.origin}).lon);
		
	}

	return;
}