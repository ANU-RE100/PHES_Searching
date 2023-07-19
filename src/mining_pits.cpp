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
	vector<ArrayCoordinate> polygon_boundary = find_edge(polygon_points, true);

	// Find the point with the largest distance between all boundary points
	double max_clearance = 0;
	ArrayCoordinate pole_point = polygon_points[0];
	for (ArrayCoordinate fill_point : polygon_points){
		
		double clearance = INT_MAX;
		for (ArrayCoordinate boundary_point : polygon_boundary) {
			double distance = find_distance(fill_point,boundary_point)*1000; // km to m
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
	// Pit lake volumes are calculated by assuming that they are conical structures
	// Find height of cone with pit_lake surface as base
	double pit_lake_surface_r = sqrt(10000*pit_lake_max_area / pi); // Area from Ha to m^2

	// Find volume of cone for test fill depth
	double pit_lake_test_r = pit_lake_test_depth * (pit_lake_surface_r / pit_lake_max_depth);
	pit_lake_test_area = pi * pit_lake_test_r * pit_lake_test_r / 10000; // Area from m^2 to Ha
	double test_cone_volume = pi * (pit_lake_test_r * pit_lake_test_r) * pit_lake_test_depth / 3;

	// Estimate volume of pit_lake
	pit_lake_test_volume = (test_cone_volume) / pow(10,6); // Volume from m^3 to GL

	return;
}

void model_pit_lakes(BulkPit &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, 
						Model<bool> *seen_pl, Model<bool> *seen_d, vector<ArrayCoordinate> &individual_pit_lake_points, Model<short> *DEM){
	int seed_row = pit.seed_point.row;
	int seed_col = pit.seed_point.col;	

	double pit_lake_max_area = pit_area_calculator(seed_row, seed_col, pit_lake_mask, depression_mask, seen_pl, seen_d, individual_pit_lake_points, pit.overlap, pit.seed_point);
	pit.pit_lake_area = pit_lake_max_area;

	// Estimate the lowest point of the pit lake
	Circle pole_of_inaccessibility = find_pole_of_inaccessibility(individual_pit_lake_points);
	ArrayCoordinate pit_lake_lowest_point = pole_of_inaccessibility.centre_point;
	double pole_clearance = pole_of_inaccessibility.radius;

	// Estimate maximum depth of the pit lake
	double pit_lake_min_elevation = DEM->get(seed_row,seed_col) - pit_lake_relative_depth * 2*pole_clearance; // Relative depth assumpion x diameter of circle formed by POI clearance. DEM has equal elevation across entire pit lake.
	if (pit_lake_min_elevation < pit.min_elevation){
		pit.lowest_point = pit_lake_lowest_point;
		pit.min_elevation = pit_lake_min_elevation;
	}
	
	// Estimate the pit volumes and areas at different depths
	double pit_lake_max_depth = DEM->get(seed_row,seed_col) - pit_lake_min_elevation;
	uint depth_tests = pit.fill_elevations.size();
	
	if (pit.overlap){
		depth_tests = (double)depth_tests/2+0.5; // Round up if odd dam_wall_heights.size()
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
	}	

	// Add maximum pit lake volume and depth to all depression volumes
	if (pit.overlap)
		for(uint i=(depth_tests); i < pit.fill_elevations.size(); i++){
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
	
	std::vector<GeographicCoordinate> depression_polygon = convert_poly(find_edge(individual_depression_points, true));
	
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
	uint depth_tests = 0;
	
	if (pit.overlap){
		depth_tests = (double)pit.fill_elevations.size()/2+0.5; // Round up if odd dam_wall_heights.size()
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
	for(uint i=(depth_tests); i < pit.fill_elevations.size(); i++) {
		pit.fill_elevations[i] = depression_min_elevation + ((i-depth_tests) / (double)(pit.fill_elevations.size() - depth_tests)) * MAX(lowest_edge_elevation - depression_min_elevation,0);
	}
	
	// Calculate volume and area of depression at different fill elevations
	for (ArrayCoordinate point : individual_depression_points) {
		for(uint i=(depth_tests); i < pit.fill_elevations.size(); i++) {
			if (DEM->get(point.row,point.col) <= pit.fill_elevations[i]) {
				pit.areas[i] += find_area(point); // Ha
				pit.volumes[i] += 0.01*find_area(point) * (pit.fill_elevations[i] - DEM->get(point.row,point.col)); // Convert to GL
			}
		}		
	}

	// Determine each of the test depths and fix areas for pit lakes that are almost entirely filled with water
	for(uint i=(depth_tests); i < pit.fill_elevations.size(); i++) {
		if(pit.overlap)
			pit.areas[i] = MAX(pit.areas[i], pit.areas[depth_tests-1]);
		pit.fill_depths[i] += pit.fill_elevations[i] - depression_min_elevation;
	}

	return;
}