#ifndef MINING_PITS_H
#define MINING_PITS_H

#include <bits/stdc++.h>
#include "coordinates.h"
#include "model2D.h"

using namespace std;

struct BulkPit {
    std::vector<double> areas;
	int min_elevation;
	std::vector<double> volumes;
	std::vector<int> fill_elevations;
    std::vector<int> fill_depths;
	double circularity;
	ArrayCoordinate lowest_point;
	ArrayCoordinate seed_point;
	ArrayCoordinate overlap_point;
	std::string res_identifier;
	std::vector<GeographicCoordinate> brownfield_polygon;		
	bool overlap;	
    double pit_lake_area;

    BulkPit(int row, int col, GeographicCoordinate origin) 
        : areas(10,0), volumes(10,0), fill_elevations(10,0), fill_depths(10,0)
    {
        min_elevation = INT_MAX;
        circularity = 0;
        lowest_point = {-1, -1, origin};
        seed_point = {row,col,origin};
        res_identifier = "unassigned";
        brownfield_polygon = {};		
        overlap = false;
        pit_lake_area = 0;	
    }
};

double pit_area_calculator(int row, int col, Model<bool> *pit_mask, Model<bool> *seen, Model<bool> *individual_pit_mask, std::vector<GeographicCoordinate> &brownfield_polygon);
ArrayCoordinate find_lowest_point_pit_lake(Model<bool> *individual_pit_mask);
double find_volume_pit_lake(double pit_area, int pit_depth);
double determine_circularity(std::vector<ArrayCoordinate> individual_pit_points, ArrayCoordinate lowest_point, double pit_area);
void model_pit_lakes(BulkPit &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, Model<bool> *seen_pl, vector<ArrayCoordinate> &individual_pit_lake_points, Model<short> *DEM);
void model_depression(BulkPit &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, Model<bool> *seen_d, std::vector<ArrayCoordinate> &individual_depression_points, Model<short> *DEM);

#endif