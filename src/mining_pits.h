#ifndef MINING_PITS_H
#define MINING_PITS_H

#include <bits/stdc++.h>
#include "coordinates.h"
#include "model2D.h"

using namespace std;

struct PitCharacteristics {
    double pit_area;
	int pit_min_elevation;
	double pit_volume;
	int pit_depth;
	double pit_circularity;
	ArrayCoordinate lowest_point;
	ArrayCoordinate seed_point;
	ArrayCoordinate overlap_point;
	std::string res_identifier;
	std::vector<GeographicCoordinate> brownfield_polygon;		
	bool pit_overlap;	
    double pit_lake_area;

    PitCharacteristics(int row, int col, GeographicCoordinate origin) {
        pit_area = 0;
        pit_min_elevation = INT_MAX;
        pit_volume = 0;
        pit_depth = 0;
        pit_circularity = 0;
        lowest_point = {-1, -1, origin};
        seed_point = {row,col,origin};
        res_identifier = "unassigned";
        brownfield_polygon = {};		
        pit_overlap = false;
        pit_lake_area = 0;	
    }
};

double pit_area_calculator(int row, int col, Model<bool> *pit_mask, Model<bool> *seen, Model<bool> *individual_pit_mask, std::vector<GeographicCoordinate> &brownfield_polygon);
ArrayCoordinate find_lowest_point_pit_lake(Model<bool> *individual_pit_mask);
double find_volume_pit_lake(double pit_area, int pit_depth);
double determine_circularity(Model<bool> *individual_pit_mask, ArrayCoordinate lowest_point, double pit_area);
void model_pit_lakes(PitCharacteristics &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, Model<bool> *seen_pl, Model<bool> *individual_pit_mask, Model<short> *DEM);
void model_depression(PitCharacteristics &pit, Model<bool> *pit_lake_mask, Model<bool> *depression_mask, Model<bool> *seen_d, Model<bool> *individual_pit_mask, Model<short> *DEM);

#endif