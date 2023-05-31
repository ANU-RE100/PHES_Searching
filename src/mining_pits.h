#ifndef MINING_PITS_H
#define MINING_PITS_H

#include <bits/stdc++.h>
#include "model2D.h"
#include "coordinates.h"

using namespace std;

void depression_volume_finding(Model<short>* DEM);
double pit_area_calculator(int row, int col, Model<bool> *pit_mask, Model<bool> *seen, Model<bool> *individual_pit_mask, std::vector<GeographicCoordinate> &brownfield_polygon);
ArrayCoordinate find_lowest_point_pit_lake(Model<bool> *individual_pit_mask);
double find_volume_pit_lake(double pit_area, int pit_depth);
void find_depression_attributes(Model<bool> *individual_pit_mask, Model<short> *DEM, ArrayCoordinate &lowest_point, double &pit_volume, int &depression_elevation, vector<GeographicCoordinate> brownfield_polygon);
double determine_circularity(Model<bool> *individual_pit_mask, ArrayCoordinate lowest_point, double pit_area);

#endif