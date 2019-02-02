#include <bits/stdc++.h>
using namespace std;

#include "phes_base.h"

void update_reservoir_boundary(vector<array<ArrayCoordinate, directions.size()> > &dam_shape_bounds, ArrayCoordinate point, int elevation_above_pp)
{
	for (uint ih=0; ih<dam_wall_heights.size(); ih++) {
		int dam_height =  dam_wall_heights[ih];
		if (dam_height >= elevation_above_pp)
			for (uint i = 0; i<directions.size(); i++) {
				if ( (directions[i].row*point.row + directions[i].col * point.col) >
				     (directions[i].row*dam_shape_bounds[ih][i].row + directions[i].col*dam_shape_bounds[ih][i].col) ){
					dam_shape_bounds[ih][i] = point;
				}
			}
	}
}

void update_reservoir_boundary(array<ArrayCoordinate, directions.size()> &dam_shape_bounds, ArrayCoordinate point)
{
	for (uint i = 0; i<directions.size(); i++) {
		if ( (directions[i].row*point.row + directions[i].col * point.col) >
		     (directions[i].row*dam_shape_bounds[i].row + directions[i].col*dam_shape_bounds[i].col) ){
			dam_shape_bounds[i].row = point.row;
			dam_shape_bounds[i].col = point.col;
		}
	}
}

RoughReservoir RoughReservoir_init(ArrayCoordinate pour_point, int elevation)
{
	RoughReservoir reservoir;
	reservoir.elevation = elevation;
	GeographicCoordinate geo_coordinate = convert_coordinates(pour_point);
	reservoir.latitude = geo_coordinate.lat;
	reservoir.longitude = geo_coordinate.lon;
	reservoir.pour_point = pour_point;
	reservoir.max_dam_height = max_wall_height;

	// initialize bounds
	for (uint ih=0; ih < dam_wall_heights.size(); ih++) {
		array<ArrayCoordinate, directions.size()> temp_array;
		for (uint idir=0; idir < directions.size(); idir++) {
			temp_array[idir].row = -100000 * directions[idir].row;
			temp_array[idir].col = -100000 * directions[idir].col;
		}
		reservoir.shape_bound.push_back(temp_array);
	}

	return reservoir;
}

Reservoir Reservoir_init(ArrayCoordinate pour_point, int elevation)
{
	Reservoir reservoir;
	reservoir.elevation = elevation;
	GeographicCoordinate geo_coordinate = convert_coordinates(pour_point);
	reservoir.latitude = geo_coordinate.lat;
	reservoir.longitude = geo_coordinate.lon;
	reservoir.pour_point = pour_point;
	for (uint idir=0; idir < directions.size(); idir++) {
		reservoir.shape_bound[idir] = pour_point;
	}
	return reservoir;
}
