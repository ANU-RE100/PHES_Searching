#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "phes_base.h"

struct RoughReservoir{
	string identifier;
	double latitude;
	double longitude;
	int elevation;
	ArrayCoordinate pour_point;
	vector<double> volumes;
	vector<double> dam_volumes;
	vector<double> areas;
	vector<double> water_rocks;
	double watershed_area;
	double max_dam_height;
	vector<array<ArrayCoordinate, directions.size()> > shape_bound;
	bool operator<(const RoughReservoir &o) const
	    {
		return elevation > o.elevation;
	    }
};

struct Reservoir{
	string identifier;
	double latitude;
	double longitude;
	int elevation;
	ArrayCoordinate pour_point;
	double volume;
	double dam_volume;
	double dam_length;
	double area;
	double water_rock;
	double watershed_area;
	double average_water_depth;
	double dam_height;
	double max_dam_height;
	string country;
	array<ArrayCoordinate, directions.size()> shape_bound;
	bool operator<(const Reservoir &o) const
	    {
		return elevation > o.elevation;
	    }
};

struct Pair {
	Reservoir upper;
	Reservoir lower;
	string identifier;
	double distance;
	double pp_distance;
	double slope;
	double required_volume;
	double volume;
	double FOM;
	char category;
	double water_rock;
	int energy_capacity;
	int storage_time;
	int head;
	int non_overlap;
	string country;
	bool operator<(const Pair &o) const
	    {
		return FOM < o.FOM;
	    }
};

void update_reservoir_boundary(vector<array<ArrayCoordinate, directions.size()> > &dam_shape_bounds, ArrayCoordinate point, int elevation_above_pp);
void update_reservoir_boundary(array<ArrayCoordinate, directions.size()> &dam_shape_bounds, ArrayCoordinate point);
RoughReservoir RoughReservoir_init(ArrayCoordinate pour_point, int elevation);
Reservoir Reservoir_init(ArrayCoordinate pour_point, int elevation);

#endif
