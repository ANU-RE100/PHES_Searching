#ifndef KML_H
#define KML_H

#include "phes_base.h"

struct KML_Holder{
	vector<string> uppers;
	vector<string> lowers;
	vector<string> upper_dams;
	vector<string> lower_dams;
	vector<string> lines;
	vector<string> points;
};

struct Reservoir_KML_Coordinates{
	string reservoir;
	vector<string> dam;
	bool is_turkeys_nest;
};

struct Pair_KML{
	Reservoir_KML_Coordinates upper;
	Reservoir_KML_Coordinates lower;
	string point;
	string line;
};

string output_kml(KML_Holder* kml_holder, string square, Test test);
void update_kml_holder(KML_Holder* kml_holder, Pair* pair, Pair_KML* pair_kml, bool keep_upper, bool keep_lower);
string get_reservoir_geometry(Reservoir_KML_Coordinates coordinates);
string get_dam_geometry(Reservoir_KML_Coordinates coordinates);
string get_dam_kml(Reservoir* reservoir, Reservoir_KML_Coordinates coordinates);
extern string kml_start;
extern string kml_end;
// void write_fusion_csv_header(FILE *csv_file);
// void write_fusion_csv(FILE *csv_file, Pair *pair, Pair_KML* pair_kml);

#endif
