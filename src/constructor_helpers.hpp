#ifndef CONSTRUCTOR_HELPER_H
#define CONSTRUCTOR_HELPER_H

#include "phes_base.h"
#include "kml.h"

vector<double> find_polygon_intersections(double lat, vector<GeographicCoordinate> &polygon);
bool check_within(GeographicCoordinate point, vector<vector<GeographicCoordinate>> polygons);
vector<vector<vector<GeographicCoordinate>>> read_countries(string filename, vector<string>& country_names);
ArrayCoordinate* get_adjacent_cells(ArrayCoordinate point1, ArrayCoordinate point2);
bool is_edge(ArrayCoordinate point1, ArrayCoordinate point2, Model<char>* model, ArrayCoordinate offset, int threshold);
bool is_dam_wall(ArrayCoordinate point1, ArrayCoordinate point2, Model<short>* DEM, ArrayCoordinate offset, double wall_elevation);

vector<ArrayCoordinate> convert_to_polygon(Model<char>* model, ArrayCoordinate offset, ArrayCoordinate pour_point, int threshold);
vector<ArrayCoordinate> convert_to_polygon(Model<bool>* model, ArrayCoordinate offset, ArrayCoordinate pour_point, int threshold);
vector<GeographicCoordinate> convert_poly(vector<ArrayCoordinate> polygon);
vector<GeographicCoordinate> corner_cut_poly(vector<GeographicCoordinate> polygon);
vector<GeographicCoordinate> compress_poly(vector<GeographicCoordinate> polygon);
string str(vector<GeographicCoordinate> polygon, double elevation);
bool model_reservoir(Reservoir *reservoir,
                     Reservoir_KML_Coordinates *coordinates, Model<bool> *seen,
                     bool *non_overlap, vector<ArrayCoordinate> *used_points,
                     BigModel big_model, Model<char> *full_cur_model,
                     vector<vector<vector<GeographicCoordinate>>> &countries,
                     vector<string> &country_names);
bool model_bulk_pit(Reservoir *reservoir, Reservoir_KML_Coordinates *coordinates,
                     vector<vector<vector<GeographicCoordinate>>> &countries,
                     vector<string> &country_names, std::vector<PitCharacteristics> pit_shapes);

#endif
