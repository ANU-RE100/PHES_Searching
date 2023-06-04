#ifndef POLYGON_H
#define POLYGON_H

#include "phes_base.h"

vector<double> find_polygon_intersections(int row, vector<GeographicCoordinate> &polygon, Model<bool>* filter);
void polygon_to_raster(vector<GeographicCoordinate> &polygon, Model<bool>* raster);
void read_shp_filter(string filename, Model<bool>* filter);
ArrayCoordinate find_edge(ArrayCoordinate seed_point, Model<bool> *mask);

#endif
