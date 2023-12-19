#ifndef JSON_H
#define JSON_H

#include "model2D.h"
#include <vector>
#include <string>

void clean_geojson(string &geojson_raw);
std::vector<std::string> split(const std::string &s, const std::string &delim);
std::vector<std::string> split(const std::string &s, char delim);
GeographicCoordinate parseCoordinates(const std::string &str);
vector<GeographicCoordinate> parseList(const std::string &str);
vector<vector<GeographicCoordinate>> parseList2D(const std::string &str);
vector<vector<vector<GeographicCoordinate>>> parseList3D(const std::string &str);

#endif