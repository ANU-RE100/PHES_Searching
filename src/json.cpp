#include "model2D.h"
#include <vector>
#include <string>

void clean_geojson(string &geojson_raw) {
  std::string start_ss = "{coordinates: ";
  std::string end_ss1 = ", type: Polygon}";
  std::string end_ss2 = ", type: MultiPolygon}";

  std::string::size_type start_i = geojson_raw.find(start_ss);
  if (start_i != std::string::npos)
    geojson_raw.erase(start_i, start_ss.length());

  std::string::size_type end1_i = geojson_raw.find(end_ss1);
  if (end1_i != std::string::npos)
    geojson_raw.erase(end1_i, end_ss1.length());

  std::string::size_type end2_i = geojson_raw.find(end_ss2);  
  if (end2_i != std::string::npos)
    geojson_raw.erase(end2_i, end_ss2.length());

  geojson_raw.erase(remove(geojson_raw.begin(), geojson_raw.end(), ' '), geojson_raw.end());

  return;
}

std::vector<std::string> split(const std::string &s, const std::string &delim) {
    std::vector<std::string> elems;
    size_t pos = 0;
    size_t lastPos = 0;
    while((pos = s.find(delim, lastPos)) != std::string::npos) {
        elems.push_back(s.substr(lastPos, pos - lastPos)+delim);
        lastPos = pos + delim.size();
    }
    // Add the last token
    elems.push_back(s.substr(lastPos));
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


GeographicCoordinate parseCoordinates(const std::string &str) {
    std::string cleaned;
    if(str[0] == ',')
        cleaned = str.substr(2, str.size() - 3);  // Remove outer brackets and comma
    else
        cleaned = str.substr(1, str.size() - 2);  // Remove outer brackets
    std::vector<std::string> parts = split(cleaned, ',');
    return { std::stod(parts[1]), std::stod(parts[0]) };
}

vector<GeographicCoordinate> parseList(const std::string &str) {
    std::string cleaned;
    if(str[0] == ',')
        cleaned = str.substr(2, str.size() - 2);  // Remove outer brackets and comma
    else
        cleaned = str.substr(1, str.size() - 2);  // Remove outer brackets
    std::vector<std::string> parts = split(cleaned, ']');
    vector<GeographicCoordinate> result;
    for (const std::string &part : parts) {
        if (part.empty()) {  // Skip commas and empty strings
            continue;
        }
        result.push_back(parseCoordinates(part + "]"));
    }
    return result;
}

vector<vector<GeographicCoordinate>> parseList2D(const std::string &str) {
    std::string cleaned;
    if(str[0] == ',')
        cleaned = str.substr(2, str.size() - 2);  // Remove outer brackets and comma
    else
        cleaned = str.substr(1, str.size() - 2);  // Remove outer brackets
    std::vector<std::string> parts = split(cleaned, "]]");
    vector<vector<GeographicCoordinate>> result;
    for (const std::string &part : parts) {
        if (part.empty()) {  // Skip non-lists
            continue;
        }
        result.push_back(parseList(part));
    }
    return result;
}

vector<vector<vector<GeographicCoordinate>>> parseList3D(const std::string &str) {
    std::string cleaned;
    if(str[0] == ',')
        cleaned = str.substr(2, str.size() - 2);  // Remove outer brackets and comma
    else
        cleaned = str.substr(1, str.size() - 2);  // Remove outer brackets
    std::vector<std::string> parts = split(cleaned, "]]]");
    vector<vector<vector<GeographicCoordinate>>> result;
    for (const std::string &part : parts) {
        if (part.empty()) {  // Skip non-lists
            continue;
        }
        result.push_back(parseList2D(part));
    }
    return result;
}
