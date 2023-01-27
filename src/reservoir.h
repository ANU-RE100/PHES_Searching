#ifndef RESERVOIR_H
#define RESERVOIR_H

#include "phes_base.h"

class RoughReservoir{
  public:
    string identifier;
    bool brownfield = false;
    bool ocean;
    bool pit;
    bool turkey;
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
    int bottom_elevation;

    RoughReservoir() {};
    virtual ~RoughReservoir() = default;
    RoughReservoir(ArrayCoordinate pour_point, int elevation): brownfield(false),
    pit(false), elevation(elevation), pour_point(pour_point), max_dam_height(max_wall_height), bottom_elevation(elevation){
      GeographicCoordinate geo_coordinate = convert_coordinates(pour_point);
      this->latitude = geo_coordinate.lat;
      this->longitude = geo_coordinate.lon;

    }

    bool operator<(const RoughReservoir &o) const
        {
      return elevation > o.elevation;
        }
  private:
};

class RoughGreenfieldReservoir : public RoughReservoir {
public:
  vector<array<ArrayCoordinate, directions.size()>> shape_bound;

  RoughGreenfieldReservoir(RoughReservoir r)
      : RoughReservoir(r) {
    for (uint ih = 0; ih < dam_wall_heights.size(); ih++) {
      array<ArrayCoordinate, directions.size()> temp_array;
      for (uint idir = 0; idir < directions.size(); idir++) {
        temp_array[idir].row = pour_point.row;
        temp_array[idir].col = pour_point.col;
      }
      this->shape_bound.push_back(temp_array);
    }
  }
};

class RoughBfieldReservoir : public RoughReservoir {
public:
  vector<ArrayCoordinate> shape_bound;
  RoughBfieldReservoir() {};
  RoughBfieldReservoir(RoughReservoir r) : RoughReservoir(r) {}
};

struct ExistingReservoir {
  string identifier;
  double latitude;
  double longitude;
  int elevation;
  int bottom_elevation;
  double volume;
  vector<GeographicCoordinate> polygon;
};

struct AltitudeVolumePair {
  int altitude;
  double volume;
  bool operator<(const AltitudeVolumePair &o) const { return altitude < o.altitude; }
};

struct ExistingPit {
  ExistingReservoir reservoir;
  vector<AltitudeVolumePair> volumes;
};

struct Reservoir {
  string identifier;
  bool brownfield;
  bool pit;
  bool ocean;
  double latitude;
  double longitude;
  int elevation;
  int bottom_elevation;
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
  bool operator<(const Reservoir &o) const { return elevation > o.elevation; }
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
  double energy_capacity;
  int storage_time;
  int head;
  int non_overlap;
  string country;
  bool operator<(const Pair &o) const { return FOM < o.FOM; }
};

void update_reservoir_boundary(vector<array<ArrayCoordinate, directions.size()>> &dam_shape_bounds,
                               ArrayCoordinate point, int elevation_above_pp);
void update_reservoir_boundary(array<ArrayCoordinate, directions.size()> &dam_shape_bounds,
                               ArrayCoordinate point);
Reservoir Reservoir_init(ArrayCoordinate pour_point, int elevation);
ExistingReservoir ExistingReservoir_init(string identifier, double latitude, double longitude,
                                         int elevation, double volume);
ExistingPit ExistingPit_init(ExistingReservoir reservoir);
GridSquare get_square_coordinate(ExistingReservoir reservoir);

#endif
