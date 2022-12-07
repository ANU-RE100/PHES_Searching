#ifndef COORDINATES_H
#define COORDINATES_H

#include <string>

using namespace std;

struct GeographicCoordinate;
struct ArrayCoordinate;
template <class T> class Model;

struct GridSquare {
  int lat, lon;
};

struct ArrayCoordinateWithHeight {
  short row, col;
  double h;
  bool operator<(const ArrayCoordinateWithHeight &o) const { return h > o.h; }
  bool operator==(const ArrayCoordinateWithHeight &o) const { return (row == o.row) && (col == o.col); }
};

GeographicCoordinate GeographicCoordinate_init(double latitude,
                                               double longitude);
ArrayCoordinate ArrayCoordinate_init(int row, int col);
ArrayCoordinate ArrayCoordinate_init(int row, int col,
                                     GeographicCoordinate origin);
GridSquare GridSquare_init(int latitude, int longitude);
ArrayCoordinateWithHeight ArrayCoordinateWithHeight_init(int row, int col,
                                                         double h);
GeographicCoordinate get_origin(GridSquare square, int border);
bool check_within(ArrayCoordinateWithHeight c, int shape[2]);
bool check_within(ArrayCoordinate c, int shape[2]);
bool check_within(GeographicCoordinate gc, GridSquare gs);
string str(GridSquare square);
double find_area(ArrayCoordinate c);
double find_distance(ArrayCoordinate c1, ArrayCoordinate c2);
double find_distance(ArrayCoordinate c1, ArrayCoordinate c2, double coslat);
double find_distance_sqd(ArrayCoordinate c1, ArrayCoordinate c2);
double find_distance_sqd(ArrayCoordinate c1, ArrayCoordinate c2, double coslat);
double find_distance(GeographicCoordinate c1, GeographicCoordinate c2);
double find_distance(GeographicCoordinate c1, GeographicCoordinate c2,
                     double coslat);
double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2);
double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2,
                         double coslat);
GeographicCoordinate convert_coordinates(ArrayCoordinate c, double offset=0.5);
ArrayCoordinate convert_coordinates(GeographicCoordinate c,
                                    GeographicCoordinate origin);
ArrayCoordinate convert_coordinates(GeographicCoordinate c,
                                    GeographicCoordinate origin, double lat_res,
                                    double lon_res);
double find_orthogonal_nn_distance(ArrayCoordinate c1, ArrayCoordinate c2);

#endif
