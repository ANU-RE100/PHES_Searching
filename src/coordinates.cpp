#include "model2D.h"
#include "phes_base.h"

GeographicCoordinate GeographicCoordinate_init(double latitude, double longitude)
{
	GeographicCoordinate geographic_coordinate;
	geographic_coordinate.lat = latitude;
	geographic_coordinate.lon = longitude;
	return geographic_coordinate;
}

GeographicCoordinate get_origin(GridSquare square, int border)
{
	GeographicCoordinate geographic_coordinate;
	geographic_coordinate.lat = square.lat+(1+((border-tile_overlap)+0.5)/(double((model_size-tile_overlap))));
	geographic_coordinate.lon = square.lon-((border+0.5)/(double((model_size-tile_overlap))));
	return geographic_coordinate;
}

ArrayCoordinate ArrayCoordinate_init(int row, int col, GeographicCoordinate origin) //Ditch at some point?
{
	ArrayCoordinate array_coordinate;
	array_coordinate.row = row;
	array_coordinate.col = col;
	array_coordinate.origin = origin;
	return array_coordinate;
}

ArrayCoordinateWithHeight ArrayCoordinateWithHeight_init(int row, int col, double h)
{
	ArrayCoordinateWithHeight array_coordinate;
	array_coordinate.row = (short)row;
	array_coordinate.col = (short)col;
	array_coordinate.h = h;
	return array_coordinate;
}

bool check_within(ArrayCoordinateWithHeight c, int shape[2])
{
	if(c.row>=0 && c.col>=0 && c.row<shape[0] && c.col<shape[1]){
        	return true;
	}
        return false;
}

bool check_within(ArrayCoordinate c, int shape[2])
{
	if(c.row>=0 && c.col>=0 && c.row<shape[0] && c.col<shape[1]){
        	return true;
	}
        return false;
}

bool check_within(GeographicCoordinate gc, GridSquare gs){
  return convert_to_int(FLOOR(gc.lat)) == gs.lat && convert_to_int(FLOOR(gc.lon)) == gs.lon;
}

bool check_strictly_within(ArrayCoordinate c, int shape[2])
{
	if(c.row>0 && c.col>0 && c.row<shape[0]-1 && c.col<shape[1]-1){
        	return true;
	}
        return false;
}


GridSquare GridSquare_init(int latitude, int longitude)
{
	GridSquare grid_square;
	grid_square.lat = latitude;
	grid_square.lon = longitude;
	return grid_square;
}

string str(GridSquare square)
{
	char buf[24];
	square.lon = (square.lon+180)%360-180;
	char c1 = (square.lat<0)?'s':'n';
	int lat = abs(square.lat);
	char c2 = (square.lon<0)?'w':'e';
	int lon = abs(square.lon);
	sprintf(buf, "%c%02d_%c%03d", c1, lat, c2, lon);
	string to_return(buf);
	return to_return;
}

string str_fabdem(GridSquare square)
{
	char buf[24];
	square.lon = (square.lon+180)%360-180;
	char c1 = (square.lat<0)?'S':'N';
	int lat = abs(square.lat);
	char c2 = (square.lon<0)?'W':'E';
	int lon = abs(square.lon);
	sprintf(buf, "%c%02d%c%03d", c1, lat, c2, lon);
	string to_return(buf);
	return to_return;
}

// area of single cell in ha
double find_area(ArrayCoordinate c)
{
	GeographicCoordinate p = convert_coordinates(c);
	return (0.0001*resolution*resolution)*COS(RADIANS(p.lat));
}


double find_distance(ArrayCoordinate c1, ArrayCoordinate c2)
{
	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return find_distance(p1, p2);
}

double find_distance(ArrayCoordinate c1, ArrayCoordinate c2, double coslat)
{
	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return find_distance(p1, p2, coslat);
}

double find_distance_sqd(ArrayCoordinate c1, ArrayCoordinate c2)
{
	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return find_distance_sqd(p1, p2);
}

double find_distance_sqd(ArrayCoordinate c1, ArrayCoordinate c2, double coslat)
{
  if(c1.origin.lat==c2.origin.lat && c1.origin.lon==c2.origin.lon)
    return (SQ(c2.row-c1.row)*coslat + SQ(c2.col-c1.col))*SQ(resolution*0.001);
	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return find_distance_sqd(p1, p2, coslat);
}

double find_distance(GeographicCoordinate c1, GeographicCoordinate c2)
{
	return SQRT(find_distance_sqd(c1, c2));
}

double find_distance(GeographicCoordinate c1, GeographicCoordinate c2, double coslat)
{
	return SQRT(find_distance_sqd(c1, c2, coslat));
}

double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2)
{
	return (SQ(c2.lat-c1.lat)+SQ((c2.lon-c1.lon)*COS(RADIANS(0.5*(c1.lat+c2.lat)))))*SQ((model_size-tile_overlap)*resolution*0.001);
}

double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2, double coslat)
{
	return (SQ(c2.lat-c1.lat)+SQ((c2.lon-c1.lon)*coslat))*SQ((model_size-tile_overlap)*resolution*0.001);
}

ArrayCoordinate convert_coordinates(GeographicCoordinate c, GeographicCoordinate origin)
{
	return ArrayCoordinate_init(convert_to_int((origin.lat-c.lat)*(model_size-tile_overlap)-0.5), convert_to_int((c.lon-origin.lon)*(model_size-tile_overlap)-0.5), origin);
}

ArrayCoordinate convert_coordinates(GeographicCoordinate c, GeographicCoordinate origin, double lat_res, double lon_res){
	return ArrayCoordinate_init(convert_to_int((c.lat-origin.lat)/lat_res-0.5), convert_to_int((c.lon-origin.lon)/lon_res-0.5), origin);
}

GeographicCoordinate convert_coordinates(ArrayCoordinate c, double offset)
{
	return GeographicCoordinate_init(c.origin.lat-(c.row+offset)/(double((model_size-tile_overlap))), c.origin.lon+(c.col+offset)/(double((model_size-tile_overlap))));
}

double find_orthogonal_nn_distance(ArrayCoordinate c1, ArrayCoordinate c2)
{
	if (c1.col == c2.col)
		return resolution;

	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return (COS(RADIANS(0.5*(p1.lat+p2.lat)))*resolution);
}

void pushback_non_duplicate_points(std::vector<ArrayCoordinate> &main_vector, std::vector<ArrayCoordinate> check_vector) {
	for(const auto& point2 : check_vector) {
        auto iterator = std::find_if(main_vector.begin(), main_vector.end(),
            [&point2](const ArrayCoordinate& point1) {
                return point1 == point2;
            });

        // If the point was not found in the first vector, add it
        if(iterator == main_vector.end()) {
            main_vector.push_back(point2);
        }
    }

	return;
}

