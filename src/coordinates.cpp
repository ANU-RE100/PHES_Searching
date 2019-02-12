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
	geographic_coordinate.lat = square.lat+(1+((border-1)+0.5)/3600.0);
	geographic_coordinate.lon = square.lon-((border+0.5)/3600.0);
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

bool flows_to(ArrayCoordinate c1, ArrayCoordinate c2, Model<char>* flow_directions) {
	return ( ( c1.row + directions[flow_directions->get(c1.row,c1.col)].row == c2.row ) &&
		 ( c1.col + directions[flow_directions->get(c1.row,c1.col)].col == c2.col ) );
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
	return (SQ(c2.lat-c1.lat)+SQ((c2.lon-c1.lon)*COS(RADIANS(0.5*(c1.lat+c2.lat)))))*SQ(3600*resolution*0.001);
}

double find_distance_sqd(GeographicCoordinate c1, GeographicCoordinate c2, double coslat)
{
	return (SQ(c2.lat-c1.lat)+SQ((c2.lon-c1.lon)*coslat))*SQ(3600*resolution*0.001);
}

ArrayCoordinate convert_coordinates(GeographicCoordinate c, GeographicCoordinate origin)
{
	return ArrayCoordinate_init(convert_to_int((origin.lat-c.lat)*3600-0.5), convert_to_int((c.lon-origin.lon)*3600-0.5), origin);
}

ArrayCoordinate convert_coordinates(GeographicCoordinate c, GeographicCoordinate origin, double lat_res, double lon_res){
	return ArrayCoordinate_init(convert_to_int((c.lat-origin.lat)/lat_res-0.5), convert_to_int((c.lon-origin.lon)/lon_res-0.5), origin);
}

GeographicCoordinate convert_coordinates(ArrayCoordinate c)
{
	return GeographicCoordinate_init(c.origin.lat-(c.row+0.5)/3600.0, c.origin.lon+(c.col+0.5)/3600.0);
}


double find_orthogonal_nn_distance(ArrayCoordinate c1, ArrayCoordinate c2)
{
	if (c1.col == c2.col)
		return resolution;

	GeographicCoordinate p1 = convert_coordinates(c1);
	GeographicCoordinate p2 = convert_coordinates(c2);
	return (COS(RADIANS(0.5*(p1.lat+p2.lat)))*resolution);
}