#include "phes_base.h"

int display = false;

struct Reservoir_KML_Coordinates {
  string reservoir;
  vector<string> dam;
  bool is_turkeys_nest;
};

// Returns a tuple with the two cells adjacent to an edge defined by two points
ArrayCoordinate *get_adjacent_cells(ArrayCoordinate point1, ArrayCoordinate point2) {
  double average_row = (point1.row + point2.row) / 2.0;
  double average_col = (point1.col + point2.col) / 2.0;

  ArrayCoordinate *points = (ArrayCoordinate *)malloc(2 * sizeof(ArrayCoordinate));

  if (average_row - (int)average_row > 0.9 || average_row - (int)average_row < 0.1) {
    points[0] = ArrayCoordinate_init((int)(average_row + EPS), (int)(average_col - 0.5 + EPS),
                                     point1.origin);
    points[1] = ArrayCoordinate_init((int)(average_row - 1 + EPS), (int)(average_col - 0.5 + EPS),
                                     point1.origin);
  } else {
    points[0] = ArrayCoordinate_init((int)(average_row - 0.5 + EPS), (int)(average_col + EPS),
                                     point1.origin);
    points[1] = ArrayCoordinate_init((int)(average_row - 0.5 + EPS), (int)(average_col - 1 + EPS),
                                     point1.origin);
  }
  return points;
}

// Determines if two points give the edge of a reservoir given its raster model
bool is_edge(ArrayCoordinate point1, ArrayCoordinate point2, Model<char> *model,
             ArrayCoordinate offset, int threshold) {

  if (point1.row + offset.row < 0 || point1.col + offset.col < 0 ||
      point1.row + offset.row > model->nrows() || point1.col + offset.col > model->ncols())
    return false;
  if (point2.row + offset.row < 0 || point2.col + offset.col < 0 ||
      point2.row + offset.row > model->nrows() || point2.col + offset.col > model->ncols())
    return false;

  ArrayCoordinate *to_check = get_adjacent_cells(point1, point2);

  if ((!model->check_within(to_check[0].row + offset.row, to_check[0].col + offset.col) &&
       model->get(to_check[1].row + offset.row, to_check[1].col + offset.col) >= threshold) ||
      (!model->check_within(to_check[1].row + offset.row, to_check[1].col + offset.col) &&
       model->get(to_check[0].row + offset.row, to_check[0].col + offset.col) >= threshold))
    return true;

  if (model->get(to_check[0].row + offset.row, to_check[0].col + offset.col) < threshold &&
      model->get(to_check[1].row + offset.row, to_check[1].col + offset.col) >= threshold)
    return true;
  if (model->get(to_check[1].row + offset.row, to_check[1].col + offset.col) < threshold &&
      model->get(to_check[0].row + offset.row, to_check[0].col + offset.col) >= threshold)
    return true;

  return false;
}

// Determines if an edge of a reservoir between two points requires a dam wall
bool is_dam_wall(ArrayCoordinate point1, ArrayCoordinate point2, Model<short> *DEM,
                 ArrayCoordinate offset, double wall_elevation) {
  if (point1.row < 0 || point1.col < 0 || point1.row > DEM->nrows() || point1.col > DEM->ncols())
    return false;
  if (point2.row < 0 || point2.col < 0 || point2.row > DEM->nrows() || point2.col > DEM->ncols())
    return false;

  double average_row = (point1.row + point2.row) / 2.0;
  double average_col = (point1.col + point2.col) / 2.0;

  ArrayCoordinate *to_check = get_adjacent_cells(point1, point2);

  if (average_row - (int)average_row > 0.9 || average_row - (int)average_row < 0.1) {
    to_check[0] = ArrayCoordinate_init((int)(average_row + EPS), (int)(average_col - 0.5 + EPS),
                                       point1.origin);
    to_check[1] = ArrayCoordinate_init((int)(average_row - 1 + EPS), (int)(average_col - 0.5 + EPS),
                                       point1.origin);
  } else {
    to_check[0] = ArrayCoordinate_init((int)(average_row - 0.5 + EPS), (int)(average_col + EPS),
                                       point1.origin);
    to_check[1] = ArrayCoordinate_init((int)(average_row - 0.5 + EPS), (int)(average_col - 1 + EPS),
                                       point1.origin);
  }

  if (!DEM->check_within(to_check[0].row, to_check[0].col) ||
      !DEM->check_within(to_check[1].row, to_check[1].col))
    return false;

  if (DEM->get(to_check[0].row + offset.row, to_check[0].col + offset.col) < wall_elevation and
      DEM->get(to_check[1].row + offset.row, to_check[1].col + offset.col) < wall_elevation)
    return true;
  return false;
}

// Converts a raster model to a polygon given a raster model and a point on the interior edge of the
// polygon
vector<ArrayCoordinate> convert_to_polygon(Model<char> *model, ArrayCoordinate offset,
                                           ArrayCoordinate pour_point, int threshold) {

  vector<ArrayCoordinate> to_return;
  int dir_def[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
  int dir_to_do[4][3] = {{2, 0, 3}, {3, 1, 2}, {1, 2, 0}, {0, 3, 1}};

  int tests[] = {1, 2, 0, 3};
  ArrayCoordinate test_coordinates[] = {
      ArrayCoordinate_init(pour_point.row + 1, pour_point.col + 1, pour_point.origin),
      ArrayCoordinate_init(pour_point.row + 1, pour_point.col, pour_point.origin),
      ArrayCoordinate_init(pour_point.row, pour_point.col, pour_point.origin),
      ArrayCoordinate_init(pour_point.row, pour_point.col + 1, pour_point.origin)};

  for (int i = 0; i < 4; i++) {
    vector<ArrayCoordinate> temp_to_return;
    int last_dir = tests[i];
    ArrayCoordinate initial = test_coordinates[i];
    ArrayCoordinate last = initial;
    while (true) {
      for (int id = 0; id < 3; id++) {
        int d = dir_to_do[last_dir][id];
        ArrayCoordinate next = ArrayCoordinate_init(last.row + dir_def[d][0],
                                                    last.col + dir_def[d][1], pour_point.origin);
        if (is_edge(last, next, model, offset, threshold)) {
          temp_to_return.push_back(next);
          last = next;
          last_dir = d;
          break;
        }
      }
      if (last.row == initial.row && last.col == initial.col)
        break;
    }
    if (temp_to_return.size() > to_return.size()) {
      to_return.clear();
      for (uint j = 0; j < temp_to_return.size(); j++) {
        to_return.push_back(temp_to_return[j]);
      }
    }
  }
  return to_return;
}

vector<GeographicCoordinate> convert_poly(vector<ArrayCoordinate> polygon) {
  vector<GeographicCoordinate> to_return;
  for (uint i = 0; i < polygon.size(); i++) {
    to_return.push_back(convert_coordinates(polygon[i]));
  }
  return to_return;
}

// "Smooths" a polygon by cutting the corners
vector<GeographicCoordinate> corner_cut_poly(vector<GeographicCoordinate> polygon) {
  vector<GeographicCoordinate> to_return;
  for (uint i = 0; i < polygon.size(); i++) {
    to_return.push_back(GeographicCoordinate_init(
        (polygon[i].lat * 2 + polygon[(i + 1) % polygon.size()].lat * 2) / 4.0,
        (polygon[i].lon * 2 + polygon[(i + 1) % polygon.size()].lon * 2) / 4.0));
  }
  to_return.push_back(to_return[0]);
  return to_return;
}

// Compresses a polygon by removing collinear points
vector<GeographicCoordinate> compress_poly(vector<GeographicCoordinate> polygon) {
  vector<GeographicCoordinate> to_return;
  to_return.push_back(polygon[0]);
  for (uint i = 1; i < polygon.size() - 1; i++) {
    if ((polygon[i + 1].lon - polygon[i].lon) / (polygon[i + 1].lat - polygon[i].lat) <
            (polygon[i].lon - polygon[i - 1].lon) / (polygon[i].lat - polygon[i - 1].lat) - EPS ||
        (polygon[i + 1].lon - polygon[i].lon) / (polygon[i + 1].lat - polygon[i].lat) >
            (polygon[i].lon - polygon[i - 1].lon) / (polygon[i].lat - polygon[i - 1].lat) + EPS) {
      to_return.push_back(polygon[i]);
    }
  }
  to_return.push_back(polygon[polygon.size() - 1]);
  return to_return;
}

string str(vector<GeographicCoordinate> polygon, double elevation) {
  string to_return = "";
  for (uint i = 0; i < polygon.size(); i++) {
    to_return +=
        dtos(polygon[i].lon, 5) + "," + dtos(polygon[i].lat, 5) + "," + dtos(elevation, 1) + " ";
  }
  return to_return;
}

string kml_start = "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
                   "  <Document id=\"Layers\">\n"
                   "    <Style id=\"reservoir_style\">\n"
                   "      <LineStyle>\n"
                   "        <width>0</width>\n"
                   "      </LineStyle>\n"
                   "    </Style>\n"
                   "    <Style id=\"wall_style\">\n"
                   "      <LineStyle>\n"
                   "        <width>0</width>\n"
                   "      </LineStyle>\n"
                   "      <PolyStyle>\n"
                   "        <color>ff999999</color>\n"
                   "      </PolyStyle>\n"
                   "    </Style>\n"
                   "    <Style id=\"pipe_style\">\n"
                   "      <LineStyle>\n"
                   "        <color>88ffffff</color>\n"
                   "        <width>2</width>\n"
                   "      </LineStyle>\n"
                   "    </Style>\n"
                   "    <Style id=\"pin_style\">\n"
                   "      <LabelStyle>\n"
                   "        <color>00ffffff</color>\n"
                   "        <scale>0.6</scale>\n"
                   "      </LabelStyle>\n"
                   "    </Style>\n";

string kml_end = "  </Document>\n"
                 "</kml>\n";
string get_html(Reservoir *reservoir) {
  return "              <html>\n"
         "              <head><META http-equiv=\"Content-Type\" content=\"text/html\"><meta "
         "http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\"></head>\n"
         "              <body style=\"margin:0px 0px 0px 0px;overflow:auto;background:#FFFFFF;\">\n"
         "              <table "
         "style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-"
         "collapse:collapse;padding:3px 3px 3px 3px\">\n"
         "              <tr style=\"text-align:center;font-weight:bold;background:#9CBCE2\"><td>" +
         reservoir->identifier +
         "</td></tr>\n"
         "              <tr><td>\n"
         "              <table "
         "style=\"font-family:Arial,Verdana,Times;font-size:12px;text-align:left;width:100%;border-"
         "spacing:0px; padding:3px 3px 3px 3px\">\n"
         "              <tr><td>Reservoir Ref.</td><td>" +
         reservoir->identifier +
         "</td></tr>\n"
         "              <tr bgcolor=\"#D4E4F3\"><td>Elevation</td><td>" +
         to_string(reservoir->elevation) +
         "</td></tr>\n"
         "              <tr><td>Latitude</td><td>" +
         dtos(reservoir->latitude, 4) +
         "</td></tr>\n"
         "              <tr bgcolor=\"#D4E4F3\"><td>Longitude</td><td>" +
         dtos(reservoir->longitude, 4) +
         "</td></tr>\n"
         "              <tr><td>Area (ha)</td><td>" +
         dtos(reservoir->area, 0) +
         "</td></tr>\n"
         "              <tr bgcolor=\"#D4E4F3\"><td>Volume (GL)</td><td>" +
         dtos(reservoir->volume, 1) +
         "</td></tr>\n"
         "              <tr><td>Dam Wall Height (m)</td><td>" +
         dtos(reservoir->dam_height, 1) +
         "</td></tr>\n"
         "              <tr bgcolor=\"#D4E4F3\"><td>Dam Length (m)</td><td>" +
         dtos(reservoir->dam_length, 0) +
         "</td></tr>\n"
         "              <tr><td>Dam Volume (GL)</td><td>" +
         dtos(reservoir->dam_volume, 1) +
         "</td></tr>\n"
         "              <tr bgcolor=\"#D4E4F3\"><td>Water/Rock Ratio</td><td>" +
         dtos(reservoir->water_rock, 1) +
         "</td></tr>\n"
         "              </table></td></tr></table></body></html>\n";
}
string get_reservoir_geometry(Reservoir_KML_Coordinates coordinates) {
  return "        <MultiGeometry>\n"
         "          <Polygon>\n"
         "            <extrude>1</extrude>\n"
         "            <altitudeMode>absolute</altitudeMode>\n"
         "            <outerBoundaryIs><LinearRing>\n"
         "              <coordinates>" +
         coordinates.reservoir +
         "</coordinates>\n"
         "            </LinearRing></outerBoundaryIs>\n"
         "          </Polygon>\n"
         "        </MultiGeometry>\n";
}
string get_dam_geometry(Reservoir_KML_Coordinates coordinates) {
  string to_return = "        <MultiGeometry>\n";
  if (coordinates.is_turkeys_nest) {
    to_return += "          <Polygon>\n"
                 "            <extrude>1</extrude>\n"
                 "            <altitudeMode>absolute</altitudeMode>\n"
                 "            <outerBoundaryIs><LinearRing>\n"
                 "              <coordinates>" +
                 coordinates.dam[0] +
                 "</coordinates>\n"
                 "            </LinearRing></outerBoundaryIs>\n"
                 "            <innerBoundaryIs><LinearRing>\n"
                 "              <coordinates>" +
                 coordinates.reservoir +
                 "</coordinates>\n"
                 "            </LinearRing></innerBoundaryIs>\n"
                 "          </Polygon>\n";
  } else {
    for (uint i = 0; i < coordinates.dam.size(); i++) {
      to_return += "          <Polygon>\n"
                   "            <extrude>1</extrude>\n"
                   "            <altitudeMode>absolute</altitudeMode>\n"
                   "            <outerBoundaryIs><LinearRing>\n"
                   "              <coordinates>" +
                   coordinates.dam[i] +
                   "</coordinates>\n"
                   "            </LinearRing></outerBoundaryIs>\n"
                   "          </Polygon>\n";
    }
  }
  return to_return + "        </MultiGeometry>\n";
}
string get_reservoir_kml(Reservoir *reservoir, string colour,
                         Reservoir_KML_Coordinates coordinates) {
  return "      <Placemark>\n"
         "        <name><![CDATA[" +
         reservoir->identifier +
         "]]></name>\n"
         "        <styleUrl>#reservoir_style</styleUrl>\n"
         "        <Style>\n"
         "          <PolyStyle>\n"
         "            <color>" +
         colour +
         "</color>\n"
         "          </PolyStyle>\n"
         "        </Style>\n" +
         get_reservoir_geometry(coordinates) +
         "        <description>\n"
         "           <![CDATA[\n" +
         get_html(reservoir) +
         "           ]]>\n"
         "        </description>\n"
         "      </Placemark>\n";
}

string get_dam_kml(Reservoir *reservoir, Reservoir_KML_Coordinates coordinates) {
  return "      <Placemark>\n"
         "        <name><![CDATA[" +
         reservoir->identifier +
         " Dam]]></name>\n"
         "        <styleUrl>#wall_style</styleUrl>\n" +
         get_dam_geometry(coordinates) + "      </Placemark>\n";
}

string output_kml(Reservoir *reservoir, Reservoir_KML_Coordinates coordinates) {
  string to_return;
  to_return += kml_start;
  to_return += get_reservoir_kml(reservoir, upper_colour, coordinates);
  to_return += get_dam_kml(reservoir, coordinates);
  to_return += kml_end;
  return to_return;
}

int main(int nargs, char **argv) {
  GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));
  display = true;

  printf("Reservoir constructor started for %s\n", convert_string(str(square_coordinate)));

  GDALAllRegister();
  unsigned long t_usec = walltime_usec();
  parse_variables(convert_string(file_storage_location + "variables"));

  BigModel big_model = BigModel_init(square_coordinate);
  Model<char> *full_cur_model =
      new Model<char>(big_model.DEM->nrows(), big_model.DEM->ncols(), MODEL_SET_ZERO);

  vector<RoughReservoir> reservoirs = read_rough_reservoir_data(
      convert_string(file_storage_location + "processing_files/reservoirs/" +
                     str(square_coordinate) + "_reservoirs_data.csv"));
  if (display)
    printf("Read in %zu reservoirs\n", reservoirs.size());

  string rs(argv[3]);

  for (uint i = 0; i < reservoirs.size(); i++) {
    if (reservoirs[i].identifier == rs) {

      ofstream kml_file(
          convert_string(file_storage_location + "output/" + reservoirs[i].identifier + ".kml"),
          ios::out);

      Reservoir *reservoir = new Reservoir();
      reservoir->identifier = reservoirs[i].identifier;
      reservoir->latitude = reservoirs[i].latitude;
      reservoir->longitude = reservoirs[i].longitude;
      reservoir->elevation = reservoirs[i].elevation;
      reservoir->dam_height = atoi(argv[4]);
      Reservoir_KML_Coordinates *coordinates = new Reservoir_KML_Coordinates();

      Model<short> *DEM = big_model.DEM;
      Model<char> *flow_directions = big_model.flow_directions[0];

      for (int i = 0; i < 9; i++)
        if (big_model.neighbors[i].lat == (int)(FLOOR(reservoirs[i].latitude) - EPS) &&
            big_model.neighbors[i].lon == (int)(FLOOR(reservoirs[i].longitude) + EPS))
          flow_directions = big_model.flow_directions[i];

      reservoirs[i].pour_point = convert_coordinates(
          GeographicCoordinate_init(reservoirs[i].latitude, reservoirs[i].longitude),
          flow_directions->get_origin());
      ArrayCoordinate offset = convert_coordinates(
          convert_coordinates(ArrayCoordinate_init(0, 0, flow_directions->get_origin())),
          DEM->get_origin());
      ArrayCoordinate reservoir_big_ac =
          convert_coordinates(convert_coordinates(reservoirs[i].pour_point), DEM->get_origin());

      reservoir->volume = 0;
      reservoir->area = 0;

      queue<ArrayCoordinate> q;
      q.push(reservoirs[i].pour_point);
      while (!q.empty()) {
        ArrayCoordinate p = q.front();
        q.pop();
        full_cur_model->set(p.row + offset.row, p.col + offset.col, 1);

        ArrayCoordinate full_big_ac = {p.row + offset.row, p.col + offset.col, DEM->get_origin()};

        reservoir->volume +=
            (reservoir->dam_height - (DEM->get(full_big_ac.row, full_big_ac.col) -
                                      DEM->get(reservoir_big_ac.row, reservoir_big_ac.col))) *
            find_area(full_big_ac) / 100;
        reservoir->area += find_area(p);
        update_reservoir_boundary(reservoir->shape_bound, p);

        for (uint d = 0; d < directions.size(); d++) {
          ArrayCoordinate neighbor = {p.row + directions[d].row, p.col + directions[d].col,
                                      p.origin};
          if (flow_directions->check_within(neighbor.row, neighbor.col) &&
              flows_to(neighbor, p, flow_directions) &&
              ((DEM->get(neighbor.row + offset.row, neighbor.col + offset.col) -
                DEM->get(reservoir_big_ac.row, reservoir_big_ac.col)) < reservoir->dam_height)) {
            q.push(neighbor);
          }
        }
      }

      vector<ArrayCoordinate> reservoir_polygon =
          convert_to_polygon(full_cur_model, offset, reservoirs[i].pour_point, 1);

      // DAM
      reservoir->dam_volume = 0;
      reservoir->dam_length = 0;
      bool is_turkeys_nest = true;
      vector<vector<ArrayCoordinate>> dam_polygon;
      bool last = false;
      vector<bool> polygon_bool;
      for (uint i = 0; i < reservoir_polygon.size(); i++) {
        ArrayCoordinate point1 = reservoir_polygon[i];
        ArrayCoordinate point2 = reservoir_polygon[(i + 1) % reservoir_polygon.size()];
        if (is_dam_wall(point1, point2, DEM, offset,
                        reservoir->elevation + reservoir->dam_height)) {
          polygon_bool.push_back(true);

          ArrayCoordinate *adjacent = get_adjacent_cells(point1, point2);
          double average_height =
              (DEM->get(adjacent[0].row + offset.row, adjacent[0].col + offset.col) +
               DEM->get(adjacent[1].row + offset.row, adjacent[1].col + offset.col)) /
              2.0;
          if (full_cur_model->get(adjacent[0].row + offset.row, adjacent[0].col + offset.col) == 0)
            full_cur_model->set(adjacent[0].row + offset.row, adjacent[0].col + offset.col, 2);
          else
            full_cur_model->set(adjacent[1].row + offset.row, adjacent[1].col + offset.col, 2);
          if (!last) {
            vector<ArrayCoordinate> temp;
            temp.push_back(point1);
            dam_polygon.push_back(temp);
          }
          double height = reservoir->elevation + reservoir->dam_height - average_height;
          double length = find_distance(point1, point2) * 1000;
          reservoir->dam_volume += convert_to_dam_volume(height, length);
          reservoir->dam_length += length;
          dam_polygon[dam_polygon.size() - 1].push_back(point2);
          last = true;
        } else {
          is_turkeys_nest = false;
          polygon_bool.push_back(false);
          last = false;
        }
      }
      if (polygon_bool[0] && polygon_bool[dam_polygon.size() - 1] && !is_turkeys_nest) {
        for (uint i = 0; i < dam_polygon[dam_polygon.size() - 1].size(); i++) {
          dam_polygon[0].push_back(dam_polygon[dam_polygon.size() - 1][i]);
        }
        dam_polygon.pop_back();
      }

      if (is_turkeys_nest) {
        ArrayCoordinate *adjacent = get_adjacent_cells(dam_polygon[0][0], dam_polygon[0][1]);
        ArrayCoordinate to_check = adjacent[1];
        if (full_cur_model->get(adjacent[0].row + offset.row, adjacent[0].col + offset.col) == 2)
          to_check = adjacent[0];
        vector<GeographicCoordinate> polygon = compress_poly(
            corner_cut_poly(convert_poly(convert_to_polygon(full_cur_model, offset, to_check, 1))));
        string polygon_string =
            str(polygon, reservoir->elevation + reservoir->dam_height + freeboard);
        coordinates->dam.push_back(polygon_string);
      } else {
        for (uint i = 0; i < dam_polygon.size(); i++) {
          ArrayCoordinate *adjacent = get_adjacent_cells(dam_polygon[i][0], dam_polygon[i][1]);
          ArrayCoordinate to_check = adjacent[1];
          if (full_cur_model->get(adjacent[0].row + offset.row, adjacent[0].col + offset.col) == 2)
            to_check = adjacent[0];
          vector<GeographicCoordinate> polygon = compress_poly(corner_cut_poly(
              convert_poly(convert_to_polygon(full_cur_model, offset, to_check, 2))));
          string polygon_string =
              str(polygon, reservoir->elevation + reservoir->dam_height + freeboard);
          coordinates->dam.push_back(polygon_string);
        }
      }

      reservoir->volume += (reservoir->dam_volume) / 2;
      reservoir->water_rock = reservoir->volume / reservoir->dam_volume;
      reservoir->average_water_depth = reservoir->volume / reservoir->area;

      string polygon_string = str(compress_poly(corner_cut_poly(convert_poly(reservoir_polygon))),
                                  reservoir->elevation + reservoir->dam_height);
      coordinates->reservoir = polygon_string;

      coordinates->is_turkeys_nest = is_turkeys_nest;

      kml_file << output_kml(reservoir, *coordinates);
      kml_file.close();
    }
  }

  printf("Reservoir constructor finished for %s. Runtime: %.2f sec\n",
         convert_string(str(square_coordinate)), 1.0e-6 * (walltime_usec() - t_usec));
}

// void update_kml_holder(KML_Holder* kml_holder, Pair* pair, Pair_KML* pair_kml){
// 	kml_holder->points.push_back(get_point_kml(pair, pair_kml->point));
// 	kml_holder->lines.push_back(get_line_kml(pair, pair_kml->line));
// 	kml_holder->uppers.push_back(get_reservoir_kml(&pair->upper, upper_colour,
// pair_kml->upper)); 	kml_holder->upper_dams.push_back(get_dam_kml(&pair->upper, pair_kml->upper));
// 	kml_holder->lowers.push_back(get_reservoir_kml(&pair->lower, lower_colour,
// pair_kml->lower)); 	kml_holder->lower_dams.push_back(get_dam_kml(&pair->lower, pair_kml->lower));
// }
