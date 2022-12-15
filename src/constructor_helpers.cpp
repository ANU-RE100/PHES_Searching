#include "constructor_helpers.hpp"
#include "coordinates.h"
#include "fsolve.hpp"

/*
 * Returns sorted vector containing the longitude of all polgon boundary interections at certain
 * latitude. Assumes last coordinate of polygon is the same as the first
 */
vector<double> find_polygon_intersections(double lat, vector<GeographicCoordinate> &polygon) {
  vector<double> to_return;
  for (uint i = 0; i < polygon.size() - 1; i++) {
    GeographicCoordinate line[2] = {polygon[i], polygon[(i + 1)]};
    if ((line[0].lat < lat && line[1].lat >= lat) || (line[0].lat >= lat && line[1].lat < lat)) {
      to_return.push_back(line[0].lon + (lat - line[0].lat) / (line[1].lat - line[0].lat) *
                                            (line[1].lon - line[0].lon));
    }
  }
  sort(to_return.begin(), to_return.end());
  return to_return;
}

/*
 * Check if a geographic point is inside a polygon
 */
bool check_within(GeographicCoordinate point, vector<vector<GeographicCoordinate>> polygons){
    for(vector<GeographicCoordinate> polygon:polygons){
        vector<double> polygon_intersections = find_polygon_intersections(point.lat, polygon);
        for(uint j = 0; j<polygon_intersections.size()/2;j++)
            if(polygon_intersections[2*j]<=point.lon && point.lon<=polygon_intersections[2*j+1])
                return true;
    }
    return false;
}

/*
 * Read countries from custom country file in /input/countries/countries.txt
 */
vector<vector<vector<GeographicCoordinate>>> read_countries(string filename, vector<string>& country_names){
    char *shp_filename = new char[filename.length() + 1];
    strcpy(shp_filename, filename.c_str());
    if(!file_exists(shp_filename)){
        cout << "No file: "+filename << endl;
        throw(1);
    }
    vector<vector<vector<GeographicCoordinate>>> relevant_polygons;

    ifstream in(filename);
    string line;
    while(getline(in, line)){
        if(line[0]=='@'){
            line.erase(remove(line.begin(), line.end(), '@'), line.end());
            line.erase(line.find_last_not_of(" \n\r\t")+1);
            country_names.push_back(line);
            vector<vector<GeographicCoordinate>> temp;
            relevant_polygons.push_back(temp);
        }else{
            vector<string> coordinates = read_from_csv_file(line);
            vector<GeographicCoordinate> temp;
            relevant_polygons[relevant_polygons.size()-1].push_back(temp);
            for(uint i = 0; i<coordinates.size()/2; i++){
                GeographicCoordinate c = {stod(coordinates[2*i+1]), stod(coordinates[2*i])};
                relevant_polygons[relevant_polygons.size()-1][relevant_polygons[relevant_polygons.size()-1].size()-1].push_back(c);
            }
        }
    }
    return relevant_polygons;
}

/*
 * Returns a tuple with the two cells adjacent to an edge defined by two points
 */
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

/*
 * Determine if two points define an edge of a reservoir given its raster model
 */
bool is_edge(ArrayCoordinate point1, ArrayCoordinate point2, Model<char>* model, ArrayCoordinate offset, int threshold){

    if (point1.row+offset.row<0 || point1.col+offset.col<0 || point1.row+offset.row>model->nrows() || point1.col+offset.col>model->ncols())
        return false;
    if (point2.row+offset.row<0 || point2.col+offset.col<0 || point2.row+offset.row>model->nrows() || point2.col+offset.col>model->ncols())
        return false;

    ArrayCoordinate* to_check = get_adjacent_cells(point1, point2);

    if((!model->check_within(to_check[0].row+offset.row,to_check[0].col+offset.col) && model->get(to_check[1].row+offset.row,to_check[1].col+offset.col)>=threshold)||
        (!model->check_within(to_check[1].row+offset.row,to_check[1].col+offset.col) && model->get(to_check[0].row+offset.row,to_check[0].col+offset.col)>=threshold))
        return true;

    if (model->get(to_check[0].row+offset.row,to_check[0].col+offset.col)<threshold && model->get(to_check[1].row+offset.row,to_check[1].col+offset.col)>=threshold)
        return true;
    if (model->get(to_check[1].row+offset.row,to_check[1].col+offset.col)<threshold && model->get(to_check[0].row+offset.row,to_check[0].col+offset.col)>=threshold)
        return true;

    return false;
}

/*
 * Determines if an edge of a reservoir between two points requires a dam wall
 */
bool is_dam_wall(ArrayCoordinate point1, ArrayCoordinate point2, Model<short>* DEM, ArrayCoordinate offset, double wall_elevation){
    if (point1.row<0 || point1.col<0 || point1.row>DEM->nrows() || point1.col>DEM->ncols())
        return false;
    if (point2.row<0 || point2.col<0 || point2.row>DEM->nrows() || point2.col>DEM->ncols())
        return false;

    double average_row = (point1.row+point2.row)/2.0;
    double average_col = (point1.col+point2.col)/2.0;

    ArrayCoordinate* to_check = get_adjacent_cells(point1, point2);

    if(average_row-(int)average_row>0.9 || average_row-(int)average_row<0.1){
        to_check[0] = ArrayCoordinate_init((int)(average_row+EPS),(int)(average_col-0.5+EPS), point1.origin);
    	to_check[1] = ArrayCoordinate_init((int)(average_row-1+EPS),(int)(average_col-0.5+EPS), point1.origin);
    }else{
    	to_check[0] = ArrayCoordinate_init((int)(average_row-0.5+EPS),(int)(average_col+EPS), point1.origin);
    	to_check[1] = ArrayCoordinate_init((int)(average_row-0.5+EPS),(int)(average_col-1+EPS), point1.origin);
    }

    if(!DEM->check_within(to_check[0].row, to_check[0].col) || !DEM->check_within(to_check[1].row, to_check[1].col))
        return false;

    if(DEM->get(to_check[0].row+offset.row, to_check[0].col+offset.col)<wall_elevation and DEM->get(to_check[1].row+offset.row, to_check[1].col+offset.col)<wall_elevation)
        return true;
    return false;
}

int dir_def[4][2] = {{-1,0},{1,0},{0,-1},{0,1}};
int dir_to_do[4][3] = {{2,0,3},{3,1,2},{1,2,0},{0,3,1}};
int testsa[] = {1, 2, 0, 3};

/*
 * Converts a raster model to a polygon given a raster model and a point on the interior edge of the polygon
 */
vector<ArrayCoordinate> convert_to_polygon(Model<char>* model, ArrayCoordinate offset, ArrayCoordinate pour_point, int threshold){

    vector<ArrayCoordinate> to_return;

    ArrayCoordinate test_coordinates[] = {
        ArrayCoordinate_init(pour_point.row+1, pour_point.col+1, pour_point.origin),
        ArrayCoordinate_init(pour_point.row+1, pour_point.col, pour_point.origin),
        ArrayCoordinate_init(pour_point.row, pour_point.col, pour_point.origin),
        ArrayCoordinate_init(pour_point.row, pour_point.col+1, pour_point.origin)};

    bool succesful_path = false;
    for(int i = 0; i<4; i++){
        vector<ArrayCoordinate> temp_to_return;
        int last_dir = testsa[i];
        ArrayCoordinate initial = test_coordinates[i];
        ArrayCoordinate last = initial;
        bool found_path;

        while(true){
            found_path = false;
            for(int id = 0; id<3; id++){
                int d = dir_to_do[last_dir][id];
                ArrayCoordinate next = ArrayCoordinate_init(last.row+dir_def[d][0],last.col+dir_def[d][1], pour_point.origin);
                if(is_edge(last, next, model, offset, threshold)){
                    temp_to_return.push_back(next);
                    last = next;
                    last_dir = d;
                    found_path = true;
                    break;
                }
            }
            if(!found_path)
                break;
            if(last.row==initial.row && last.col == initial.col){
                succesful_path = true;
                break;
            }
        }
        if(temp_to_return.size()>to_return.size()){
            to_return.clear();
            for(uint j = 0; j<temp_to_return.size(); j++){
                to_return.push_back(temp_to_return[j]);
            }
        }
    }
    if(!succesful_path){
      search_config.logger.error("Could not find a succesful path around the polygon.");
      throw 1;
    }
    return to_return;
}

/*
 * Convert polygon of grid vertices to polygon of geographic coordinates
 */
vector<GeographicCoordinate> convert_poly(vector<ArrayCoordinate> polygon){
    vector<GeographicCoordinate> to_return;
    for(uint i = 0; i<polygon.size(); i++){
    	to_return.push_back(convert_coordinates(polygon[i], 0));
    }
    return to_return;
}

/*
 * "Smooth" a polygon by cutting the corners (take midpoint of each side)
 */
vector<GeographicCoordinate> corner_cut_poly(vector<GeographicCoordinate> polygon){
    vector<GeographicCoordinate> to_return;
    for(uint i = 0; i<polygon.size(); i++){
    	to_return.push_back(GeographicCoordinate_init((polygon[i].lat*2+polygon[(i+1)%polygon.size()].lat*2)/4.0, (polygon[i].lon*2+polygon[(i+1)%polygon.size()].lon*2)/4.0));
    }
    to_return.push_back(to_return[0]);
    return to_return;
}

/*
 * Compresses a polygon by removing collinear points
 */
vector<GeographicCoordinate> compress_poly(vector<GeographicCoordinate> polygon){
    vector<GeographicCoordinate> to_return;
    to_return.push_back(polygon[0]);
    for(uint i = 1; i<polygon.size()-1; i++){
        if ((polygon[i+1].lon-polygon[i].lon)/(polygon[i+1].lat-polygon[i].lat) < (polygon[i].lon-polygon[i-1].lon)/(polygon[i].lat-polygon[i-1].lat)-EPS
        	|| (polygon[i+1].lon-polygon[i].lon)/(polygon[i+1].lat-polygon[i].lat) > (polygon[i].lon-polygon[i-1].lon)/(polygon[i].lat-polygon[i-1].lat)+EPS ){
            to_return.push_back(polygon[i]);
        }
    }
    to_return.push_back(polygon[polygon.size()-1]);
    return to_return;
}

/*
 * Convert geographic polygon to KML coordinate string
 */
string str(vector<GeographicCoordinate> polygon, double elevation) {
  string to_return = "";
  for (uint i = 0; i < polygon.size(); i++) {
    to_return += dtos(polygon[i].lon, 5) + "," + dtos(polygon[i].lat, 5) + "," +
                 dtos(elevation, 1) + " ";
  }
  return to_return;
}

/*
* System of equations solved to find the optimal dam wall height of a turkey nest site
*/
double TN_nonlinear_system(double dam_height[],
                           vector<ArrayCoordinateWithHeight> reservoir_points,
                           vector<ArrayCoordinateWithHeight> dam_points, double req_volume, double reservoir_area) {
  vector<double> reservoir_height_diffs;
  vector<double> dam_height_sqdiffs;
  double volume_original = 0;
  double volume_earthwork = 0;
  double reservoir_diffs_sum = 0;
  double dam_diffs2_sum = 0;
  double dam_length_at_height = 0;

  // Calculate the length of the dam wall for the specified dam wall height
  dam_length_at_height = turkey_dam_length(dam_points, dam_height[0]);

  // Set up the vectors of elevation differences
  for (uint point_index = 0; point_index < reservoir_points.size(); point_index++) {
    reservoir_height_diffs.push_back(max(0.0, dam_height[0] - reservoir_points[point_index].h));
  }
  for (uint point_index = 0; point_index < dam_points.size(); point_index++) {
    dam_height_sqdiffs.push_back(max(0.0, dam_height[0] - dam_points[point_index].h + freeboard) * (cwidth + dambatter * (dam_height[0] - dam_points[point_index].h + freeboard)));
  }

  // Volume of original storage capacity (without earthwork)
  reservoir_diffs_sum = accumulate(reservoir_height_diffs.begin(), reservoir_height_diffs.end(), 0.0);
  volume_original = 10000 * reservoir_area * reservoir_diffs_sum / reservoir_height_diffs.size();

  // Volume of dam i.e. required earthwork
  dam_diffs2_sum = accumulate(dam_height_sqdiffs.begin(), dam_height_sqdiffs.end(), 0.0);
  volume_earthwork = dam_length_at_height * dam_diffs2_sum / dam_height_sqdiffs.size();

  return volume_earthwork / 2 + volume_original - (req_volume * 1000000); // Convert GL to m3
}

/*
 * Accurately model a single reservoir, determining optimal dam wall height for given volume.
 *
 * Pass negative reservoir volume to model single dam wall height
 */
bool model_reservoir(Reservoir *reservoir, Reservoir_KML_Coordinates *coordinates,
                     Model<bool> *seen, bool *non_overlap, vector<ArrayCoordinate> *used_points,
                     BigModel big_model, Model<char> *full_cur_model,
                     vector<vector<vector<GeographicCoordinate>>> &countries,
                     vector<string> &country_names) {
  
  ///////////////////////////////////////////
  /* int count1 = 0;
  for (int row = 0; row < full_cur_model->nrows(); row++) {
    for (int col = 0; col < full_cur_model->ncols(); col++) {
      if (full_cur_model->get(row, col) == 1)
        count1++;
    }
  }  
  printf("Count1: %d\n", count1); */
  ///////////////////////////////////////////

  Model<short> *DEM = big_model.DEM;
  Model<char> *flow_directions = big_model.flow_directions[0];

  for (int i = 0; i < 9; i++)
    if (big_model.neighbors[i].lat == convert_to_int(FLOOR(reservoir->latitude + EPS)) &&
        big_model.neighbors[i].lon == convert_to_int(FLOOR(reservoir->longitude + EPS)))
      flow_directions = big_model.flow_directions[i];

  ArrayCoordinate offset = convert_coordinates(
      convert_coordinates(ArrayCoordinate_init(0, 0, flow_directions->get_origin())),
      DEM->get_origin());

  ArrayCoordinate reservoir_big_ac =
      convert_coordinates(convert_coordinates(reservoir->pour_point), DEM->get_origin());

  double req_volume = reservoir->volume;
  reservoir->volume = 0;
  reservoir->area = 0;
  vector<ArrayCoordinate> temp_used_points;

  // GREENFIELD RESERVOIR
  if (!reservoir->turkey) {
    char last_dir = 'd';
    while ((req_volume > 0 &&
            (reservoir->volume * (1 + 0.5 / reservoir->water_rock) <
                (1 - volume_accuracy) * req_volume ||
            reservoir->volume * (1 + 0.5 / reservoir->water_rock) >
                (1 + volume_accuracy) * req_volume)) ||
          reservoir->volume == 0) {
      temp_used_points.clear();
      reservoir->volume = 0;
      reservoir->area = 0;

      queue<ArrayCoordinate> q;
      q.push(reservoir->pour_point);
      while (!q.empty()) {
        ArrayCoordinate p = q.front();
        q.pop();

        if (full_cur_model != NULL)
          full_cur_model->set(p.row + offset.row, p.col + offset.col, 1);

        ArrayCoordinate full_big_ac = {p.row + offset.row, p.col + offset.col, DEM->get_origin()};

        temp_used_points.push_back(full_big_ac);
        if (seen != NULL && seen->get(full_big_ac.row, full_big_ac.col)){
          if (non_overlap != NULL){
            *non_overlap = false;
          } else {
            return false;
          }
        }
        if (DEM->get(full_big_ac.row, full_big_ac.col) < -2000)
          return false;

        reservoir->volume +=
            (reservoir->dam_height -
            (DEM->get(full_big_ac.row, full_big_ac.col) -
              DEM->get(reservoir_big_ac.row, reservoir_big_ac.col))) *
            find_area(full_big_ac) / 100;
        reservoir->area += find_area(p);
        update_reservoir_boundary(reservoir->shape_bound, p);

        for (uint d = 0; d < directions.size(); d++) {
          ArrayCoordinate neighbor = {p.row + directions[d].row,
                                      p.col + directions[d].col, p.origin};
          if (flow_directions->check_within(neighbor.row, neighbor.col) &&
              flow_directions->flows_to(neighbor, p) &&
              ((DEM->get(neighbor.row + offset.row, neighbor.col + offset.col) -
                DEM->get(reservoir_big_ac.row, reservoir_big_ac.col)) <
              reservoir->dam_height)) {
            q.push(neighbor);
          }
        }
      }

      if (req_volume >0 && reservoir->volume * (1 + 0.5 / reservoir->water_rock) <
          (1 - volume_accuracy) * req_volume) {
        reservoir->dam_height += dam_wall_height_resolution;
        if (reservoir->dam_height > reservoir->max_dam_height)
          return false;
        last_dir = 'u';
      }

      if (req_volume >0 && reservoir->volume * (1 + 0.5 / reservoir->water_rock) >
          (1 + volume_accuracy) * req_volume) {
        if (last_dir == 'u')
          return false;
        reservoir->dam_height -= dam_wall_height_resolution;
        last_dir = 'd';
      }
    } 
  } else { // TURKEY NEST RESERVOIR
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
    queue<ArrayCoordinate> q;
    queue<ArrayCoordinate> empty;
    ArrayCoordinate c = reservoir->pour_point;
    ArrayCoordinate neighbor;
    ArrayCoordinate dam_neighbor;
    ArrayCoordinateWithHeight dam_neighbor_with_height;
    vector<ArrayCoordinateWithHeight> dam_neighbors;
    ArrayCoordinateWithHeight c_with_height = ArrayCoordinateWithHeight_init(c.row,c.col,DEM->get(c.row,c.col));
    ArrayCoordinateWithHeight neighbor_with_height;
    vector<ArrayCoordinateWithHeight> dam_points;
    vector<double> dam_elevations;
    vector<ArrayCoordinateWithHeight> reservoir_points; 
    Reservoir* test_reservoir;
    bool continue_test = false;
    bool growing_adjustment = false;

    reservoir->area = 0;
    reservoir->volume = 0;
    reservoir->dam_volume = 0;
    reservoir->water_rock = 0;
    reservoir->area += find_area(c);
    reservoir_points.push_back(c_with_height); 

    q.push(c);
    reservoir_points.push_back(c_with_height);
    temp_used_points.clear();
    temp_used_points.push_back(c);
    
    if (full_cur_model != NULL)
                full_cur_model->set(c.row, c.col, 1);

    // Model the rough reservoir that has the required volume
    while (!q.empty() && reservoir->volume < req_volume) {
      c = q.front();
      q.pop();

      for (uint d = 0; d < directions.size(); d++) {
        if (directions[d].row * directions[d].col == 0) {
          neighbor = ArrayCoordinate_init(c.row + directions[d].row, c.col + directions[d].col, DEM->get_origin());
          neighbor_with_height = ArrayCoordinateWithHeight_init(neighbor.row,neighbor.col,DEM->get(neighbor.row,neighbor.col));
          
          // && filter->get(neighbor.row, neighbor.col) == false
          if (DEM->check_within(c.row + directions[d].row, c.col + directions[d].col) && !seen->get(neighbor.row, neighbor.col) && DEM->get_slope(neighbor.row, neighbor.col) <= TN_elevation_tolerance) {
            temp_used_points.push_back(neighbor);

            dam_points.push_back(neighbor_with_height);
            reservoir_points.push_back(neighbor_with_height);
 
            reservoir->area += find_area(neighbor);
            
            if (full_cur_model != NULL)
                full_cur_model->set(neighbor.row, neighbor.col, 1);

            q.push(neighbor);
          }  
        }              
      }
    }

    reservoir = update_TN_volumes(dam_points, reservoir_points, reservoir);
    update_reservoir_boundary(reservoir->shape_bound, dam_points);

    // Modify the boundaries of the reservoir to increase the water-to-rock ratio
    swap(q, empty);

    for (uint point_index = 0; point_index < dam_points.size(); point_index++)
      q.push(ArrayCoordinate_init(dam_points[point_index].row, dam_points[point_index].col, get_origin(search_config.grid_square, border)));

    while(!q.empty()) {
      c = q.front();
      q.pop();
      test_reservoir = reservoir;

      // Remove the point from the dam points
      for (uint point_index = 0; point_index < dam_points.size(); point_index++) {
        if (dam_points[point_index] == c_with_height) {
          dam_points.erase(dam_points.begin() + point_index);
          break;
        }
      }

      for (uint d = 0; d < directions.size(); d++) {
        if (directions[d].row * directions[d].col == 0) {
          // Determine the test point
          neighbor = ArrayCoordinate_init(c.row + directions[d].row, c.col + directions[d].col, DEM->get_origin());
          neighbor_with_height = ArrayCoordinateWithHeight_init(neighbor.row,neighbor.col,DEM->get(neighbor.row,neighbor.col));

          // Make sure the new point maintains dam integrity
          dam_neighbors.clear();
          for (uint d = 0; d < directions.size(); d++) {
            dam_neighbor = ArrayCoordinate_init(c.row + directions[d].row, c.col + directions[d].col, DEM->get_origin());
            dam_neighbor_with_height = ArrayCoordinateWithHeight_init(dam_neighbor.row, dam_neighbor.col, DEM->get(dam_neighbor.row, dam_neighbor.col));

            // Find all neighboring points that are part of the dam wall
            for (uint point_index = 0; point_index < dam_points.size(); point_index++) {
              if (dam_points[point_index] == dam_neighbor_with_height) {
                dam_neighbors.push_back(dam_neighbor_with_height);                  
              }
            }
          }

          // Check that at least 2 of the dam wall points maintain integrity of the wall with the neighbor point
          for (uint point_index1 = 0; point_index1 < dam_neighbors.size(); point_index1++) {
            for (uint point_index2 = 0; point_index2 < dam_neighbors.size(); point_index2++) {
              if (dam_neighbors[point_index1].row == dam_neighbors[point_index2].row + 2 || dam_neighbors[point_index1].row == dam_neighbors[point_index2].row - 2 ||dam_neighbors[point_index1].col == dam_neighbors[point_index2].col + 2 ||dam_neighbors[point_index1].col == dam_neighbors[point_index2].col - 2 ) {
                continue_test = true;
                break;
              }
            }
            if (continue_test)
              break;
          }

          // Add the new point to the dam vector, either add or remove to the vector
          dam_points.push_back(neighbor_with_height);

          for (uint point_index1 = 0; point_index1 < reservoir_points.size(); point_index1++) {
            if(reservoir_points[point_index1] == neighbor_with_height) {
              for (uint point_index2 = 0; point_index2 < reservoir_points.size(); point_index2++) {
                if(reservoir_points[point_index2] == c_with_height) {
                  reservoir_points.erase(reservoir_points.begin() + point_index2);
                  test_reservoir->area -= find_area(c);
                  break;
                }
              }

              for (uint point_index2 = 0; point_index2 < temp_used_points.size(); point_index2++) {
                if(temp_used_points[point_index2] == c) {
                  temp_used_points.erase(temp_used_points.begin() + point_index2);
                  break;
                }
              }

              growing_adjustment = false;
              break;
            } else {
              continue;
            }

            reservoir_points.push_back(neighbor_with_height);
            test_reservoir->area += find_area(neighbor);
            temp_used_points.push_back(neighbor);
            growing_adjustment = true;                 
          }
          
          // Update the dam wall height
          // Set up the inputs for the floating point solver to find the optimal dam height
          dam_elevations.clear();
          for (uint point_index = 0; point_index < dam_points.size(); point_index++) {
            dam_elevations.push_back(dam_points[point_index].h);
          }

          double* fx;
          int info;
          int lwa;
          int n = 1;
          double tol = 0.00001;
          double* wa;
          double optimal_elevation;
          double* dam_height;

          lwa = (n * (3 * n + 13)) / 2;
          wa = new double[lwa];
          dam_height = new double[n];
          fx = new double[n];

          vector<double>::iterator max_dam_elevation = max_element(dam_elevations.begin(), dam_elevations.end());
          vector<double>::iterator min_dam_elevation = min_element(dam_elevations.begin(), dam_elevations.end());

          dam_height[0] = *max_dam_elevation; // Initial guess at maximum dam ground elevation

          // Declare a std::function wrapper for the lambda wrapper for the TN_nonlinear_system function
          function<void(int, double *, double *)> nonlinear_function{
              [&](int n1, double x1[], double fx1[]) {
                fx1[0] =
                    TN_nonlinear_system(x1, reservoir_points, dam_points, req_volume, test_reservoir->area);
                return;
              }};

          nonlinear_function(n, dam_height, fx);

          // Run the floating point solver to find the roots of the non-linear equations
          info = fsolve(nonlinear_function, n, dam_height, fx, tol, wa, lwa);

          optimal_elevation = dam_height[0];
          test_reservoir->dam_height = optimal_elevation - *min_dam_elevation;

          // Calculate the water-rock ratio of the test reservoir
          test_reservoir = update_TN_volumes(dam_points, reservoir_points, test_reservoir);

          // Check if the new reservoir has a higher water/rock ratio. If not, reject the change
          if (test_reservoir->water_rock > reservoir->water_rock) {
            reservoir = test_reservoir;
            q.push(neighbor);

            if (!growing_adjustment && full_cur_model != NULL) {
              full_cur_model->set(c.row, c.col, 0);
            } else if (growing_adjustment && full_cur_model != NULL) {
              full_cur_model->set(neighbor.row, neighbor.col, 1);
            }

          } else {
            // Remove the test point from the dam vectors
            dam_points.pop_back();

            // Add the original point back to the dam
            dam_points.push_back(neighbor_with_height);

            // Add the original point or remove the test point to the reservoir, depending upon if the test reservoir grew or shrunk
            if (growing_adjustment) {
              reservoir_points.pop_back();
              temp_used_points.pop_back();
            } else {
              reservoir_points.push_back(c_with_height);
              temp_used_points.push_back(c);
            }
          }         
        }
      }
    }

    update_reservoir_boundary(reservoir->shape_bound, dam_points);

    printf("Point: %d %d %d %.2f\n", int(reservoir_points.size()) + int(dam_points.size()), int(reservoir_points.size()), int(dam_points.size()), reservoir->volume);
  }

  vector<ArrayCoordinate> reservoir_polygon;
  reservoir_polygon = convert_to_polygon(full_cur_model, offset, reservoir->pour_point, 1);

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

  // full_cur_model->write("out1.tif", GDT_Byte);
  // DAM WALL
  reservoir->dam_volume = 0;
  reservoir->dam_length = 0;
  bool is_turkeys_nest = true;
  vector<vector<ArrayCoordinate>> dam_polygon;
  bool last = false;
  vector<bool> polygon_bool;
  for (uint i = 0; i < reservoir_polygon.size(); i++) {
    ArrayCoordinate point1 = reservoir_polygon[i];
    ArrayCoordinate point2 =
        reservoir_polygon[(i + 1) % reservoir_polygon.size()];
    if (is_dam_wall(point1, point2, DEM, offset,
                    reservoir->elevation + reservoir->dam_height)) {
      polygon_bool.push_back(true);

      ArrayCoordinate *adjacent = get_adjacent_cells(point1, point2);
      double average_height = (DEM->get(adjacent[0].row + offset.row,
                                        adjacent[0].col + offset.col) +
                               DEM->get(adjacent[1].row + offset.row,
                                        adjacent[1].col + offset.col)) /
                              2.0;
      if (full_cur_model->get(adjacent[0].row + offset.row, adjacent[0].col + offset.col) != 1) {
        full_cur_model->set(adjacent[0].row + offset.row, adjacent[0].col + offset.col, 2);

      } else {
        full_cur_model->set(adjacent[1].row + offset.row, adjacent[1].col + offset.col, 2);
      }
      if (!last) {
        vector<ArrayCoordinate> temp;
        temp.push_back(point1);
        dam_polygon.push_back(temp);
      }
      double height =
          reservoir->elevation + reservoir->dam_height - average_height;
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

  //full_cur_model->write("out2.tif", GDT_Byte);
  if (polygon_bool[0] && polygon_bool[dam_polygon.size() - 1] &&
      !is_turkeys_nest && dam_polygon.size() > 1) {
    for (uint i = 0; i < dam_polygon[dam_polygon.size() - 1].size(); i++) {
      dam_polygon[0].push_back(dam_polygon[dam_polygon.size() - 1][i]);
    }
    dam_polygon.pop_back();
  }

  for (uint i = 0; i < dam_polygon.size(); i++) {
    ArrayCoordinate *adjacent = get_adjacent_cells(dam_polygon[i][0], dam_polygon[i][1]);
    ArrayCoordinate to_check = adjacent[1];
    if (full_cur_model->get(adjacent[0].row + offset.row, adjacent[0].col + offset.col) == 2)
      to_check = adjacent[0];
    vector<GeographicCoordinate> polygon;
    try {
      polygon = compress_poly(
          corner_cut_poly(convert_poly(convert_to_polygon(full_cur_model, offset, to_check, 2))));
    } catch (int e) {
      return false;
    }
    string polygon_string =
        str(polygon, reservoir->elevation + reservoir->dam_height + freeboard);
    coordinates->dam.push_back(polygon_string);
  }

  reservoir->volume += (reservoir->dam_volume) / 2;
  reservoir->water_rock = reservoir->volume / reservoir->dam_volume;
  reservoir->average_water_depth = reservoir->volume / reservoir->area;

  string polygon_string =
      str(compress_poly(corner_cut_poly(convert_poly(reservoir_polygon))),
          reservoir->elevation + reservoir->dam_height);
  coordinates->reservoir = polygon_string;

  coordinates->is_turkeys_nest = is_turkeys_nest;

  queue<ArrayCoordinate> q;
  q.push(reservoir->pour_point);
  while (!q.empty()) {
    ArrayCoordinate p = q.front();
    q.pop();
    full_cur_model->set(p.row + offset.row, p.col + offset.col, 0);
    for (uint d = 0; d < directions.size(); d++) {
      ArrayCoordinate neighbor = {p.row + directions[d].row,
                                  p.col + directions[d].col, p.origin};
      if (full_cur_model->get(neighbor.row + offset.row,
                              neighbor.col + offset.col) != 0) {
        full_cur_model->set(neighbor.row + offset.row,
                            neighbor.col + offset.col, 0);
        q.push(neighbor);
      }
    }
  }
  for (uint i = 0; i < countries.size(); i++) {
    if (check_within(convert_coordinates(reservoir->pour_point),
                     countries[i])) {
      reservoir->country = country_names[i];
      break;
    }
  }

  return true;
}
