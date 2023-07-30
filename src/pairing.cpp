#include "coordinates.h"
#include "phes_base.h"
#include "reservoir.h"
#include "search_config.hpp"

vector<ExistingPit> pit_details;
ExistingPit single_pit_details;

vector<int> pairs;

vector<GeographicCoordinate> find_points_to_test(RoughReservoir* &reservoir,
                                                 double &wall_height, ArrayCoordinate &pour_point) {
  vector<GeographicCoordinate> bound;
  if (RoughGreenfieldReservoir *gr = dynamic_cast<RoughGreenfieldReservoir *>(reservoir)) {
    array<ArrayCoordinate, directions.size()> one_point = {pour_point, pour_point, pour_point,
                                                           pour_point, pour_point, pour_point,
                                                           pour_point, pour_point};
    int i = 0;
    while (dam_wall_heights[i] < wall_height) {
      i += 1;
    }
    int lower_wall_height = (i) ? dam_wall_heights[i - 1] : 0;
    array<ArrayCoordinate, directions.size()> lower_shape =
        (i) ? gr->shape_bound[i - 1] : one_point;
    double inv_wall_height_interval = 0.1;
    for (uint j = 0; j < directions.size(); j++) {
      GeographicCoordinate point1 = convert_coordinates(lower_shape[j]);
      GeographicCoordinate point2 = convert_coordinates(gr->shape_bound[i][j]);
      bound.push_back((GeographicCoordinate){
          point1.lat + (point2.lat - point1.lat) * (wall_height - lower_wall_height) *
                           inv_wall_height_interval,
          point1.lon + (point2.lon - point1.lon) * (wall_height - lower_wall_height) *
                           inv_wall_height_interval});
    }
    bound.push_back(convert_coordinates(pour_point));
  } else {
    RoughBfieldReservoir *br = dynamic_cast<RoughBfieldReservoir *>(reservoir);
    for (ArrayCoordinate c : br->shape_bound)
      bound.push_back(convert_coordinates(c));
  }
  return bound;
}

double find_least_distance_sqd(RoughReservoir* upper, RoughReservoir* &lower,
                               double upper_wall_height, double lower_wall_height,
                               ArrayCoordinate* upper_pour_point, ArrayCoordinate* lower_pour_point) {
  double mindist2 = INF;
  vector<GeographicCoordinate> upper_points =
      find_points_to_test(upper, upper_wall_height, *upper_pour_point);
  vector<GeographicCoordinate> lower_points =
      find_points_to_test(lower, lower_wall_height, *lower_pour_point);

  for (uint iu = 0; iu < upper_points.size(); iu++) {
    GeographicCoordinate p1 = upper_points[iu];
    for (uint il = 0; il < lower_points.size(); il++) {
      GeographicCoordinate p2 = lower_points[il];
      if (mindist2 > find_distance_sqd(p1, p2)){
        mindist2 = find_distance_sqd(p1, p2);
        *upper_pour_point = convert_coordinates(p1,upper_pour_point->origin);
        *lower_pour_point = convert_coordinates(p2,lower_pour_point->origin);
      }
    }
  }

  return mindist2;
}

int max_altitude(vector<AltitudeVolumePair> pairs) {
  return pairs[pairs.size() - 1].altitude;
}

vector<int> get_altitudes(ExistingPit &pit) {
  vector<int> to_return;
  for (AltitudeVolumePair pair : pit.volumes)
    to_return.push_back(pair.altitude);
  return to_return;
}

vector<double> get_volumes(ExistingPit &pit) {
  vector<double> to_return;
  for (AltitudeVolumePair pair : pit.volumes)
    to_return.push_back(pair.volume);
  return to_return;
}

vector<double> int_to_double_vector(vector<int> int_vector) {
  vector<double> to_return(int_vector.begin(), int_vector.end());
  return to_return;
}

double pit_volume(ExistingPit &pit, int bottom_elevation, int top_elevation) {
  return linear_interpolate(top_elevation,
                            int_to_double_vector(get_altitudes(pit)),
                            get_volumes(pit)) -
         linear_interpolate(bottom_elevation,
                            int_to_double_vector(get_altitudes(pit)),
                            get_volumes(pit));
}

bool determine_pit_elevation_and_volume(RoughReservoir* &upper,
                                        RoughReservoir* &lower,
                                        double energy_capacity,
                                        ExistingPit &pit_details_single,
                                        double &required_volume, int &head) {
  RoughReservoir* greenfield = upper;
  RoughReservoir* pit = lower;
  if (upper->brownfield) {
    greenfield = lower;
    pit = upper;
  }

  while (pit->elevation < max_altitude(pit_details_single.volumes)) {
    pit->max_dam_height = max_altitude(pit_details_single.volumes) - pit->elevation;
    int pit_depth = 0;
    while (pit_depth < pit->max_dam_height) {
      pit_depth += 1;
      double volume =
          pit_volume(pit_details_single, pit->elevation, pit->elevation + pit_depth);
      double greenfield_wall_height =
          linear_interpolate(volume, greenfield->volumes, dam_wall_heights);
      head = convert_to_int(ABS(((0.5 * (double)greenfield_wall_height +
                        (double)greenfield->elevation) -
                       (0.5 * (double)pit_depth + (double)pit->elevation))));
      if (head < min_head || head > max_head)
        continue;
      double head_ratio =
          (head + 0.5 * (greenfield_wall_height + (double)pit_depth)) /
          (head - 0.5 * (greenfield_wall_height + (double)pit_depth));
      // cout << volume << " " << greenfield_wall_height << " " <<
      // greenfield->elevation << " " << pit_depth << " " << pit->elevation << " "
      // << head << " " << head_ratio << "\n";

      if (head_ratio > (1 + max_head_variability)) {
        break;
      }

      if (volume < find_required_volume(energy_capacity, head)) {
        continue;
      }
      required_volume = volume;
      return true;
    }

    pit->elevation += pit_height_resolution;
  }
  return false;
}

bool determine_pit_elevation_and_volume(RoughReservoir &pit,
                                        RoughReservoir &greenfield,
                                        double energy_capacity,
                                        double &required_volume, int &head) {

  int max_fill_depth = *max_element(pit.fill_depths.begin(), pit.fill_depths.end());
 
  vector<double> p_fill_depth_doubles;
  vector<double> g_fill_depth_doubles;
  for (uint i=0; i< pit.fill_depths.size(); i++){
    p_fill_depth_doubles.push_back((double)pit.fill_depths[i]);
    g_fill_depth_doubles.push_back((double)greenfield.fill_depths[i]);
  }
    

  while (pit.elevation < pit.bottom_elevation + max_fill_depth) {
    pit.max_dam_height = pit.bottom_elevation + max_fill_depth - pit.elevation;
    pit.fill_depth_from_MOL = 0;
    while (pit.fill_depth_from_MOL < pit.max_dam_height) {
      pit.fill_depth_from_MOL += 1;
      double volume = linear_interpolate(pit.fill_depth_from_MOL, p_fill_depth_doubles, pit.volumes);
      double greenfield_wall_height = linear_interpolate(volume, greenfield.volumes, g_fill_depth_doubles);
      
      double max_head_pair = MAX(MAX(greenfield.elevation + greenfield_wall_height - pit.elevation, greenfield.elevation - (pit.elevation + pit.fill_depth_from_MOL)),0);
      double min_head_pair = MAX(MIN(greenfield.elevation + greenfield_wall_height - pit.elevation, greenfield.elevation - (pit.elevation + pit.fill_depth_from_MOL)),0);

      head = convert_to_int(ABS(((0.5 * (double)greenfield_wall_height +
                        (double)greenfield.elevation) -
                       (0.5 * (double)pit.fill_depth_from_MOL + (double)pit.elevation))));
      if (head < min_head || head > max_head)
        continue;
      double head_ratio = max_head_pair != 0 ? (max_head_pair - min_head_pair) / max_head_pair : (double)INT_MAX;
      // cout << volume << " " << greenfield_wall_height << " " <<
      // greenfield->elevation << " " << pit_depth << " " << pit->elevation << " "
      // << head << " " << head_ratio << "\n";

      if (volume < find_required_volume(energy_capacity, head)) {
        continue;
      }
      required_volume = volume;

      if (head_ratio > max_head_variability) {
        break;
      }

      if(greenfield.pit)
        greenfield.fill_depth_from_MOL = linear_interpolate(required_volume, greenfield.volumes, int_to_double_vector(greenfield.fill_depths));
      
      return true;
    }

    pit.elevation += pit_height_resolution;
  }
  return false;
}

Pair *check_good_pair(RoughReservoir* upper, RoughReservoir* lower,
                      double energy_capacity, int storage_time, Pair *pair,
                      int max_FOM) {
  int head = upper->elevation - lower->elevation;
  double required_volume = find_required_volume(energy_capacity, head);
  if ((max(upper->volumes) < required_volume) ||
      (max(lower->volumes) < required_volume * (lower->river ? 5 : 1))) {
    return NULL;
  }
  ExistingPit single_pit;
  RoughReservoir greenfield = *upper;
  RoughReservoir pit = *lower;
  if (upper->pit) {
      greenfield = *lower;
      pit = *upper;
    }

  if (search_config.search_type == SearchType::SINGLE_PIT) {
    single_pit = single_pit_details;
  }

  if (search_config.search_type == SearchType::SINGLE_PIT) {
    if (!determine_pit_elevation_and_volume(upper, lower, energy_capacity,
                                          single_pit, required_volume, head)) {
      return NULL;
    }
  } else if (search_config.search_type == SearchType::BULK_PIT) {
    if (!determine_pit_elevation_and_volume(pit, greenfield, energy_capacity, required_volume, head)) {
      return NULL;
    }
  }

  double upper_dam_wall_height = 0;
  double lower_dam_wall_height = 0;
  double upper_water_rock_estimate = INF;
  double lower_water_rock_estimate = INF;

  if (!upper->brownfield) {
    upper_dam_wall_height =
        linear_interpolate(required_volume, upper->volumes, dam_wall_heights);
    upper_water_rock_estimate =
        required_volume / linear_interpolate(upper_dam_wall_height,
                                             dam_wall_heights,
                                             upper->dam_volumes);
  } else {
    if (search_config.search_type == SearchType::SINGLE_PIT)
      upper_dam_wall_height = linear_interpolate(required_volume +
                                 pit_volume(single_pit,
                                            single_pit.reservoir.elevation,
                                            upper->elevation),
                             get_volumes(single_pit),
                             int_to_double_vector(get_altitudes(single_pit))) -
          upper->elevation;
    else
      upper_dam_wall_height = dam_wall_heights[0];
    upper_water_rock_estimate = INF;
  }

  if (!lower->brownfield && !lower->ocean) {
    lower_dam_wall_height =
        linear_interpolate(required_volume, lower->volumes, dam_wall_heights);
    lower_water_rock_estimate =
        required_volume / linear_interpolate(lower_dam_wall_height,
                                             dam_wall_heights,
                                             lower->dam_volumes);
  } else {
    if (search_config.search_type == SearchType::SINGLE_PIT)
      lower_dam_wall_height =
          linear_interpolate(required_volume +
                                 pit_volume(single_pit,
                                            single_pit.reservoir.elevation,
                                            lower->elevation),
                             get_volumes(single_pit),
                             int_to_double_vector(get_altitudes(single_pit))) -
          lower->elevation;
    else
      lower_dam_wall_height = dam_wall_heights[0];
    lower_water_rock_estimate = INF;
  }

  if ((!upper->brownfield && upper_dam_wall_height > upper->max_dam_height) ||
      (!lower->brownfield && !lower->ocean &&
       lower_dam_wall_height > lower->max_dam_height)) {
    return NULL;
  }

  if ((upper_water_rock_estimate * lower_water_rock_estimate) <
      min_pair_water_rock *
          (upper_water_rock_estimate + lower_water_rock_estimate)) {
    return NULL;
  }

  ArrayCoordinate upper_coordinates = upper->pour_point;
  ArrayCoordinate lower_coordinates = lower->pour_point;

  double least_distance = find_least_distance_sqd(
      upper, lower, upper_dam_wall_height,
      lower_dam_wall_height, &upper_coordinates, &lower_coordinates);

  if (SQ(head * 0.001) < least_distance * SQ(min_slope)) {
    return NULL;
  }

  if(search_config.search_type==SearchType::OCEAN){
    lower->pour_point=lower_coordinates;
  }

  Reservoir upper_reservoir = Reservoir_init(upper->pour_point, upper->pit ? pit.elevation : upper->elevation, upper->elevation);
  Reservoir lower_reservoir = Reservoir_init(lower->pour_point, (upper->pit && lower->pit) ? greenfield.elevation : (lower->pit ? pit.elevation : lower->elevation), lower->elevation);

  upper_reservoir.identifier = upper->identifier;
  upper_reservoir.volume = required_volume;

  if (!upper->brownfield) {
    upper_reservoir.dam_volume = linear_interpolate(
        upper_dam_wall_height, dam_wall_heights, upper->dam_volumes);
    upper_reservoir.area = linear_interpolate(upper_dam_wall_height,
                                              dam_wall_heights, upper->areas);
  } else if (upper->pit) {
    upper_reservoir.area = linear_interpolate(required_volume, upper->volumes, upper->areas);
  } else {
    upper_reservoir.area = upper->areas[0];
  }
  upper_reservoir.water_rock = upper_water_rock_estimate;
  upper_reservoir.dam_height = upper_dam_wall_height;
  upper_reservoir.max_dam_height = upper->max_dam_height;
  if(!upper->pit)
    upper_reservoir.fill_depth = linear_interpolate(required_volume, upper->volumes, int_to_double_vector(upper->fill_depths));
  else
    upper_reservoir.fill_depth = pit.fill_depth_from_MOL;
  upper_reservoir.brownfield = upper->brownfield;
  upper_reservoir.river = upper->river;
  upper_reservoir.pit = upper->pit;

  lower_reservoir.identifier = lower->identifier;
  lower_reservoir.volume = required_volume;

  if (!lower->brownfield) {
    lower_reservoir.dam_volume = linear_interpolate(
        lower_dam_wall_height, dam_wall_heights, lower->dam_volumes);
    lower_reservoir.area = linear_interpolate(lower_dam_wall_height,
                                              dam_wall_heights, lower->areas);
  } else if (lower->pit) {
    lower_reservoir.area = linear_interpolate(required_volume, lower->volumes, lower->areas);
  } else {
    lower_reservoir.area = lower->areas[0];
  }
  lower_reservoir.water_rock = lower_water_rock_estimate;
  lower_reservoir.dam_height = lower_dam_wall_height;
  lower_reservoir.max_dam_height = lower->max_dam_height;
  if(!lower->pit)
    lower_reservoir.fill_depth = linear_interpolate(required_volume, lower->volumes, int_to_double_vector(lower->fill_depths));
  else
    lower_reservoir.fill_depth = upper->pit ? greenfield.fill_depth_from_MOL : pit.fill_depth_from_MOL;
  lower_reservoir.brownfield = lower->brownfield;
  lower_reservoir.river = lower->river;
  lower_reservoir.pit = lower->pit;
  lower_reservoir.ocean = lower->ocean;

  pair->identifier = upper->identifier + " & " + lower->identifier;
  pair->upper = upper_reservoir;
  pair->lower = lower_reservoir;
  pair->head = head;
  pair->distance = SQRT(least_distance);
  pair->pp_distance =
      find_distance(pair->upper.pour_point, pair->lower.pour_point);
  pair->energy_capacity = energy_capacity;
  pair->storage_time = storage_time;
  pair->required_volume = required_volume;
  pair->slope = pair->head / (pair->distance) * 0.001;
  pair->water_rock =
      1 / (1 / pair->upper.water_rock + 1 / pair->lower.water_rock);

  set_FOM(pair);

  if (pair->FOM > max_FOM)
    return NULL;
  return pair;
}

void pairing(vector<unique_ptr<RoughReservoir>> &upper_reservoirs,
             vector<unique_ptr<RoughReservoir>> &lower_reservoirs, FILE *csv_file,
             FILE *csv_data_file, bool existing_existing_allowed=false) {
  vector<set<Pair>> temp_pairs;
  for (uint itest = 0; itest < tests.size(); itest++) {
    pairs.push_back(0);
    set<Pair> a;
    temp_pairs.push_back(a);
  }

  for (uint iupper = 0; iupper < upper_reservoirs.size(); iupper++) {
    RoughReservoir* upper_reservoir = upper_reservoirs[iupper].get();
    double coslat = COS(RADIANS(upper_reservoir->latitude));
    for (uint ilower = 0; ilower < lower_reservoirs.size(); ilower++) {
      RoughReservoir* lower_reservoir = lower_reservoirs[ilower].get();
      int head = upper_reservoir->elevation - lower_reservoir->elevation;
      if (!upper_reservoir->river && !lower_reservoir->river)
        if (head < min_head || head > max_head)
          continue;

      if (!existing_existing_allowed && upper_reservoir->brownfield && lower_reservoir->brownfield)
        continue;

      // Pour point separation
      double min_dist_sqd = find_distance_sqd(
          upper_reservoir->pour_point, lower_reservoir->pour_point, coslat);

      if(upper_reservoir->brownfield){
        min_dist_sqd = INF;
        RoughBfieldReservoir* br = static_cast<RoughBfieldReservoir*>(upper_reservoir);
        for(size_t i = 0; i<br->shape_bound.size(); i++){
          ArrayCoordinate ac = br->shape_bound[i];
          double dist_sqd = find_distance_sqd(ac, lower_reservoir->pour_point, coslat);
          if(dist_sqd < min_dist_sqd){
            min_dist_sqd = dist_sqd;
          }
        }
      }
      if(lower_reservoir->brownfield || lower_reservoir->ocean){
        min_dist_sqd = INF;
        RoughBfieldReservoir* lr = static_cast<RoughBfieldReservoir*>(lower_reservoir);

        int idx = 0;
        for(size_t i = 0; i<lr->shape_bound.size(); i++){
          ArrayCoordinate ac = lr->shape_bound[i];
          double dist_sqd = find_distance_sqd(ac, upper_reservoir->pour_point, coslat);
          if(dist_sqd < min_dist_sqd){
            idx = i;
            min_dist_sqd = dist_sqd;
          }
        }
        if(lower_reservoir->river){
          lower_reservoir->elevation = lr->elevations[idx];
          lower_reservoir->pour_point = lr->shape_bound[idx];
        }
      }

      if(upper_reservoir->river)
        continue;

      head = upper_reservoir->elevation - lower_reservoir->elevation;
      if (head < min_head || head > max_head)
        continue;

      if (SQ(head * 0.001) <= min_dist_sqd * SQ(min_pp_slope))
        continue;


      for (uint itest = 0; itest < tests.size(); itest++) {
        Pair temp_pair;
        int max_FOM =
            (category_cutoffs[0].storage_cost * tests[itest].storage_time +
             category_cutoffs[0].power_cost) *
            (1 + tolerance_on_FOM);

        if (check_good_pair(upper_reservoir, lower_reservoir,
                            tests[itest].energy_capacity,
                            tests[itest].storage_time, &temp_pair, max_FOM)) {
          temp_pairs[itest].insert(temp_pair);

          if ((int)temp_pairs[itest].size() > max_lowers_per_upper ||
              ((search_config.search_type == SearchType::BULK_PIT || search_config.search_type == SearchType::SINGLE_PIT)&&
              temp_pairs[itest].size() > 1))
            temp_pairs[itest].erase(prev(temp_pairs[itest].end()));
        }
      }
    }

    for (uint itest = 0; itest < tests.size(); itest++) {
      for (Pair pair : temp_pairs[itest]) {
        write_rough_pair_csv(csv_file, &pair);
        write_rough_pair_data(csv_data_file, &pair);
        pairs[itest]++;
      }
      temp_pairs[itest].clear();
    }
  }
}

int main(int nargs, char **argv) {
  search_config = SearchConfig(nargs, argv);
  cout << "Pairing started for " << search_config.filename() << endl;

  unsigned long t_usec = walltime_usec();
  parse_variables(convert_string("storage_location"));
  parse_variables(convert_string(file_storage_location + "variables"));

  vector<unique_ptr<RoughReservoir>> upper_reservoirs;
  vector<unique_ptr<RoughReservoir>> lower_reservoirs;
  
  if (search_config.search_type.existing()) {
    if (search_config.search_type.single())
      search_config.grid_square = get_square_coordinate(get_existing_reservoir(search_config.name));
    upper_reservoirs = read_rough_reservoir_data(
        convert_string(file_storage_location + "processing_files/reservoirs/" +
                       search_config.filename() + "_reservoirs_data.csv"));
    
    if (search_config.search_type == SearchType::BULK_EXISTING || search_config.search_type == SearchType::BULK_PIT){
      
      vector<unique_ptr<RoughReservoir>> greenfield_reservoirs = read_rough_reservoir_data(
          convert_string(file_storage_location + "processing_files/reservoirs/" +
                         str(search_config.grid_square) + "_reservoirs_data.csv"));
      for(size_t i = 0; i<greenfield_reservoirs.size(); i++){
        upper_reservoirs.push_back(std::move(greenfield_reservoirs[i]));
      }
    }
    if (search_config.search_type == SearchType::SINGLE_PIT) {
      single_pit_details = get_pit_details(search_config.name);
    }
  } else
    upper_reservoirs = read_rough_reservoir_data(
        convert_string(file_storage_location + "processing_files/reservoirs/" +
                       str(search_config.grid_square) + "_reservoirs_data.csv"));
  
  GridSquare neighbors[9] = {
      (GridSquare){search_config.grid_square.lat, search_config.grid_square.lon},
      (GridSquare){search_config.grid_square.lat + 1, search_config.grid_square.lon - 1},
      (GridSquare){search_config.grid_square.lat + 1, search_config.grid_square.lon},
      (GridSquare){search_config.grid_square.lat + 1, search_config.grid_square.lon + 1},
      (GridSquare){search_config.grid_square.lat, search_config.grid_square.lon + 1},
      (GridSquare){search_config.grid_square.lat - 1, search_config.grid_square.lon + 1},
      (GridSquare){search_config.grid_square.lat - 1, search_config.grid_square.lon},
      (GridSquare){search_config.grid_square.lat - 1, search_config.grid_square.lon - 1},
      (GridSquare){search_config.grid_square.lat, search_config.grid_square.lon - 1}};

  set<string> lower_ids;
  for (int i = 0; i < 9; i++) {
    try {
      vector<unique_ptr<RoughReservoir>> temp = read_rough_reservoir_data(convert_string(
          file_storage_location + "processing_files/reservoirs/" +
          search_config.search_type.lowers_prefix() + str(neighbors[i]) + "_reservoirs_data.csv"));

      for (uint j = 0; j < temp.size(); j++) {
        if ((search_config.search_type == SearchType::BULK_EXISTING || search_config.search_type == SearchType::BULK_PIT) &&
            lower_ids.contains(temp[j]->identifier))
          continue;

        lower_ids.insert(temp[j]->identifier);
        lower_reservoirs.push_back(std::move(temp[j]));
      }
    } catch (int e) {
      search_config.logger.debug("Could not import reservoirs from " + file_storage_location +
                                 "processing_files/reservoirs/" +
                                 search_config.search_type.lowers_prefix() + str(neighbors[i]) +
                                 "_reservoirs_data.csv");
    }
  }
  search_config.logger.debug("Read in "+to_string(upper_reservoirs.size())+" uppers");
  search_config.logger.debug("Read in " + to_string(lower_reservoirs.size()) + " lowers");

  mkdir(convert_string(file_storage_location + "output/pairs"), 0777);
  FILE *csv_file =
      fopen(convert_string(file_storage_location + "output/pairs/" +
                           search_config.filename() + "_rough_pairs.csv"),
            "w");
  if (!csv_file) {
    fprintf(stderr, "Failed to open reservoir pair CSV file\n");
    exit(1);
  }
  write_rough_pair_csv_header(csv_file);

  mkdir(convert_string(file_storage_location + "processing_files/pairs"), 0777);
  FILE *csv_data_file =
      fopen(convert_string(file_storage_location + "processing_files/pairs/" +
                           search_config.filename() + "_rough_pairs_data.csv"),
            "w");
  if (!csv_data_file) {
    fprintf(stderr, "Failed to open reservoir pair CSV data file\n");
    exit(1);
  }
  write_rough_pair_data_header(csv_data_file);

  pairing(upper_reservoirs, lower_reservoirs, csv_file, csv_data_file, true);
  if (search_config.search_type.existing())
    pairing(lower_reservoirs, upper_reservoirs, csv_file, csv_data_file);

  int total = 0;
  for (uint itest = 0; itest < tests.size(); itest++) {
    search_config.logger.debug(to_string(pairs[itest]) + " " +
                               to_string(tests[itest].energy_capacity) +
                               "GWh " + to_string(tests[itest].storage_time) +
                               "h pairs");
    total += pairs[itest];
  }

  fclose(csv_file);
  fclose(csv_data_file);

  cout << "Pairing finished for " << search_config.filename() << ". Found "
       << total << " pairs. Runtime: " << 1.0e-6 * (walltime_usec() - t_usec)
       << " sec" << endl;
}
