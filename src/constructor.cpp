#include "coordinates.h"
#include "csv.h"
#include "mining_pits.h"
#include "model2D.h"
#include "phes_base.h"
#include "constructor_helpers.hpp"
#include "kml.h"
#include "search_config.hpp"

vector<vector<Pair>> pairs;

bool model_existing_reservoir(Reservoir* reservoir, Reservoir_KML_Coordinates* coordinates, vector<vector<vector<vector<GeographicCoordinate>>>>& countries, vector<string>& country_names){
  if(!reservoir->river){
    ExistingReservoir r;
    if (use_tiled_bluefield)
      r = get_existing_tiled_reservoir(reservoir->identifier, reservoir->latitude,
          reservoir->longitude);
    else
      r = get_existing_reservoir(reservoir->identifier);
    reservoir->volume = r.volume;
    string polygon_string = str(compress_poly(corner_cut_poly(r.polygon)), reservoir->pit ? r.elevation+reservoir->dam_height : r.elevation+5);
    coordinates->reservoir = polygon_string;

    if (search_config.search_type.single()) 
      search_config.grid_square = get_square_coordinate(get_existing_reservoir(search_config.name));
    GeographicCoordinate origin = get_origin(search_config.grid_square, border);
    reservoir->shape_bound.clear();
    for(GeographicCoordinate p : r.polygon){
      reservoir->shape_bound.push_back(convert_coordinates(p, origin));
      reservoir->reservoir_polygon.push_back(convert_coordinates(p, origin));
    }

    //KML
    reservoir->country = find_country(GeographicCoordinate_init(reservoir->latitude, reservoir->longitude), countries, country_names);
  }
    return true;
}

bool model_pair(Pair *pair, Pair_KML *pair_kml, Model<bool> *seen,
                bool *non_overlap, int max_FOM, BigModel big_model,
                Model<char> *full_cur_model,
                vector<vector<vector<vector<GeographicCoordinate>>>> &countries,
                vector<string> &country_names, vector<vector<GeographicCoordinate>> &seen_polygons) {

  vector<ArrayCoordinate> used_points;
  *non_overlap = true;

  // Remove pairs not within grid_square. Required for overlap analysis of existing reservoirs/pits
  if(search_config.search_type.existing()){
    if((pair->upper.brownfield && !check_within(convert_coordinates(pair->upper.pour_point), search_config.grid_square)) || (pair->lower.brownfield && !check_within(convert_coordinates(pair->lower.pour_point), search_config.grid_square)))
      return false;
  }

  if (pair->upper.brownfield && (search_config.search_type != SearchType::BULK_PIT)) {
    if (!model_existing_reservoir(&pair->upper, &pair_kml->upper, countries,
                                  country_names))
      return false;
  } else if (pair->upper.brownfield && (search_config.search_type == SearchType::BULK_PIT)) {
    if (!model_bulk_pit(&pair->upper, &pair_kml->upper,
                              countries, country_names)){
      return false;
    }
  } else if (!model_reservoir(&pair->upper, &pair_kml->upper, seen, non_overlap,
                              &used_points, big_model, full_cur_model,
                              countries, country_names))
    return false;

  if (pair->lower.brownfield && (search_config.search_type != SearchType::BULK_PIT)) {
    if (!model_existing_reservoir(&pair->lower, &pair_kml->lower, countries,
                                  country_names))
      return false;
  } else if (pair->lower.brownfield && (search_config.search_type == SearchType::BULK_PIT)){
    if (!model_bulk_pit(&pair->lower, &pair_kml->lower, 
                              countries, country_names)){
      return false;
    }
  
  } else if (!pair->lower.ocean &&
             !model_reservoir(&pair->lower, &pair_kml->lower, seen, non_overlap,
                              &used_points, big_model, full_cur_model,
                              countries, country_names))
    return false;

  pair->country = pair->upper.country;

  ArrayCoordinate upper_closest_point = pair->upper.pour_point;
  ArrayCoordinate lower_closest_point = pair->lower.pour_point;
  double mindist = find_distance(upper_closest_point, lower_closest_point);
  for (ArrayCoordinate u_bound : pair->upper.shape_bound)
    for (ArrayCoordinate l_bound : pair->lower.shape_bound)
      if (find_distance(u_bound, l_bound) < mindist) {
        upper_closest_point = u_bound;
        lower_closest_point = l_bound;
        mindist = find_distance(u_bound, l_bound);
      }

  pair->distance = mindist;
  pair->slope = pair->head / (pair->distance * 1000);
  pair->volume = min(pair->upper.volume, pair->lower.volume);
  pair->water_rock =
      1 / ((1 / pair->upper.water_rock) + (1 / pair->lower.water_rock));
  set_FOM(pair);
  if (pair->FOM > max_FOM || pair->category == 'Z') {
    return false;
  }

  // Check overlap between reservoirs during pit and existing reservoir constructor
  if (pair->upper.brownfield || pair->lower.brownfield) {
    *non_overlap = true;

    // Overlap between upper and lower reservoirs, and overlap with other sites
    for (ArrayCoordinate ac : pair->lower.reservoir_polygon){
      if(!*non_overlap)
        break;

      if(check_within(convert_coordinates(ac), compress_poly(corner_cut_poly(convert_poly(pair->upper.reservoir_polygon)))))
        return false;

      *non_overlap = (!check_within(convert_coordinates(ac), seen_polygons) && *non_overlap);
    }
    for (ArrayCoordinate ac : pair->upper.reservoir_polygon){
      if(!*non_overlap)
        break;

      if(check_within(convert_coordinates(ac), compress_poly(corner_cut_poly(convert_poly(pair->lower.reservoir_polygon)))))
        return false;

      *non_overlap = (!check_within(convert_coordinates(ac), seen_polygons) && *non_overlap);
    }

    if(pair->upper.pit)
      seen_polygons.push_back(compress_poly(corner_cut_poly(convert_poly(pair->upper.reservoir_polygon))));

    if(pair->lower.pit)
      seen_polygons.push_back(compress_poly(corner_cut_poly(convert_poly(pair->lower.reservoir_polygon))));
  }

  if (*non_overlap) {
    for (uint i = 0; i < used_points.size(); i++) {
      seen->set(used_points[i].row, used_points[i].col, true);
    }
    pair->non_overlap = 1;
  } else {
    pair->non_overlap = 0;
  }

    GeographicCoordinate average = GeographicCoordinate_init((convert_coordinates(upper_closest_point).lat+convert_coordinates(lower_closest_point).lat)/2,
    	((convert_coordinates(upper_closest_point).lon+convert_coordinates(lower_closest_point).lon)/2));
    pair_kml->point = dtos(average.lon,5)+","+dtos(average.lat,5)+",0";
    pair_kml->line = dtos(convert_coordinates(upper_closest_point).lon,5)+","+dtos(convert_coordinates(upper_closest_point).lat,5)+",0 "+dtos(convert_coordinates(lower_closest_point).lon,5)+","+dtos(convert_coordinates(lower_closest_point).lat,5)+",0";
	return true;
}

int main(int nargs, char **argv)
{
    search_config = SearchConfig(nargs, argv);

    cout << "Constructor started for " << search_config.filename() << endl;

    GDALAllRegister();
    parse_variables(convert_string("storage_location"));
    parse_variables(convert_string(file_storage_location+"variables"));
    unsigned long t_usec = walltime_usec();

    pairs = read_rough_pair_data(convert_string(file_storage_location+"processing_files/pretty_set_pairs/"+search_config.filename()+"_rough_pretty_set_pairs_data.csv"));

    mkdir(convert_string(file_storage_location+"output/final_output_classes"), 0777);
    mkdir(convert_string(file_storage_location+"output/final_output_classes/"+search_config.filename()),0777);
    mkdir(convert_string(file_storage_location+"output/final_output_FOM"), 0777);
    mkdir(convert_string(file_storage_location+"output/final_output_FOM/"+search_config.filename()),0777);

    FILE *total_csv_file_classes = fopen(convert_string(file_storage_location+"output/final_output_classes/"+search_config.filename()+"/"+search_config.filename()+"_total.csv"), "w");
    write_summary_csv_header(total_csv_file_classes);
    FILE *total_csv_file_FOM = fopen(convert_string(file_storage_location+"output/final_output_FOM/"+search_config.filename()+"/"+search_config.filename()+"_total.csv"), "w");
    write_summary_csv_header(total_csv_file_FOM);

    // Read in all rough pairs from neighbours. Only keep pairs that have a pit/existing reservoir in the centre grid square
    // Rough pairs CSV will not match constructor pairs for a given grid square, but this is required for overlap analysis
    if(search_config.search_type.existing()){
      for (uint d=0; d<directions.size(); d++){   
        GridSquare neighbor_gs = GridSquare_init(search_config.grid_square.lat + directions[d].row, search_config.grid_square.lon+ directions[d].col);

        vector<vector<Pair>> neighbor_pairs;
        try {
          neighbor_pairs = read_rough_pair_data(convert_string(file_storage_location+"processing_files/pretty_set_pairs/"+search_config.search_type.prefix()+str(neighbor_gs) +"_rough_pretty_set_pairs_data.csv"));
          for(uint i=0; i<neighbor_pairs.size(); i++){
            pairs.push_back(neighbor_pairs[i]);
          }
        } catch(int e) {
          search_config.logger.debug("Could not import pretty set pairs from " + file_storage_location +
                                 "processing_files/pretty_set_pairs/" +
                                 search_config.search_type.prefix()+str(neighbor_gs) +"_rough_pretty_set_pairs_data.csv");
        }
      }
    }

    uint total_pairs = 0;
	for(uint i = 0; i<pairs.size(); i++)
		total_pairs += pairs[i].size();

	if (total_pairs == 0) {
        for(uint i = 0; i<tests.size(); i++){
            write_summary_csv(total_csv_file_classes, str(search_config.grid_square), str(tests[i]),
                                0, 0, 0);
            write_summary_csv(total_csv_file_FOM, str(search_config.grid_square), str(tests[i]),
                                0, 0, 0);
        }
        write_summary_csv(total_csv_file_classes, str(search_config.grid_square), "TOTAL", 0, -1, 0);
        write_summary_csv(total_csv_file_FOM, str(search_config.grid_square), "TOTAL", 0, -1, 0);
        fclose(total_csv_file_classes);
        fclose(total_csv_file_FOM);
        cout << "Constructor finished for " << convert_string(search_config.filename()) << ". Found no pairs. Runtime: " << 1.0e-6*(walltime_usec() - t_usec) << " sec" << endl;
        return 0;
    }

    if(search_config.search_type.single())
        search_config.grid_square = get_square_coordinate(get_existing_reservoir(search_config.name));

    BigModel big_model = BigModel_init(search_config.grid_square);
    vector<string> country_names;
    vector<vector<vector<vector<GeographicCoordinate>>>> countries = read_countries(file_storage_location+"input/countries/countries.txt", country_names);

    Model<bool>* seen = new Model<bool>(big_model.DEM->nrows(), big_model.DEM->nrows(), MODEL_SET_ZERO);
    seen->set_geodata(big_model.DEM->get_geodata());
    Model<char>* full_cur_model = new Model<char>(big_model.DEM->nrows(), big_model.DEM->ncols(), MODEL_SET_ZERO);
    full_cur_model->set_geodata(big_model.DEM->get_geodata());
    vector<vector<GeographicCoordinate>> seen_polygons;

    int total_count = 0;
    int total_capacity = 0;
    for(uint i = 0; i<tests.size(); i++){
      int count = 0;
      int non_overlapping_count = 0;
      if (pairs[i].size() != 0) {
        // Extract bulk pit shape bounds found during screening
        if (search_config.search_type == SearchType::BULK_PIT) {  
          read_pit_polygons(convert_string(file_storage_location + "processing_files/reservoirs/pit_" +
                          str(search_config.grid_square) + "_reservoirs_data.csv"), pairs[i], search_config.grid_square);

          for (uint d=0; d<directions.size(); d++){   
            GridSquare neighbor_gs = GridSquare_init(search_config.grid_square.lat + directions[d].row, search_config.grid_square.lon+ directions[d].col);
            read_pit_polygons(convert_string(file_storage_location + "processing_files/reservoirs/pit_" +
                          str(neighbor_gs) + "_reservoirs_data.csv"), pairs[i], neighbor_gs);
          }
        }

        FILE *csv_file_classes = fopen(convert_string(file_storage_location+"output/final_output_classes/"+search_config.filename()+"/"+search_config.filename()+"_"+str(tests[i])+".csv"), "w");
        write_pair_csv_header(csv_file_classes, false);
        FILE *csv_file_FOM = fopen(convert_string(file_storage_location+"output/final_output_FOM/"+search_config.filename()+"/"+search_config.filename()+"_"+str(tests[i])+".csv"), "w");
        write_pair_csv_header(csv_file_FOM, true);

        ofstream kml_file_classes(convert_string(file_storage_location+"output/final_output_classes/"+search_config.filename()+"/"+search_config.filename()+"_"+str(tests[i])+".kml"), ios::out);
        ofstream kml_file_FOM(convert_string(file_storage_location+"output/final_output_FOM/"+search_config.filename()+"/"+search_config.filename()+"_"+str(tests[i])+".kml"), ios::out);
        KML_Holder kml_holder;

        sort(pairs[i].begin(), pairs[i].end());
	      set<string> lowers;
        set<string> uppers;
        bool keep_lower;
        bool keep_upper;
        for(uint j=0; j<pairs[i].size(); j++){
            Pair_KML pair_kml;
            bool non_overlap;
            int max_FOM = category_cutoffs[0].storage_cost*tests[i].storage_time+category_cutoffs[0].power_cost;

            if(model_pair(&pairs[i][j], &pair_kml, seen, &non_overlap, max_FOM, big_model, full_cur_model, countries, country_names, seen_polygons)){
                write_pair_csv(csv_file_classes, &pairs[i][j], false);
                write_pair_csv(csv_file_FOM, &pairs[i][j], true);
                keep_lower = !lowers.contains(pairs[i][j].lower.identifier);
                keep_upper = !uppers.contains(pairs[i][j].upper.identifier);
                lowers.insert(pairs[i][j].lower.identifier);
                uppers.insert(pairs[i][j].upper.identifier);
                update_kml_holder(&kml_holder, &pairs[i][j], &pair_kml, keep_upper, keep_lower);
                count++;
                if(non_overlap){
                    non_overlapping_count++;
                    total_count++;
                    total_capacity+=tests[i].energy_capacity;
                }
            }
        }
        kml_file_classes << output_kml(&kml_holder, search_config.filename(), tests[i]);
        kml_file_FOM << output_kml(&kml_holder, search_config.filename(), tests[i]);
        search_config.logger.debug(to_string(count) + " " + to_string(tests[i].energy_capacity) + "GWh "+to_string(tests[i].storage_time) + "h Pairs");
        kml_file_classes.close();
        kml_file_FOM.close();
        fclose(csv_file_classes);
        fclose(csv_file_FOM);
      }
      write_summary_csv(total_csv_file_classes, str(search_config.grid_square), str(tests[i]),
                        non_overlapping_count, count, count*tests[i].energy_capacity);
      write_summary_csv(total_csv_file_FOM, str(search_config.grid_square), str(tests[i]),
                        non_overlapping_count, count, count*tests[i].energy_capacity);
    }
    write_summary_csv(total_csv_file_classes, str(search_config.grid_square), "TOTAL", total_count, -1, total_capacity);
    write_summary_csv(total_csv_file_FOM, str(search_config.grid_square), "TOTAL", total_count, -1, total_capacity);
    fclose(total_csv_file_classes);
    fclose(total_csv_file_FOM);
    cout << "Constructor finished for " << convert_string(search_config.filename()) << ". Found " << total_count << " non-overlapping pairs with a total of " << total_capacity << "GWh. Runtime: " << 1.0e-6*(walltime_usec() - t_usec) << " sec" << endl;
}
