#include "coordinates.h"
#include "kml.h"
#include "phes_base.h"
#include "reservoir.h"
#include "mining_pits.h"
#include "search_config.hpp"

string ReplaceAll(string str, const string &from, const string &to) {
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos +=
        to.length(); // Handles case where 'to' is a substring of 'from'
  }
  return str;
}

void write_to_csv_file(FILE *csv_file, vector<string> cols) {
  for (uint i = 0; i < cols.size(); i++) {
    if (cols[i].find(',') != std::string::npos) {
      cols[i] = ReplaceAll(cols[i], string("\""), string("\"\""));
      cols[i] = ReplaceAll(cols[i], string("  "), string(""));
      cols[i] = ReplaceAll(cols[i], string("\n"), string(""));
      cols[i] = '"' + cols[i] + '"';
    }
    fprintf(csv_file, "%s", cols[i].c_str());
    if (i != cols.size() - 1)
      fprintf(csv_file, ",");
  }
  fprintf(csv_file, "\n");
}

vector<string> read_from_csv_file(string line) {
  return read_from_csv_file(line, ',');
}

vector<string> read_from_csv_file(string line, char delimeter) {
  line.erase(remove(line.begin(), line.end(), '\n'), line.end());
  line.erase(remove(line.begin(), line.end(), '\r'), line.end());
  vector<string> cols;
  string cur = "";
  bool instring = false;
  for (uint i = 0; i < line.length(); i++) {
    if (line[i] == '"')
      instring = !instring;
    else if (!instring && line[i] == delimeter) {
      cols.push_back(cur);
      cur = "";
    } else
      cur.push_back(line[i]);
  }
  cols.push_back(cur);
  return cols;
}

vector<ExistingReservoir> read_existing_reservoir_data(const char *filename) {
  vector<ExistingReservoir> reservoirs;
  ifstream inputFile(filename);
  bool header = true;
  string s;
  while (getline(inputFile, s)) {
    if (header) {
      header = false;
      continue;
    }
    vector<string> line = read_from_csv_file(s);
    ExistingReservoir reservoir = ExistingReservoir_init(
        line[0], stod(line[1]), stod(line[2]), stod(line[3]), stod(line[4]));
    reservoirs.push_back(reservoir);
  }
  if (header) {
    cout << "CSV file " << filename << " is empty." << endl;
    throw 1;
  }

  return reservoirs;
}

vector<ExistingPit> read_existing_pit_data(char *filename) {
  vector<ExistingPit> pits;
  ifstream inputFile(filename);
  bool header = true;
  string s;
  while (getline(inputFile, s)) {
    if (header) {
      header = false;
      continue;
    }
    vector<string> line = read_from_csv_file(s);
    //printf("%s\n",convert_string(line[0]));
    ExistingReservoir reservoir = ExistingReservoir_init(
        line[0], stod(line[1]), stod(line[2]), stod(line[3]), stod(line[4]));
    ExistingPit pit = ExistingPit_init(reservoir);
    AltitudeVolumePair bottom;
    bottom.altitude = reservoir.elevation;
    bottom.volume = 0;
    pit.volumes.push_back(bottom);
    for (int i = 0; i < num_altitude_volume_pairs; i++) {
      AltitudeVolumePair pair;
      pair.altitude = stoi(line[5 + 2 * i]);
      pair.volume = stod(line[6 + 2 * i]);
      pit.volumes.push_back(pair);
    }
    sort(pit.volumes.begin(), pit.volumes.end());
    pits.push_back(pit);
  }
  if (header) {
    throw 1;
  }

  return pits;
}

vector<string> read_names(char *filename) {
  vector<string> names;
  ifstream inputFile(filename);
  bool header = true;
  string s;
  while (getline(inputFile, s)) {
    if (header) {
      header = false;
      continue;
    }
    names.push_back(read_from_csv_file(s)[0]);
  }
  if (header) {
    cout << "CSV file " << filename << " is empty." << endl;
    throw 1;
  }
  return names;
}

void write_rough_reservoir_csv_header(FILE *csv_file) {
  vector<string> header = {"Identifier",         "Latitude",
                           "Longitude",          "Elevation (m)",
                           "Maximum dam height", "Watershed area (ha)"};
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back(dtos(dam_wall_heights[i], 0) + "m reservoir volume (GL)");
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back(dtos(dam_wall_heights[i], 0) + "m reservoir area (ha)");
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back(dtos(dam_wall_heights[i], 0) + "m dam volume (GL)");
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back(dtos(dam_wall_heights[i], 0) + "m water to rock");
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back("Test " + str(i+1) + " fill depth");
  write_to_csv_file(csv_file, header);
}

void write_rough_reservoir_data_header(FILE *csv_file) {
  vector<string> header = {"Identifier",         "Latitude",
                           "Longitude",          "Elevation (m)",
                           "Maximum dam height", "Watershed area (ha)"};
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back(dtos(dam_wall_heights[i], 0) + "m reservoir volume (GL)");
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back(dtos(dam_wall_heights[i], 0) + "m reservoir area (ha)");
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back(dtos(dam_wall_heights[i], 0) + "m dam volume (GL)");
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    header.push_back("Test " + str(i+1) + " fill depth (m)");
  header.push_back("Brownfield");
  header.push_back("Ocean");
  header.push_back("Turkey");
  header.push_back("Boundary");
  write_to_csv_file(csv_file, header);
}

void write_rough_reservoir_csv(FILE *csv_file, RoughReservoir *reservoir) {
  vector<string> line = {reservoir->identifier,
                         dtos(reservoir->latitude, 4),
                         dtos(reservoir->longitude, 4),
                         to_string(reservoir->elevation),
                         dtos(reservoir->max_dam_height, 0),
                         dtos(reservoir->watershed_area, 1)};
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    line.push_back(dtos(reservoir->volumes[i], 2));
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    line.push_back(dtos(reservoir->areas[i], 2));
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    line.push_back(dtos(reservoir->dam_volumes[i], 2));
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    line.push_back(dtos(reservoir->water_rocks[i], 1));
  for (uint i = 0; i < dam_wall_heights.size(); i++) {
    line.push_back(to_string(reservoir->fill_depths[i]));
  }
    
  write_to_csv_file(csv_file, line);
}

void write_rough_reservoir_data(FILE *csv_file, RoughReservoir *reservoir) {
  vector<string> line = {reservoir->identifier,
                         dtos(reservoir->latitude, 5),
                         dtos(reservoir->longitude, 5),
                         to_string(reservoir->elevation),
                         dtos(reservoir->max_dam_height, 1),
                         dtos(reservoir->watershed_area, 2)};
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    line.push_back(dtos(reservoir->volumes[i], 5));
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    line.push_back(dtos(reservoir->areas[i], 2));
  for (uint i = 0; i < dam_wall_heights.size(); i++)
    line.push_back(dtos(reservoir->dam_volumes[i], 5));
  for (uint i = 0; i < dam_wall_heights.size(); i++) {
    line.push_back(to_string(reservoir->fill_depths[i]));
  }
  if (reservoir->river)
    line.push_back("3");
  else if (reservoir->pit)
    line.push_back("2");
  else if (reservoir->brownfield)
    line.push_back("1");
  else
    line.push_back("0");
  line.push_back(to_string(reservoir->ocean));
  line.push_back(to_string(reservoir->turkey));
  if(RoughGreenfieldReservoir* gr = dynamic_cast<RoughGreenfieldReservoir*>(reservoir))
    for (uint ih = 0; ih < dam_wall_heights.size(); ih++) {
      for (uint idir = 0; idir < directions.size(); idir++) {
        line.push_back(to_string(gr->shape_bound[ih][idir].row));
        line.push_back(to_string(gr->shape_bound[ih][idir].col));
      }
    }
  if(RoughBfieldReservoir* br = dynamic_cast<RoughBfieldReservoir*>(reservoir)){
    line.push_back(to_string(br->shape_bound.size()));
    for(ArrayCoordinate c : br->shape_bound){
      line.push_back(to_string(c.row));
      line.push_back(to_string(c.col));
    }
    if(br->river)
      for(int e : br->elevations)
        line.push_back(to_string(e));
  }
  write_to_csv_file(csv_file, line);
}

vector<unique_ptr<RoughReservoir>> read_rough_reservoir_data(char *filename) {
  vector<unique_ptr<RoughReservoir>> reservoirs;
  ifstream inputFile(filename);
  string s;
  bool header = true;
  bool compressed_format = false;
  while (getline(inputFile, s)) {
    if (header) {
      header = false;
      vector<string> line = read_from_csv_file(s);
      compressed_format = line[6+4*dam_wall_heights.size()+1] == "Ocean";
      continue;
    }
    vector<string> line = read_from_csv_file(s);
    GeographicCoordinate gc =
        GeographicCoordinate_init(stod(line[1]), stod(line[2]));
    GeographicCoordinate origin = get_origin(
        GridSquare_init(convert_to_int(FLOOR(gc.lat)), convert_to_int(FLOOR(gc.lon))),
        border);
    unique_ptr<RoughReservoir> reservoir(new RoughReservoir(convert_coordinates(gc, origin), stoi(line[3])));
    for (uint i = 0; i < dam_wall_heights.size(); i++)
      reservoir->volumes.push_back(stod(line[6 + i]));
    for (uint i = 0; i < dam_wall_heights.size(); i++)
      reservoir->areas.push_back(stod(line[6 + dam_wall_heights.size() + i]));
    for (uint i = 0; i < dam_wall_heights.size(); i++)
      reservoir->dam_volumes.push_back(
          stod(line[6 + 2 * dam_wall_heights.size() + i]));
    for (uint i = 0; i < dam_wall_heights.size(); i++)
      reservoir->fill_depths.push_back(stoi(line[6 + 3 * dam_wall_heights.size() + i]));
    reservoir->max_dam_height = stod(line[4]);
    reservoir->watershed_area = stod(line[5]);
    reservoir->identifier = line[0];
    if(compressed_format){
      reservoir->brownfield =
          stoi(line[6 + 4 * dam_wall_heights.size()]) > 0;
      reservoir->pit =
          stoi(line[6 + 4 * dam_wall_heights.size()]) == 2;
      reservoir->river =
          stoi(line[6 + 4 * dam_wall_heights.size()]) == 3;
      reservoir->ocean =
          stoi(line[6 + 4 * dam_wall_heights.size() + 1]) > 0;
      reservoir->turkey =
          stoi(line[6 + 4 * dam_wall_heights.size() + 2]) > 0;
    }else{
      reservoir->brownfield =
          stoi(line[6 + 4 * dam_wall_heights.size() +
                    (dam_wall_heights.size() * directions.size()) * 2]) > 0;
      reservoir->pit =
          stoi(line[6 + 4 * dam_wall_heights.size() +
                    (dam_wall_heights.size() * directions.size()) * 2]) > 1;
      reservoir->ocean =
          stoi(line[6 + 4 * dam_wall_heights.size() +
                    (dam_wall_heights.size() * directions.size()) * 2 + 1]) > 0;
    }
    for (uint i = 0; i < dam_wall_heights.size(); i++)
      if(reservoir->river)
        reservoir->volumes.push_back(stod(line[6 + i])*60*60*24*365/1e6);
      else
        reservoir->volumes.push_back(stod(line[6 + i]));
    for (uint i = 0; i < dam_wall_heights.size(); i++)
      reservoir->areas.push_back(stod(line[6 + dam_wall_heights.size() + i]));
    for (uint i = 0; i < dam_wall_heights.size(); i++)
      reservoir->dam_volumes.push_back(
          stod(line[6 + 2 * dam_wall_heights.size() + i]));
    reservoir->max_dam_height = stod(line[4]);
    reservoir->watershed_area = stod(line[5]);
    reservoir->identifier = line[0];

    if(!compressed_format || (!reservoir->ocean && !reservoir->brownfield)){
      unique_ptr<RoughGreenfieldReservoir> greenfield_reservoir(new RoughGreenfieldReservoir(*reservoir));
      for (uint ih = 0; ih < dam_wall_heights.size(); ih++) {
        for (uint idir = 0; idir < directions.size(); idir++) {
          greenfield_reservoir->shape_bound[ih][idir].row =
              stoi(line[(compressed_format ? 9 : 6) + 4 * dam_wall_heights.size() +
                        (ih * directions.size() + idir) * 2]);
          greenfield_reservoir->shape_bound[ih][idir].col =
              stoi(line[(compressed_format ? 9 : 6) + 4 * dam_wall_heights.size() + 1 +
                        (ih * directions.size() + idir) * 2]);
          greenfield_reservoir->shape_bound[ih][idir].origin = origin;
        }
      }
      reservoirs.push_back(std::move(greenfield_reservoir));
    } else if(reservoir->ocean || reservoir->brownfield){
      int point_len = stoi(line[9+4*dam_wall_heights.size()]);
      unique_ptr<RoughBfieldReservoir> bfield_reservoir(new RoughBfieldReservoir(*reservoir));
      for(int i = 0; i<point_len; i++){
        bfield_reservoir->shape_bound.push_back(ArrayCoordinate_init(stoi(line[10+4*dam_wall_heights.size()+i*2]), stoi(line[10+4*dam_wall_heights.size()+i*2+1]), origin));
      }
      if(reservoir->river)
        for(int i = 0; i<point_len; i++){
          bfield_reservoir->elevations.push_back(stoi(line[10+3*dam_wall_heights.size()+2*point_len+i]));
        }
      reservoirs.push_back(std::move(bfield_reservoir));
    }
  }
  if (header) {
    cout << "Cannot read empty CSV " << filename << endl;
    throw 1;
  }
  return reservoirs;
}

void write_rough_pair_csv_header(FILE *csv_file) {
  vector<string> header = {"Pair Identifier",
                           "Upper Identifier",
                           "Upper latitude",
                           "Upper longitude",
                           "Upper elevation (m)",
                           "Upper dam height (m)",
                           "Upper max dam height (m)",
                           "Upper water to rock estimate",
                           "Upper area estimate (ha)",
                           "Upper fill depth (m)",
                           "Lower Identifier",
                           "Lower latitude",
                           "Lower longitude",
                           "Lower elevation (m)",
                           "Lower dam height (m)",
                           "Lower max dam height (m)",
                           "Lower water to rock estimate",
                           "Lower area estimate (ha)",
                           "Lower fill depth (m)",
                           "Head (m)",
                           "Pourpoint distance (km)",
                           "Distance (km)",
                           "Slope",
                           "Volume (GL)",
                           "Energy (GWh)",
                           "Storage time (h)",
                           "Figure of merit"};
  write_to_csv_file(csv_file, header);
}

void write_rough_pair_data_header(FILE *csv_file) {
  vector<string> header = {"Pair Identifier",
                           "Upper Identifier",
                           "Upper latitude",
                           "Upper longitude",
                           "Upper elevation (m)",
                           "Upper dam height (m)",
                           "Upper max dam height (m)",
                           "Upper water to rock estimate",
                           "Upper area estimate (ha)",
                           "Upper fill depth (m)",
                           "Upper Brownfield",
                           "Lower Identifier",
                           "Lower latitude",
                           "Lower longitude",
                           "Lower elevation (m)",
                           "Lower dam height (m)",
                           "Lower max dam height (m)",
                           "Lower water to rock estimate",
                           "Lower area estimate (ha)",
                           "Lower fill depth (m)",
                           "Lower Brownfield",
                           "Lower Ocean",
                           "Head (m)",
                           "Pourpoint separation (km)",
                           "Separation (km)",
                           "Slope",
                           "Volume (GL)",
                           "Energy (GWh)",
                           "Storage time (h)",
                           "Figure of merit"};
  write_to_csv_file(csv_file, header);
}

void write_rough_pair_csv(FILE *csv_file, Pair *pair) {
  vector<string> line = {pair->identifier,
                         pair->upper.identifier,
                         dtos(pair->upper.latitude, 4),
                         dtos(pair->upper.longitude, 4),
                         to_string(pair->upper.elevation),
                         dtos(pair->upper.dam_height, 1),
                         dtos(pair->upper.max_dam_height, 0),
                         dtos(pair->upper.water_rock, 1),
                         dtos(pair->upper.area, 1),
                         dtos(pair->upper.fill_depth,1),
                         pair->lower.identifier,
                         dtos(pair->lower.latitude, 4),
                         dtos(pair->lower.longitude, 4),
                         to_string(pair->lower.elevation),
                         dtos(pair->lower.dam_height, 1),
                         dtos(pair->lower.max_dam_height, 0),
                         dtos(pair->lower.water_rock, 1),
                         dtos(pair->lower.area, 1),
                         dtos(pair->lower.fill_depth,1),
                         to_string(pair->head),
                         dtos(pair->pp_distance, 2),
                         dtos(pair->distance, 2),
                         dtos(pair->slope, 2),
                         dtos(pair->required_volume, 2),
                         energy_capacity_to_string(pair->energy_capacity),
                         to_string(pair->storage_time),
                         dtos(pair->FOM, 1)};
  write_to_csv_file(csv_file, line);
}

void write_rough_pair_data(FILE *csv_file, Pair *pair) {
  vector<string> line = {
      pair->identifier,
      pair->upper.identifier,
      dtos(pair->upper.latitude, 6),
      dtos(pair->upper.longitude, 6),
      to_string(pair->upper.elevation),
      dtos(pair->upper.dam_height, 3),
      dtos(pair->upper.max_dam_height, 1),
      dtos(pair->upper.water_rock, 5),
      dtos(pair->upper.area, 1),
      pair->upper.river ? "3" : (pair->upper.pit ? "2" : to_string(pair->upper.brownfield)),
      pair->lower.identifier,
      dtos(pair->lower.latitude, 6),
      dtos(pair->lower.longitude, 6),
      to_string(pair->lower.elevation),
      dtos(pair->lower.dam_height, 3),
      dtos(pair->lower.max_dam_height, 1),
      dtos(pair->lower.water_rock, 5),
      dtos(pair->lower.area, 1),
      pair->lower.river ? "3" : (pair->lower.pit ? "2" : to_string(pair->lower.brownfield)),
      to_string(pair->lower.ocean),
      to_string(pair->head),
      dtos(pair->pp_distance, 5),
      dtos(pair->distance, 5),
      dtos(pair->slope, 6),
      dtos(pair->required_volume, 5),
      energy_capacity_to_string(pair->energy_capacity),
      to_string(pair->storage_time),
      dtos(pair->FOM, 3)};
  write_to_csv_file(csv_file, line);
}

vector<vector<Pair>> read_rough_pair_data(char *filename) {
  vector<vector<Pair>> pairs;
  for (uint i = 0; i < tests.size(); i++) {
    vector<Pair> t;
    pairs.push_back(t);
  }

  ifstream inputFile(filename);
  string s;
  bool header = true;
  while (getline(inputFile, s)) {
    if (header) {
      header = false;
      continue;
    }
    vector<string> line = read_from_csv_file(s);

    Pair pair;
    GeographicCoordinate gc =
        GeographicCoordinate_init(stod(line[2]), stod(line[3]));
    GeographicCoordinate origin =
        get_origin(GridSquare_init(convert_to_int(FLOOR(gc.lat + EPS)),
                                   convert_to_int(FLOOR(gc.lon + EPS))),
                   border);
    pair.upper = Reservoir_init(convert_coordinates(gc, origin), stoi(line[4]));
    gc = GeographicCoordinate_init(stod(line[12]), stod(line[13]));
    origin = get_origin(GridSquare_init(convert_to_int(FLOOR(gc.lat + EPS)),
                                        convert_to_int(FLOOR(gc.lon + EPS))),
                        border);
    pair.lower =
        Reservoir_init(convert_coordinates(gc, origin), stoi(line[14]));

    pair.identifier = line[0];

    pair.upper.identifier = line[1];
    pair.upper.dam_height = stod(line[5]);
    pair.upper.max_dam_height = stod(line[6]);
    pair.upper.water_rock = stod(line[7]);
    pair.upper.area = stod(line[8]);
    pair.upper.fill_depth = stod(line[9]);
    pair.upper.brownfield = stoi(line[10]) > 0;
    pair.upper.pit = stoi(line[10]) == 2;
    pair.upper.river = stoi(line[10]) == 3;

    pair.lower.identifier = line[11];
    pair.lower.dam_height = stod(line[15]);
    pair.lower.max_dam_height = stod(line[16]);
    pair.lower.water_rock = stod(line[17]);
    pair.lower.area = stod(line[18]);
    pair.lower.fill_depth = stod(line[19]);
    pair.lower.brownfield = stoi(line[20]) > 0;
    pair.lower.pit = stoi(line[20]) == 2;
    pair.lower.river = stoi(line[20]) == 3;
    pair.lower.ocean = stoi(line[21]) > 0;

    pair.head = stoi(line[22]);
    pair.pp_distance = stod(line[23]);
    pair.distance = stod(line[24]);
    pair.slope = stod(line[25]);
    pair.required_volume = stod(line[26]);
    pair.upper.volume = stod(line[26]);
    pair.lower.volume = stod(line[26]);
    pair.energy_capacity = stod(line[27]);
    pair.storage_time = stoi(line[28]);

    pair.FOM = stod(line[29]);

    for (uint i = 0; i < tests.size(); i++)
      if (abs(pair.energy_capacity - tests[i].energy_capacity) < EPS &&
          pair.storage_time == tests[i].storage_time)
        pairs[i].push_back(pair);
  }

  if (header)
    throw 1;

  return pairs;
}

void write_pair_csv_header(FILE *csv_file, bool output_FOM) {
  vector<string> header = {"Pair Identifier",
                           "Class",
                           "Head (m)",
                           "Separation (km)",
                           "Slope (%)",
                           "Volume (GL)",
                           "Energy (GWh)",
                           "Storage time (h)",
                           "Combined water to rock ratio",
                           "Country",
                           "Non-overlapping",
                           "Upper Identifier",
                           "Upper elevation (m)",
                           "Upper latitude",
                           "Upper longitude",
                           "Upper reservoir area (ha)",
                           "Upper reservoir volume (GL)",
                           "Upper dam height (m)",
                           "Upper dam length (m)",
                           "Upper dam volume (GL)",
                           "Upper water to rock ratio",
                           "Upper fill depth (m)",
                           "Upper country",
                           "Lower Identifier",
                           "Lower elevation (m)",
                           "Lower latitude",
                           "Lower longitude",
                           "Lower reservoir area (ha)",
                           "Lower reservoir volume (GL)",
                           "Lower dam height (m)",
                           "Lower dam length (m)",
                           "Lower dam volume (GL)",
                           "Lower water to rock ratio",
                           "Lower fill depth (m)",
                           "Lower country"};
  if (output_FOM)
    header.push_back("Figure of Merit");
  write_to_csv_file(csv_file, header);
}

void write_pair_csv(FILE *csv_file, Pair *pair, bool output_FOM) {
  vector<string> line = {
      pair->identifier,
      string(1, pair->category),
      to_string(pair->head),
      dtos(pair->distance, 2),
      dtos(pair->slope * 100, 0),
      dtos(pair->volume, 1),
      energy_capacity_to_string(pair->energy_capacity),
      to_string(pair->storage_time),
      dtos(pair->water_rock, 1),
      pair->country,
      to_string(pair->non_overlap),
      pair->upper.identifier,
      to_string(pair->upper.elevation),
      dtos(pair->upper.latitude, 4),
      dtos(pair->upper.longitude, 4),
      ((pair->upper.brownfield && !pair->upper.pit) ? "NA" : dtos(pair->upper.area, 0)),
      dtos(pair->upper.volume, 1),
      (pair->upper.brownfield
           ? "NA"
           : dtos(pair->upper.dam_height, 1)),
      (pair->upper.brownfield ? "NA" : dtos(pair->upper.dam_length, 0)),
      (pair->upper.brownfield ? "NA" : dtos(pair->upper.dam_volume, 2)),
      (pair->upper.brownfield ? "NA" : dtos(pair->upper.water_rock, 1)),
      dtos(pair->upper.fill_depth,1),
      pair->upper.country,
      pair->lower.identifier,
      to_string(pair->lower.elevation),
      dtos(pair->lower.latitude, 4),
      dtos(pair->lower.longitude, 4),
      (((pair->lower.brownfield && !pair->upper.pit) || pair->lower.ocean)
           ? "NA"
           : dtos(pair->lower.area, 0)),
      dtos(pair->lower.volume, 1),
      ((pair->lower.brownfield || pair->lower.ocean)
           ? "NA"
           : dtos(pair->lower.dam_height, 1)),
      ((pair->lower.brownfield || pair->lower.ocean)
           ? "NA"
           : dtos(pair->lower.dam_length, 0)),
      ((pair->lower.brownfield || pair->lower.ocean)
           ? "NA"
           : dtos(pair->lower.dam_volume, 2)),
      ((pair->lower.brownfield || pair->lower.ocean)
           ? "NA"
           : dtos(pair->lower.water_rock, 1)),
      dtos(pair->lower.fill_depth,1),
      pair->lower.country};
  if (output_FOM)
    line.push_back(dtos(pair->FOM, 0));
  write_to_csv_file(csv_file, line);
}

void write_summary_csv_header(FILE *csv_file) {
  vector<string> header = {"Grid Identifier", "Reservoir type",
                           "Non-overlapping paired sites",
                           "Total paired sites",
                           "Total potential capacity (GWh)"};
  write_to_csv_file(csv_file, header);
}

void write_summary_csv(FILE *csv_file, string square_name, string test,
                      int non_overlapping_sites, int num_sites,
                      int energy_capacity) {
  string str_num_sites = (num_sites == -1) ? "N/A" : to_string(num_sites);
  vector<string> line = {square_name, test, to_string(non_overlapping_sites),
                         str_num_sites, to_string(energy_capacity)};
  write_to_csv_file(csv_file, line);
}

void read_pit_polygons(std::string filename, vector<Pair> &pairs, GridSquare gs){ 
  ifstream inputFile(filename);
  string s;
  bool header = true;
  
  while (getline(inputFile, s)) {
    if (header) {
      header = false;
      vector<string> line = read_from_csv_file(s);
      continue;
    }
    vector<string> line = read_from_csv_file(s);
    GeographicCoordinate origin = get_origin(gs, border);
    
    for(Pair &pair : pairs){
      if (pair.upper.identifier == line[0]){
        pair.upper.shape_bound.clear();
        int point_len = stoi(line[9+4*dam_wall_heights.size()]);
        for(int i = 0; i<point_len; i++){
          ArrayCoordinate array_point = convert_coordinates(convert_coordinates(ArrayCoordinate_init(stoi(line[10+4*dam_wall_heights.size()+i*2]), stoi(line[10+4*dam_wall_heights.size()+i*2+1]), origin)), get_origin(search_config.grid_square, border));
          pair.upper.shape_bound.push_back(array_point);
        }   
      }
      if (pair.lower.identifier == line[0]){
        pair.lower.shape_bound.clear();
        int point_len = stoi(line[9+4*dam_wall_heights.size()]);
        for(int i = 0; i<point_len; i++){
          ArrayCoordinate array_point = convert_coordinates(convert_coordinates(ArrayCoordinate_init(stoi(line[10+4*dam_wall_heights.size()+i*2]), stoi(line[10+4*dam_wall_heights.size()+i*2+1]), origin)), get_origin(search_config.grid_square, border));
          pair.lower.shape_bound.push_back(array_point);
        }   
      }
    }    
  }

  return;
}