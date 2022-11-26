#include "constructor_helpers.hpp"
#include "kml.h"
#include "phes_base.h"
#include <gdal/gdal.h>

string output_kml(Reservoir *reservoir, Reservoir_KML_Coordinates coordinates) {
  string to_return;
  to_return += kml_start;
  to_return += get_reservoir_kml(reservoir, upper_colour, coordinates, NULL);
  to_return += get_dam_kml(reservoir, coordinates);
  to_return += kml_end;
  return to_return;
}

int main(int nargs, char **argv) {
  if (nargs < 5) {
    cout << "Not enough arguements. Need <lon> <lat> <res_id> <dam_height>"
         << endl;
    return -1;
  }
  GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));
  search_config.logger = Logger::DEBUG;

  printf("Reservoir constructor started for %s\n",
         convert_string(str(square_coordinate)));

  GDALAllRegister();
  unsigned long t_usec = walltime_usec();
  parse_variables(convert_string("storage_location"));
  parse_variables(convert_string(file_storage_location + "variables"));

  BigModel big_model = BigModel_init(square_coordinate);
  Model<char> *full_cur_model = new Model<char>(
      big_model.DEM->nrows(), big_model.DEM->ncols(), MODEL_SET_ZERO);
  full_cur_model->set_geodata(big_model.DEM->get_geodata());

  vector<RoughReservoir> reservoirs = read_rough_reservoir_data(
      convert_string(file_storage_location + "processing_files/reservoirs/" +
                     str(square_coordinate) + "_reservoirs_data.csv"));
  search_config.logger.debug("Read in " + to_string(reservoirs.size()) +
                             " reservoirs");

  vector<string> country_names;
  vector<vector<vector<GeographicCoordinate>>> countries = read_countries(
      file_storage_location + "input/countries/countries.txt", country_names);

  string rs(argv[3]);
  for (uint i = 0; i < reservoirs.size(); i++) {
    if (reservoirs[i].identifier == rs) {
      ofstream kml_file(convert_string(file_storage_location + "output/" +
                                       reservoirs[i].identifier + ".kml"),
                        ios::out);

      Reservoir *reservoir = new Reservoir();
      reservoir->identifier = reservoirs[i].identifier;
      reservoir->latitude = reservoirs[i].latitude;
      reservoir->longitude = reservoirs[i].longitude;
      reservoir->elevation = reservoirs[i].elevation;
      reservoir->pour_point = reservoirs[i].pour_point;
      reservoir->volume = -1;
      reservoir->dam_height = atoi(argv[4]);
      Reservoir_KML_Coordinates *coordinates = new Reservoir_KML_Coordinates();

      model_reservoir(reservoir, coordinates, NULL, NULL, NULL, big_model,
                      full_cur_model, countries, country_names);

      kml_file << output_kml(reservoir, *coordinates);
      kml_file.close();
    }
  }

  printf("Reservoir constructor finished for %s. Runtime: %.2f sec\n",
         convert_string(str(square_coordinate)),
         1.0e-6 * (walltime_usec() - t_usec));
}
