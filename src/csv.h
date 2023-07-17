#ifndef CSV_H
#define CSV_H

#include "phes_base.h"
#include "mining_pits.h"

void write_to_csv_file(FILE *csv_file, vector<string> cols);
vector<string> read_from_csv_file(string line);
vector<string> read_from_csv_file(string line, char delimeter);

vector<ExistingReservoir> read_existing_reservoir_data(char* filename);
vector<ExistingPit> read_existing_pit_data(char* filename);
vector<string> read_names(char* filename);

void write_rough_reservoir_csv_header(FILE *csv_file);
void write_rough_reservoir_data_header(FILE *csv_file);
void write_rough_reservoir_csv(FILE *csv_file, RoughReservoir *reservoir);
void write_rough_reservoir_data(FILE *csv_file, RoughReservoir *reservoir);
vector<unique_ptr<RoughReservoir>> read_rough_reservoir_data(char *filename);

void write_rough_pair_csv_header(FILE *csv_file);
void write_rough_pair_data_header(FILE *csv_file);
void write_rough_pair_csv(FILE *csv_file, Pair *pair);
void write_rough_pair_data(FILE *csv_file, Pair *pair);
vector<vector<Pair> > read_rough_pair_data(char* filename);

void write_pair_csv_header(FILE *csv_file, bool output_FOM);
void write_pair_csv(FILE *csv_file, Pair *pair, bool output_FOM);
void write_summary_csv_header(FILE *csv_file);
void write_summary_csv(FILE *csv_file, string square_name, string test, 
                      int non_overlapping_sites, int num_sites, 
                      int energy_capacity);
void read_pit_polygons(std::string filename, std::vector<Pair> &pairs, GridSquare gs);

#endif
