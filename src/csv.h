#ifndef CSV_H
#define CSV_H

#include "phes_base.h"

void write_to_csv_file(FILE *csv_file, vector<string> cols);
vector<string> read_from_csv_file(string line);

void write_rough_reservoir_csv_header(FILE *csv_file);
void write_rough_reservoir_data_header(FILE *csv_file);
void write_rough_reservoir_csv(FILE *csv_file, RoughReservoir reservoir);
void write_rough_reservoir_data(FILE *csv_file, RoughReservoir reservoir);
vector<RoughReservoir> read_rough_reservoir_data(char* filename);

void write_rough_pair_csv_header(FILE *csv_file);
void write_rough_pair_data_header(FILE *csv_file);
void write_rough_pair_csv(FILE *csv_file, Pair *pair);
void write_rough_pair_data(FILE *csv_file, Pair *pair);
vector<vector<Pair> > read_rough_pair_data(char* filename);

void write_pair_csv_header(FILE *csv_file);
void write_pair_csv(FILE *csv_file, Pair *pair);
void write_total_csv_header(FILE *csv_file);
void write_total_csv(FILE *csv_file, string square_name, int num_sites, int energy_capacity);

#endif
