#include "phes_base.h"
#include "kml.h"

int display = false;

int main(int nargs, char **argv)
{

	GridSquare square_coordinate = GridSquare_init(atoi(argv[2]), atoi(argv[1]));
    if(nargs>3)
        display = atoi(argv[3]);

    printf("Convert to shp started for %s\n",convert_string(str(square_coordinate)));

    GDALAllRegister();
    parse_variables(convert_string("storage_location"));
    parse_variables(convert_string(file_storage_location+"variables"));
    unsigned long t_usec = walltime_usec();

    

	// pairs = read_rough_pair_data(convert_string(file_storage_location+"processing_files/pretty_set_pairs/"+str(square_coordinate)+"_rough_pretty_set_pairs_data.csv"));
	// mkdir(convert_string(file_storage_location+"output/final_output"), 0777);
	// mkdir(convert_string(file_storage_location+"output/final_output/"+str(square_coordinate)),0777);

	// FILE *total_csv_file = fopen(convert_string(file_storage_location+"output/final_output/"+str(square_coordinate)+"/"+str(square_coordinate)+"_total.csv"), "w");
	// write_total_csv_header(total_csv_file);

	// int total_count = 0;
	// int total_capacity = 0;
	// for(uint i = 0; i<tests.size(); i++){
 //        sort(pairs[i].begin(), pairs[i].end());
        
	// 	FILE *csv_file = fopen(convert_string(file_storage_location+"output/final_output/"+str(square_coordinate)+"/"+str(square_coordinate)+"_"+str(tests[i])+".csv"), "w");
	// 	write_pair_csv_header(csv_file);

	// 	ofstream kml_file(convert_string(file_storage_location+"output/final_output/"+str(square_coordinate)+"/"+str(square_coordinate)+"_"+str(tests[i])+".kml"), ios::out);
	// 	KML_Holder kml_holder;

	// 	//FILE *fusion_csv_file = fopen(convert_string(file_storage_location+"Output/Final Output/"+str(square_coordinate)+"/"+str(square_coordinate)+"_"+str(tests[i])+"_Fusion_Table.csv"), "w");
	// 	//write_fusion_csv_header(fusion_csv_file);

	// 	sort(pairs[i].begin(), pairs[i].end());
	// 	int count = 0;
	// 	for(uint j=0; j<pairs[i].size(); j++){
	// 		Pair_KML pair_kml;
	// 		bool non_overlap;
	// 		if(model_pair(&pairs[i][j], &pair_kml, seen, &non_overlap, tests[i].max_FOM, big_model, full_cur_model)){
	// 			write_pair_csv(csv_file, &pairs[i][j]);
	// 			//write_fusion_csv(fusion_csv_file, &pairs[i][j], &pair_kml);
	// 			update_kml_holder(&kml_holder, &pairs[i][j], &pair_kml);
	// 			count++;
	// 			if(non_overlap){
	// 				total_count++;
	// 				total_capacity+=tests[i].energy_capacity;
	// 			}
	// 		}
	// 	}

	// 	kml_file << output_kml(&kml_holder, str(square_coordinate), tests[i]);
 //        if(display)
 //            printf("%d %dGWh %dh Pairs\n", count, tests[i].energy_capacity, tests[i].storage_time);
	// 	kml_file.close();
 //        fclose(csv_file);
	// }
	// write_total_csv(total_csv_file, str(square_coordinate), total_count, total_capacity);
    printf("Convert to shp finished for %s. Runtime: %.2f sec\n", convert_string(str(square_coordinate)), 1.0e-6*(walltime_usec() - t_usec) );
}
