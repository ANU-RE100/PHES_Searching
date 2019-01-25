#ifndef MODEL2D_H
#define MODEL2D_H

// double indexed arrays - ordering model[row][col]
 
struct Model_double {
	int shape[2];
	double **d;
};

struct Model_int32 {
	int shape[2];
	int **d;
};

struct Model_int16 {
	int shape[2];
	int16_t **d;
};

struct Model_int8 {
	int shape[2];
	char **d;
};


#define MODEL_UNSET     0
#define MODEL_SET_ZERO  1

Model_double *Model_double_create(int shape[2], int zero);
Model_int32 *Model_int32_create(int shape[2], int zero);
Model_int16 *Model_int16_create(int shape[2], int zero);
Model_int8 *Model_int8_create(int shape[2], int zero);

void Model_double_free(Model_double *m);
void Model_int32_free(Model_int32 *m);
void Model_int16_free(Model_int16 *m);
void Model_int8_free(Model_int8 *m);

void Model_double_print(Model_double *m);
void Model_int32_print(Model_int32 *m);
void Model_int16_print(Model_int16 *m);
void Model_int8_print(Model_int8 *m);

Model_int16 *Model_double_to_int16(Model_double *m);

#endif
