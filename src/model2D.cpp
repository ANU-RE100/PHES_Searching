#include "phes_base.h"


Model_double *Model_double_create(int shape[2], int zero)
{
	Model_double *m = (Model_double*)malloc(sizeof(Model_double));
	memcpy(m->shape, shape, 2*sizeof(int));
	double *data_blk;
	if (zero)
		data_blk = (double*)calloc(shape[0]*shape[1], sizeof(double));
	else
		data_blk = (double*)malloc(shape[0]*shape[1]*sizeof(double));
	m->d = (double**)malloc(shape[0]*sizeof(double *));
	for (int iy=0;iy<shape[0];iy++)
		m->d[iy] = data_blk + iy*shape[1];
	return m;
}

Model_int32 *Model_int32_create(int shape[2], int zero)
{
	Model_int32 *m = (Model_int32*)malloc(sizeof(Model_int32));
	memcpy(m->shape, shape, 2*sizeof(int));
	int *data_blk;
	if (zero)
		data_blk = (int*)calloc(shape[0]*shape[1], sizeof(int));
	else
		data_blk = (int*)malloc(shape[0]*shape[1]*sizeof(int));
	m->d = (int**)malloc(shape[0]*sizeof(int *));
	for (int iy=0;iy<shape[0];iy++)
		m->d[iy] = data_blk + iy*shape[1];
	return m;
}

Model_int16 *Model_int16_create(int shape[2], int zero)
{
	Model_int16 *m = (Model_int16*)malloc(sizeof(Model_int16));
	memcpy(m->shape, shape, 2*sizeof(int));
	int16_t *data_blk;
	if (zero)
		data_blk = (int16_t*)calloc(shape[0]*shape[1], sizeof(int16_t));
	else
		data_blk = (int16_t*)malloc(shape[0]*shape[1]*sizeof(int16_t));
	m->d = (int16_t**)malloc(shape[0]*sizeof(int16_t *));
	for (int iy=0;iy<shape[0];iy++)
		m->d[iy] = data_blk + iy*shape[1];
	return m;
}

Model_int8 *Model_int8_create(int shape[2], int zero)
{
	Model_int8 *m = (Model_int8*)malloc(sizeof(Model_int8));
	memcpy(m->shape, shape, 2*sizeof(int));
	char *data_blk;
	if (zero)
		data_blk = (char*)calloc(shape[0]*shape[1], sizeof(char));
	else
		data_blk = (char*)malloc(shape[0]*shape[1]*sizeof(char));
	m->d = (char**)malloc(shape[0]*sizeof(char *));
	for (int iy=0;iy<shape[0];iy++)
		m->d[iy] = data_blk + iy*shape[1];
	return m;
}

void Model_double_free(Model_double *m)
{
	free(&m->d[0][0]);
	free(m->d);
	free(m);		
}

void Model_int32_free(Model_int32 *m)
{
	free(&m->d[0][0]);
	free(m->d);
	free(m);		
}

void Model_int16_free(Model_int16 *m)
{
	free(&m->d[0][0]);
	free(m->d);
	free(m);		
}

void Model_int8_free(Model_int8 *m)
{
	free(&m->d[0][0]);
	free(m->d);
	free(m);		
}

void Model_double_print(Model_double *m)
{
	int nx16 = m->shape[1]>>4;
	int nx32 = m->shape[1]>>5;
	int ny16 = m->shape[0]>>4;
	int ny32 = m->shape[0]>>5;

	printf("       ");
	for (int i= 0; i<16; i++) {
		int ix = nx32 + i*nx16;
		printf(" %8d ", ix);
	}
	printf("\n");

	for (int j= 0; j<16; j++) {
		int iy = ny32 + j*ny16;
		printf("%4d:  ", iy);
		for (int i= 0; i<16; i++) {
			int ix = nx32 + i*nx16;
			printf(" %8.3f ", m->d[iy][ix]);
		}
		printf("\n");
	}
}

void Model_int32_print(Model_int32 *m)
{
	int nx16 = m->shape[1]>>4;
	int nx32 = m->shape[1]>>5;
	int ny16 = m->shape[0]>>4;
	int ny32 = m->shape[0]>>5;

	printf("       ");
	for (int i= 0; i<16; i++) {
		int ix = nx32 + i*nx16;
		printf(" %7d ", ix);
	}
	printf("\n");

	for (int j= 0; j<16; j++) {
		int iy = ny32 + j*ny16;
		printf("%4d:  ", iy);
		for (int i= 0; i<16; i++) {
			int ix = nx32 + i*nx16;
			printf(" %7d ", m->d[iy][ix]);
		}
		printf("\n");
	}
}

void Model_int16_print(Model_int16 *m)
{
	int nx16 = m->shape[1]>>4;
	int nx32 = m->shape[1]>>5;
	int ny16 = m->shape[0]>>4;
	int ny32 = m->shape[0]>>5;

	printf("       ");
	for (int i= 0; i<16; i++) {
		int ix = nx32 + i*nx16;
		printf(" %7d ", ix);
	}
	printf("\n");
	for (int j= 0; j<16; j++) {
		int iy = ny32 + j*ny16;
		printf("%4d:  ", iy);
		for (int i= 0; i<16; i++) {
			int ix = nx32 + i*nx16;
			printf(" %7d ", m->d[iy][ix]);
		}
		printf("\n");
	}
}

void Model_int8_print(Model_int8 *m)
{
	int nx16 = m->shape[1]>>4;
	int nx32 = m->shape[1]>>5;
	int ny16 = m->shape[0]>>4;
	int ny32 = m->shape[0]>>5;

	printf("       ");
	for (int i= 0; i<16; i++) {
		int ix = nx32 + i*nx16;
		printf(" %4d ", ix);
	}
	printf("\n");

	for (int j= 0; j<16; j++) {
		int iy = ny32 + j*ny16;
		printf("%4d:  ", iy);
		for (int i= 0; i<16; i++) {
			int ix = nx32 + i*nx16;
			printf(" %4d ", m->d[iy][ix]);
		}
		printf("\n");
	}
}

Model_int16 *Model_double_to_int16(Model_double *m)
{
	Model_int16 *to_return = Model_int16_create(m->shape, MODEL_UNSET);

	for (int row = 0; row < m->shape[0]; row++)
		for (int col = 0; col <  m->shape[1]; col++) {
			to_return->d[row][col] = convert_to_int(m->d[row][col]);
		}
	return to_return;
}
