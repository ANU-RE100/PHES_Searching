#ifndef MODEL2D_H
#define MODEL2D_H

#include <bits/stdc++.h>
using namespace std;

#define MODEL_UNSET     0
#define MODEL_SET_ZERO  1

template<class T>
class Model{
public:
	Model(){
		rows = 0;
		cols = 0;
		data = NULL;
	}
	Model(int rows, int cols){
		this->rows = rows;
		this->cols = cols;
		data = new T[rows*cols];
	}
	Model(int rows, int cols, int zero){
		this->rows = rows;
		this->cols = cols;
		if(zero){
			data = new T[rows*cols]{0};
		}else{
			data = new T[rows*cols];
		}
	}
	~Model(){
		delete[] data;
	}
	int nrows(){
		return rows;
	}
	int ncols(){
		return cols;
	}
	void print(){
		int nx16 = cols>>4;
		int nx32 = cols>>5;
		int ny16 = rows>>4;
		int ny32 = rows>>5;

		printf("       ");
		cout<<"       ";
		for (int i= 0; i<16; i++) {
			cout<<setw(7)<<" "<<nx32 + i*nx16<<" ";
			printf(" %7d ", nx32 + i*nx16);
		}
		printf("\n");
		cout<<"\n";

		for (int j= 0; j<16; j++) {
			int iy = ny32 + j*ny16;
			cout<<setw(4)<<iy<<":  ";
			printf("%4d:  ", iy);
			for (int i= 0; i<16; i++) {
				int ix = nx32 + i*nx16;
				cout<<setw(7)<<" "<<get(iy,ix)<<" ";
				printf(" %7d ", get(iy,ix));
			}
			printf("\n");
			cout<<"\n";
		}
	}
	T get(int row, int col){
		return data[row*cols+col];
	}
	void set(int row, int col, T value){
		data[row*cols+col] = value;
	}
private:
	int rows;
	int cols;
	T *data;
};
 
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
