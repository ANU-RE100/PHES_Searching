#ifndef MODEL_H
#define MODEL_H

#include "phes_base.h"

#define MODEL_UNSET     0
#define MODEL_SET_ZERO  1

struct Geodata{
	double geotransform[6];
	char* geoprojection;
};

struct GeographicCoordinate{
	double lat,lon;
};

template<class T>
class Model{
public:
	Model(string filename, GDALDataType data_type);
	void write(string filename, GDALDataType data_type);
	void print();
	Model(int rows, int cols, int zero){
		this->rows = rows;
		this->cols = cols;
		if(zero && rows*cols!=0){
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
	T get(int row, int col){
		return data[row*cols+col];
	}
	void set(int row, int col, T value){
		data[row*cols+col] = value;
	}
	Geodata get_geodata(){
		return geodata;
	}
	void set_geodata(Geodata geodata){
		this->geodata = geodata;
	}
	void set_origin(double lat, double lon){
		geodata.geotransform[0] = lon;
		geodata.geotransform[3] = lat;
	}
	bool check_within(int row, int col){
		return (row>=0 && col>=0 && row<nrows() && col<ncols());
	}
	bool check_within(GeographicCoordinate& g){
		return check_within(floor((g.lat-geodata.geotransform[3])/geodata.geotransform[5]), floor((g.lon-geodata.geotransform[0])/geodata.geotransform[1]));
	}
	T get(GeographicCoordinate g){
		return get(floor((g.lat-geodata.geotransform[3])/geodata.geotransform[5]), floor((g.lon-geodata.geotransform[0])/geodata.geotransform[1]));
	}
	void set(GeographicCoordinate& g, T value){
		set(floor((g.lat-geodata.geotransform[3])/geodata.geotransform[5]), floor((g.lon-geodata.geotransform[0])/geodata.geotransform[1]), value);
	}
	GeographicCoordinate get_origin(){
		GeographicCoordinate to_return = {geodata.geotransform[3],geodata.geotransform[0]};
		return to_return;
	}
	GeographicCoordinate get_coordinate(int row, int col){
		GeographicCoordinate to_return = {geodata.geotransform[3]+((double)row+0.5)*geodata.geotransform[5],geodata.geotransform[0]+((double)col+0.5)*geodata.geotransform[1]};
		return to_return;
	}
	vector<GeographicCoordinate> get_corners(){
		vector<GeographicCoordinate> to_return = {
			get_origin(),
			get_coordinate(0,       ncols()),
			get_coordinate(nrows(), ncols()),
			get_coordinate(nrows(), 0      )};
		return to_return;
	}
private:
	int rows;
	int cols;
	T* data;
	Geodata geodata;
};

// Check C++ Code for this
template <typename T> 
Model<T>::Model(string filename, GDALDataType data_type){ 
	char *tif_filename = new char[filename.length() + 1];
	strcpy(tif_filename, filename.c_str());
    if(!file_exists(tif_filename)){
		if(display)
			cout << "No file: "+filename;
		throw(1);
	}
	GDALDataset* Dataset = (GDALDataset *) GDALOpen(tif_filename, GA_ReadOnly);
	if( Dataset == NULL ) {
		if(display)
			cout << "Cannot Open: "+filename;
		throw(1);
	}
	if(Dataset->GetProjectionRef()  != NULL ){
		geodata.geoprojection =  const_cast<char*> (Dataset->GetProjectionRef());
	}else{
		if(display)
			cout << "Cannot get projection from: "+filename;
		throw(1);
	}
	if(Dataset->GetGeoTransform(geodata.geotransform) != CE_None ) {
		if(display)
			cout << "Cannot get transform from: "+filename;
		throw(1);
	}
	GDALRasterBand* Band = Dataset->GetRasterBand( 1 );
	rows = Band->GetYSize();
	cols = Band->GetXSize();
	data = new T[rows*cols];
	T temp_arr[rows];
	for(int row=0; row<rows; row++) {
		CPLErr err = Band->RasterIO( GF_Read, 0, row, cols, 1,
                  temp_arr, cols, 1, data_type,
                  0, 0 );
		if (err!=CPLE_None) exit(1);
		for(int col = 0; col<cols; col++){
			set(row, col, (T)temp_arr[col]);
		}
	}
}

template <typename T> 
void Model<T>::write(string filename, GDALDataType data_type){ 
	char *tif_filename = new char[filename.length() + 1];
	strcpy(tif_filename, filename.c_str());
    const char *pszFormat = "GTiff";
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	if( Driver == NULL )
		exit(1);
	GDALDataset* OutDS = Driver->Create(tif_filename, rows, cols, 1, data_type, NULL);
	OutDS->SetGeoTransform(geodata.geotransform);
	OutDS->SetProjection(geodata.geoprojection);
	GDALRasterBand* Band = OutDS->GetRasterBand(1);
	T temp_arr[rows];
    for(int row=0; row<rows; row++) {
    	for(int col = 0; col<cols; col++){
			temp_arr[col] = get(row, col);
		}
		CPLErr err = Band->RasterIO(GF_Write, 0, row, rows, 1,
					   temp_arr, rows, 1, data_type, 0, 0 );
        if (err!=CPLE_None) exit(1);
    }
	GDALClose((GDALDatasetH) OutDS );
}

template <typename T>
void Model<T>::print(){
	int nx16 = cols>>4;
	int nx32 = cols>>5;
	int ny16 = rows>>4;
	int ny32 = rows>>5;
	cout<<"       ";
	for (int i= 0; i<16; i++) {
		cout<<" "<<setw(7)<<nx32 + i*nx16<<" ";
	}
	cout<<"\n";

	for (int j= 0; j<16; j++) {
		int iy = ny32 + j*ny16;
		cout<<setw(4)<<iy<<":  ";
		for (int i= 0; i<16; i++) {
			int ix = nx32 + i*nx16;
			cout<<" "<<setw(7)<<+get(iy,ix)<<" ";
		}
		cout<<"\n";
	}
}



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
