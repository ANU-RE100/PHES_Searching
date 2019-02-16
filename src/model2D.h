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
	Model(int rows, int cols){
		this->rows = rows;
		this->cols = cols;
		data = new T[rows*cols];
	}
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
	T* get_pointer(int row, int col){
		return &data[row*cols+col];
	}
	void set(int row, int col, T value){
		data[row*cols+col] = value;
	}
	void set(T value){
		for(int row=0; row<rows; row++)
			for(int col=0; col<cols; col++)
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

template <typename T> 
Model<T>::Model(string filename, GDALDataType data_type){ 
	char *tif_filename = new char[filename.length() + 1];
	strcpy(tif_filename, filename.c_str());
    if(!file_exists(tif_filename)){
		if(display)
			cout << "No file: "+filename+"\n";
		throw(1);
	}
	GDALDataset* Dataset = (GDALDataset *) GDALOpen(tif_filename, GA_ReadOnly);
	if( Dataset == NULL ) {
		if(display)
			cout << "Cannot Open: "+filename+"\n";
		throw(1);
	}
	if(Dataset->GetProjectionRef()  != NULL ){
		geodata.geoprojection =  const_cast<char*> (Dataset->GetProjectionRef());
	}else{
		if(display)
			cout << "Cannot get projection from: "+filename+"\n";
		throw(1);
	}
	if(Dataset->GetGeoTransform(geodata.geotransform) != CE_None ) {
		if(display)
			cout << "Cannot get transform from: "+filename+"\n";
		throw(1);
	}
	GDALRasterBand* Band = Dataset->GetRasterBand( 1 );
	rows = Band->GetYSize();
	cols = Band->GetXSize();
	int temp_cols = cols;
	if(cols==1801){
		cols = 3601;
		geodata.geotransform[1] = geodata.geotransform[1]/2.0;
	}
	data = new T[rows*cols];
	T temp_arr[temp_cols];

	for(int row=0; row<rows; row++) {
		CPLErr err = Band->RasterIO( GF_Read, 0, row, temp_cols, 1,
                  temp_arr, temp_cols, 1, data_type,
                  0, 0 );
		if (err!=CPLE_None) exit(1);
		if (temp_cols==1801){
			for(int col = 0; col<temp_cols-1; col++){
				set(row, col*2, (T)temp_arr[col]);
				set(row, col*2+1, (T)temp_arr[col]);
			}
			set(row, 3600, (T)temp_arr[1800]);
		}else{
			for(int col = 0; col<temp_cols; col++)
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
		cout<<" "<<setw(8)<<nx32 + i*nx16<<" ";
	}
	cout<<"\n";

	for (int j= 0; j<16; j++) {
		int iy = ny32 + j*ny16;
		cout<<setw(4)<<iy<<":  ";
		for (int i= 0; i<16; i++) {
			int ix = nx32 + i*nx16;
			cout<<" "<<setw(8)<<+get(iy,ix)<<" ";
		}
		cout<<"\n";
	}
}

#endif
