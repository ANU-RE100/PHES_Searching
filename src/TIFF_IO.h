#ifndef TIFF_IO_H
#define TIFF_IO_H

#include "phes_base.h"

void TIFF_IO_init( );
Model_int32 *TIFF_Read_int32(char *tiff_filename, double *geotransform, char **geoprojection);
Model_int16 *TIFF_Read_int16(char *tiff_filename, double *geotransform, char **geoprojection);
Model_int8 *TIFF_Read_int8(char *tiff_filename, double *geotransform, char **geoprojection);
Model_double *TIFF_Read_double(char *tiff_filename, double *geotransform, char **geoprojection);
void TIFF_Write_int16(char *tiff_filename, double *geotransform, char *geoprojection, Model_int16 *m);
void TIFF_Write_int8(char *tiff_filename, double *geotransform, char *geoprojection, Model_int8 *m);
void TIFF_Write_int32(char *tiff_filename, double *geotransform, char *geoprojection, Model_int32 *m);
void TIFF_Write_double(char *tiff_filename, double *geotransform, char *geoprojection, Model_double *m);

template<class T>
Model<T> TIFF_Read(char *tiff_filename, double *geotransform, char **geoprojection, int data_type)
{
	if(!file_exists(tiff_filename)){
		if(display)
			fprintf(stderr, "No file %s\n", tiff_filename);
		throw(1);
	}
	GDALDatasetH hDataset = GDALOpen(tiff_filename, GA_ReadOnly );
	if( hDataset == NULL ) {
		if(display)
			fprintf(stderr, "Cannot open %s\n", tiff_filename);
		throw(1);
	}

	int nXSize = GDALGetRasterXSize( hDataset );
	int nYSize = GDALGetRasterYSize( hDataset );

	try{
		if (geoprojection)
			*geoprojection = strdup(GDALGetProjectionRef( hDataset ));
	}catch(int e){
		if(display)
			fprintf(stderr, "Cannot get projection from %s\n", tiff_filename);
		throw(1);
	}

	if (geotransform) {
		if( GDALGetGeoTransform( hDataset, geotransform ) != CE_None ) {
			fprintf(stderr, "Cannot get transform from %s\n", tiff_filename);
			throw(1);
		}
	}
	
	GDALRasterBandH hBand = GDALGetRasterBand( hDataset, 1 );
	int nBlockXSize, nBlockYSize;
	GDALGetBlockSize( hBand, &nBlockXSize, &nBlockYSize );

	if (nXSize != nBlockXSize) {
		fprintf(stderr, "Expecting row striped TIFF data - image is %d x %d, blocks are %d x %d\n", nXSize, nYSize, nBlockXSize, nBlockYSize);
		throw(1);
	}

	vector<vector<T> > data;

	for(int row=0; row<nYSize; row++) {
		long temp_arr[nXSize];
		vector<T> temp;
		data.push_back(temp);
		CPLErr err = GDALRasterIO( hBand, GF_Read, 0, row, nXSize, 1,
					   temp_arr, nXSize, 1,  data_type,
					   0, 0 );
		if (err!=CPLE_None) exit(1);
		for(int col = 0; col<nXSize; col++){
			data[row].push_back((T)temp_arr[col]);
		}
	}

	return data;
}

template <class T>
void TIFF_Write(char *tiff_filename, double *geotransform, char *geoprojection, vector<vector<T> >& m, int data_type)
{
	const char *pszFormat = "GTiff";
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	if( hDriver == NULL )
		exit( 1 );

	char **papszOptions = NULL;
	GDALDatasetH OutDS = GDALCreate( hDriver, tiff_filename, m.size(), m[0].size(), 1, data_type,
					 papszOptions );

	GDALSetGeoTransform( OutDS, geotransform );
	GDALSetProjection( OutDS, geoprojection );
	GDALRasterBandH hBand = GDALGetRasterBand( OutDS, 1 );
        for(uint y=0; y<m.size(); y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Write, 0, y, m.size(), 1,
					   m[y].data(), m.size(), 1, data_type, 0, 0 );
        if (err!=CPLE_None) exit(1);
        }

	GDALClose( OutDS );
}

template <class T>
void print(vector<vector<T> >& m)
{
	int nx16 = m[0].size()>>4;
	int nx32 = m[0].size()>>5;
	int ny16 = m.size()>>4;
	int ny32 = m.size()>>5;

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
			printf(" %7d ", (int)m[iy][ix]);
		}
		printf("\n");
	}
}
#endif