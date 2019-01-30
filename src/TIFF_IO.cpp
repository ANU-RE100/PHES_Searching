#include <gdal/gdal.h>
#include <gdal/cpl_conv.h>
#include <gdal/cpl_string.h>

#include "phes_base.h"

void TIFF_IO_init( )
{
	GDALAllRegister();
}

Model_int16 *TIFF_Read_int16(char *tiff_filename, double *geotransform, char **geoprojection)
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

	// printf( "Size is %dx%dx%d\n",
	// 	GDALGetRasterXSize( hDataset ),
	// 	GDALGetRasterYSize( hDataset ),
	// 	GDALGetRasterCount( hDataset ) );

//	int nXSize = GDALGetRasterBandXSize( hBand );
//	int nYSize = GDALGetRasterYSize( hDataset );
	int nXSize = GDALGetRasterXSize( hDataset );
	int nYSize = GDALGetRasterYSize( hDataset );

        if( GDALGetProjectionRef( hDataset ) == NULL ) {
		fprintf(stderr, "Cannot get projection from %s\n", tiff_filename);
		throw(1);
	}
	
	//printf( "Projection is `%s'\n", GDALGetProjectionRef( hDataset ) );
	if (geoprojection)
		*geoprojection = strdup(GDALGetProjectionRef( hDataset ));

	if (geotransform) {
		if( GDALGetGeoTransform( hDataset, geotransform ) != CE_None ) {
			fprintf(stderr, "Cannot get transform from %s\n", tiff_filename);
			throw(1);
		}

		// printf( "Origin = (%.6f,%.6f)\n", geotransform[0], geotransform[3] );
		// printf( "Pixel Size = (%.9f,%.9f)\n", geotransform[1], geotransform[5] );
		// printf( "Off diagonal = (%.9f,%.9f)\n", geotransform[2], geotransform[4] );
	}
	
	GDALRasterBandH hBand = GDALGetRasterBand( hDataset, 1 );
//	int   nXSize = GDALGetRasterBandXSize( hBand );
	int nBlockXSize, nBlockYSize;
	GDALGetBlockSize( hBand, &nBlockXSize, &nBlockYSize );

	if (nXSize != nBlockXSize) {
		fprintf(stderr, "Expecting row striped TIFF data - image is %d x %d, blocks are %d x %d\n", nXSize, nYSize, nBlockXSize, nBlockYSize);
		throw(1);
	}

	int shape[2] = { nYSize, nXSize};
	Model_int16 *model = Model_int16_create(shape, MODEL_UNSET);

/* #if double_SIZE == 32 */
/* 	GDALDataType model_type = GDT_double32; */
/* #else */
/* 	GDALDataType model_type = GDT_double64; */
/* #endif */
	for(int y=0; y<nYSize; y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Read, 0, y, nXSize, 1,
					   model->d[y], nXSize, 1,  GDT_Int16, //model_type,
					   0, 0 );
		if (err!=CPLE_None) exit(1);
	}

	return model;
}

Model_int8 *TIFF_Read_int8(char *tiff_filename, double *geotransform, char **geoprojection)
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

	// printf( "Size is %dx%dx%d\n",
	// 	GDALGetRasterXSize( hDataset ),
	// 	GDALGetRasterYSize( hDataset ),
	// 	GDALGetRasterCount( hDataset ) );

//	int nXSize = GDALGetRasterBandXSize( hBand );
//	int nYSize = GDALGetRasterYSize( hDataset );
	int nXSize = GDALGetRasterXSize( hDataset );
	int nYSize = GDALGetRasterYSize( hDataset );

        if( GDALGetProjectionRef( hDataset ) == NULL ) {
		fprintf(stderr, "Cannot get projection from %s\n", tiff_filename);
		throw(1);
	}
	
	//printf( "Projection is `%s'\n", GDALGetProjectionRef( hDataset ) );
	if (geoprojection)
		*geoprojection = strdup(GDALGetProjectionRef( hDataset ));

	if (geotransform) {
		if( GDALGetGeoTransform( hDataset, geotransform ) != CE_None ) {
			fprintf(stderr, "Cannot get transform from %s\n", tiff_filename);
			throw(1);
		}

		// printf( "Origin = (%.6f,%.6f)\n", geotransform[0], geotransform[3] );
		// printf( "Pixel Size = (%.9f,%.9f)\n", geotransform[1], geotransform[5] );
		// printf( "Off diagonal = (%.9f,%.9f)\n", geotransform[2], geotransform[4] );
	}
	
	GDALRasterBandH hBand = GDALGetRasterBand( hDataset, 1 );
//	int   nXSize = GDALGetRasterBandXSize( hBand );
	int nBlockXSize, nBlockYSize;
	GDALGetBlockSize( hBand, &nBlockXSize, &nBlockYSize );

	if (nXSize != nBlockXSize) {
		fprintf(stderr, "Expecting row striped TIFF data - image is %d x %d, blocks are %d x %d\n", nXSize, nYSize, nBlockXSize, nBlockYSize);
		throw(1);
	}

	int shape[2] = { nYSize, nXSize};
	Model_int8 *model = Model_int8_create(shape, MODEL_UNSET);

/* #if double_SIZE == 32 */
/* 	GDALDataType model_type = GDT_double32; */
/* #else */
/* 	GDALDataType model_type = GDT_double64; */
/* #endif */
	for(int y=0; y<nYSize; y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Read, 0, y, nXSize, 1,
					   model->d[y], nXSize, 1,  GDT_Byte, //model_type,
					   0, 0 );
		if (err!=CPLE_None) exit(1);
	}

	return model;
}

Model_int32 *TIFF_Read_int32(char *tiff_filename, double *geotransform, char **geoprojection)
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

	// printf( "Size is %dx%dx%d\n",
	// 	GDALGetRasterXSize( hDataset ),
	// 	GDALGetRasterYSize( hDataset ),
	// 	GDALGetRasterCount( hDataset ) );

//	int nXSize = GDALGetRasterBandXSize( hBand );
//	int nYSize = GDALGetRasterYSize( hDataset );
	int nXSize = GDALGetRasterXSize( hDataset );
	int nYSize = GDALGetRasterYSize( hDataset );

        if( GDALGetProjectionRef( hDataset ) == NULL ) {
		fprintf(stderr, "Cannot get projection from %s\n", tiff_filename);
		throw(1);
	}
	
	//printf( "Projection is `%s'\n", GDALGetProjectionRef( hDataset ) );
	if (geoprojection)
		*geoprojection = strdup(GDALGetProjectionRef( hDataset ));

	if (geotransform) {
		if( GDALGetGeoTransform( hDataset, geotransform ) != CE_None ) {
			fprintf(stderr, "Cannot get transform from %s\n", tiff_filename);
			throw(1);
		}

		// printf( "Origin = (%.6f,%.6f)\n", geotransform[0], geotransform[3] );
		// printf( "Pixel Size = (%.9f,%.9f)\n", geotransform[1], geotransform[5] );
		// printf( "Off diagonal = (%.9f,%.9f)\n", geotransform[2], geotransform[4] );
	}
	
	GDALRasterBandH hBand = GDALGetRasterBand( hDataset, 1 );
//	int   nXSize = GDALGetRasterBandXSize( hBand );
	int nBlockXSize, nBlockYSize;
	GDALGetBlockSize( hBand, &nBlockXSize, &nBlockYSize );

	if (nXSize != nBlockXSize) {
		fprintf(stderr, "Expecting row striped TIFF data - image is %d x %d, blocks are %d x %d\n", nXSize, nYSize, nBlockXSize, nBlockYSize);
		throw(1);
	}

	int shape[2] = { nYSize, nXSize};
	Model_int32 *model = Model_int32_create(shape, MODEL_UNSET);

/* #if double_SIZE == 32 */
/* 	GDALDataType model_type = GDT_double32; */
/* #else */
/* 	GDALDataType model_type = GDT_double64; */
/* #endif */
	for(int y=0; y<nYSize; y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Read, 0, y, nXSize, 1,
					   model->d[y], nXSize, 1,  GDT_Int32, //model_type,
					   0, 0 );
		if (err!=CPLE_None) exit(1);
	}

	return model;
}

Model_double *TIFF_Read_double(char *tiff_filename, double *geotransform, char **geoprojection)
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

	// printf( "Size is %dx%dx%d\n",
	// 	GDALGetRasterXSize( hDataset ),
	// 	GDALGetRasterYSize( hDataset ),
	// 	GDALGetRasterCount( hDataset ) );

//	int nXSize = GDALGetRasterBandXSize( hBand );
//	int nYSize = GDALGetRasterYSize( hDataset );
	int nXSize = GDALGetRasterXSize( hDataset );
	int nYSize = GDALGetRasterYSize( hDataset );

        if( GDALGetProjectionRef( hDataset ) == NULL ) {
		fprintf(stderr, "Cannot get projection from %s\n", tiff_filename);
		throw(1);
	}
	
	//printf( "Projection is `%s'\n", GDALGetProjectionRef( hDataset ) );
	if (geoprojection)
		*geoprojection = strdup(GDALGetProjectionRef( hDataset ));

	if (geotransform) {
		if( GDALGetGeoTransform( hDataset, geotransform ) != CE_None ) {
			fprintf(stderr, "Cannot get transform from %s\n", tiff_filename);
			throw(1);
		}

		// printf( "Origin = (%.6f,%.6f)\n", geotransform[0], geotransform[3] );
		// printf( "Pixel Size = (%.9f,%.9f)\n", geotransform[1], geotransform[5] );
		// printf( "Off diagonal = (%.9f,%.9f)\n", geotransform[2], geotransform[4] );
	}
	
	GDALRasterBandH hBand = GDALGetRasterBand( hDataset, 1 );
//	int   nXSize = GDALGetRasterBandXSize( hBand );
	int nBlockXSize, nBlockYSize;
	GDALGetBlockSize( hBand, &nBlockXSize, &nBlockYSize );

	if (nXSize != nBlockXSize) {
		fprintf(stderr, "Expecting row striped TIFF data - image is %d x %d, blocks are %d x %d\n", nXSize, nYSize, nBlockXSize, nBlockYSize);
		throw(1);
	}

	int shape[2] = { nYSize, nXSize};
	Model_double *model = Model_double_create(shape, MODEL_UNSET);

/* #if double_SIZE == 32 */
/* 	GDALDataType model_type = GDT_double32; */
/* #else */
/* 	GDALDataType model_type = GDT_double64; */
/* #endif */
	for(int y=0; y<nYSize; y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Read, 0, y, nXSize, 1,
					   model->d[y], nXSize, 1,  GDT_Float64, //model_type,
					   0, 0 );
		if (err!=CPLE_None) exit(1);
	}

	return model;
}


void TIFF_Write_int16(char *tiff_filename, double *geotransform, char *geoprojection, Model_int16 *m)
{
	const char *pszFormat = "GTiff";
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	if( hDriver == NULL )
		exit( 1 );

	char **papszOptions = NULL;
	GDALDatasetH OutDS = GDALCreate( hDriver, tiff_filename, m->shape[0], m->shape[1], 1, GDT_Int16,
					 papszOptions );

	GDALSetGeoTransform( OutDS, geotransform );
	GDALSetProjection( OutDS, geoprojection );
	GDALRasterBandH hBand = GDALGetRasterBand( OutDS, 1 );
        for(int y=0; y<m->shape[0]; y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Write, 0, y, m->shape[0], 1,
					   m->d[y], m->shape[0], 1, GDT_Int16, 0, 0 );
                if (err!=CPLE_None) exit(1);
        }

	/* Once we're done, close properly the dataset */
	GDALClose( OutDS );
}


void TIFF_Write_int8(char *tiff_filename, double *geotransform, char *geoprojection, Model_int8 *m)
{
	const char *pszFormat = "GTiff";
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	if( hDriver == NULL )
		exit( 1 );

	char **papszOptions = NULL;
	GDALDatasetH OutDS = GDALCreate( hDriver, tiff_filename, m->shape[0], m->shape[1], 1, GDT_Int16,
					 papszOptions );

	GDALSetGeoTransform( OutDS, geotransform );
	GDALSetProjection( OutDS, geoprojection );

	GDALRasterBandH hBand = GDALGetRasterBand( OutDS, 1 );
        for(int y=0; y<m->shape[0]; y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Write, 0, y, m->shape[0], 1,
					   m->d[y], m->shape[0], 1, GDT_Byte, 0, 0 );
                if (err!=CPLE_None) exit(1);
        }

	/* Once we're done, close properly the dataset */
	GDALClose( OutDS );
}

void TIFF_Write_int32(char *tiff_filename, double *geotransform, char *geoprojection, Model_int32 *m)
{
	const char *pszFormat = "GTiff";
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	if( hDriver == NULL )
		exit( 1 );

	char **papszOptions = NULL;
	GDALDatasetH OutDS = GDALCreate( hDriver, tiff_filename, m->shape[0], m->shape[1], 1, GDT_Int32,
					 papszOptions );

	GDALSetGeoTransform( OutDS, geotransform );
	GDALSetProjection( OutDS, geoprojection );

	GDALRasterBandH hBand = GDALGetRasterBand( OutDS, 1 );
        for(int y=0; y<m->shape[0]; y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Write, 0, y, m->shape[0], 1,
					   m->d[y], m->shape[0], 1, GDT_Int32, 0, 0 );
                if (err!=CPLE_None) exit(1);
        }

	/* Once we're done, close properly the dataset */
	GDALClose( OutDS );
}

void TIFF_Write_double(char *tiff_filename, double *geotransform, char *geoprojection, Model_double *m)
{
	const char *pszFormat = "GTiff";
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	if( hDriver == NULL )
		exit( 1 );

	char **papszOptions = NULL;
	GDALDatasetH OutDS = GDALCreate( hDriver, tiff_filename, m->shape[0], m->shape[1], 1, GDT_Float64,
					 papszOptions );

	GDALSetGeoTransform( OutDS, geotransform );
	GDALSetProjection( OutDS, geoprojection );

	GDALRasterBandH hBand = GDALGetRasterBand( OutDS, 1 );
        for(int y=0; y<m->shape[0]; y++) {
		CPLErr err = GDALRasterIO( hBand, GF_Write, 0, y, m->shape[0], 1,
					   m->d[y], m->shape[0], 1, GDT_Float64, 0, 0 );
                if (err!=CPLE_None) exit(1);
        }

	/* Once we're done, close properly the dataset */
	GDALClose( OutDS );
}
