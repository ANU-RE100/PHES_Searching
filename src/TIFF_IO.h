#include "model2D.h"


void TIFF_IO_init( );
Model_int32 *TIFF_Read_int32(char *tiff_filename, double *geotransform, char **geoprojection);
Model_int16 *TIFF_Read_int16(char *tiff_filename, double *geotransform, char **geoprojection);
Model_int8 *TIFF_Read_int8(char *tiff_filename, double *geotransform, char **geoprojection);
Model_double *TIFF_Read_double(char *tiff_filename, double *geotransform, char **geoprojection);
void TIFF_Write_int16(char *tiff_filename, double *geotransform, char *geoprojection, Model_int16 *m);
void TIFF_Write_int8(char *tiff_filename, double *geotransform, char *geoprojection, Model_int8 *m);
void TIFF_Write_int32(char *tiff_filename, double *geotransform, char *geoprojection, Model_int32 *m);
void TIFF_Write_double(char *tiff_filename, double *geotransform, char *geoprojection, Model_double *m);