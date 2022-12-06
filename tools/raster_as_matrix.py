'''
Used to view the .tif DEM rasters as a matrix of values when debugging. Also allows the matrix to be written to a .csv file.
'''

import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt


def printRaster(DEM_filepath, output_filename, model_type, write_bool):
    ds = gdal.Open(DEM_filepath)

    DEM_array = np.array(ds.GetRasterBand(1).ReadAsArray())


    print("Shape: " + str(DEM_array.shape))
    print("Size: " + str(DEM_array.size))
    print("Array: \n")
    print(DEM_array)

    if model_type == bool:
        DEM_bool_scatter = [[],[]]
        for x in range(0,len(DEM_array[0])):
            for y in range(0,len(DEM_array)):
                if DEM_array[y,x] == 1:
                    DEM_bool_scatter[0].append(x)
                    DEM_bool_scatter[1].append(y)
        plt.plot(DEM_bool_scatter[0], DEM_bool_scatter[1], 's')
    
    plt.show()          

    if write_bool == True:
        np.savetxt("tools/DEM_arrays/" + output_filename, DEM_array.astype(int), fmt='%i', delimiter=",")

    return

if __name__ == '__main__':
    printRaster("debug/TN_scanning_region/s18_e144_TN_scanning_region.tif", "test.csv", bool, True)

