#!/bin/env python
"""
Description:
    Computes UPNESS index for a 'filled' DEM in tif format
References:

CreationDate:   03/05/21
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     03/05/21   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import numpy as np
from collections import defaultdict
import sys
import os
import gdal
import click
from osgeo.gdalconst import *
from scipy import interpolate
import ctypes
from ctypes import *
from numpy.ctypeslib import ndpointer
import psutil

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('Filled-DEM-tif', required=True,
                type=click.Path(exists=True))
@click.argument('Output-file-name', required=True,
                type=click.Path())
@click.option('--nodata-source', default=None,
              type=click.Path(exists=True), help='GeoTif file to source nodata values from; default is none. '
                              'Note that the input DEM raster should have nodata values for optimal performance. '
                              'In cases when the input DEM raster does not have nodata values, this secondary raster '
                              'of the same dimensions can be used to source nodata values from')
@click.option('--tolerance', default=0,
              type=float, help='Tolerance for including nodes that are slightly '
                               'lower than a given node being processed; default is 0')
@click.option('--nproc', default=-1,
              type=int, help='Number of processors to use in parallel runs; default is -1 '
                             ', which uses all available processors')
def process(filled_dem_tif, output_file_name, nodata_source, tolerance, nproc, band_idx=1):
    """
    FILLED_DEM_TIF: A filled DEM in GeoTIF format \n
    OUTPUT_FILE_NAME: Name of output GeoTIF file \n
    """

    cpu_count = psutil.cpu_count()
    if(nproc == -1): nproc = cpu_count
    else: nproc = np.min([nproc, cpu_count])

    # load c-module
    path = os.path.dirname(os.path.abspath(__file__))
    lib = None
    try:
        lib = cdll.LoadLibrary(os.path.join(path, "core_c.so"))
    except Exception as e:
        print(e)
        print('Failed to load core_c.so. Make sure the shared libray has been built via the `make` command in %s'%(path))
        exit(0)
    # end try
    compute_upness = lib.compute_upness

    compute_upness.argtypes = [ndpointer(ctypes.c_float),
                               ndpointer(ctypes.c_float),
                               ctypes.c_float,
                               ctypes.c_float,
                               ctypes.c_int,
                               ctypes.c_int,
                               ctypes.c_int]
    
    # load data   
    src_ds = gdal.Open(filled_dem_tif, gdal.GA_ReadOnly)
    band = src_ds.GetRasterBand(band_idx)
    data = np.array(band.ReadAsArray(), dtype=np.float32)
    nodataval = np.float32(band.GetNoDataValue())
    tolerance = np.float32(tolerance)

    if(nodata_source): 
        nodata_array = None
        nodata_ds = gdal.Open(nodata_source, gdal.GA_ReadOnly)
        band = nodata_ds.GetRasterBand(band_idx)
        nodataval = np.float32(band.GetNoDataValue())
        nodata_array = np.array(band.ReadAsArray(), dtype=np.float32)
        
        assert data.shape == nodata_array.shape, 'Dimensions of Input raster and nodata-source must be the same'

        # copy nodata values
        data[nodata_array == nodataval] = nodataval
    # end if
        
    result = np.zeros(data.shape, dtype=np.float32)
    compute_upness(data, result,
                   nodataval,
                   tolerance,
                   data.shape[0],
                   data.shape[1],
                   nproc)

    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.CreateCopy(output_file_name, src_ds, strict=0)
    out_band = out_ds.GetRasterBand(1)

    result[data==nodataval]=nodataval
    out_band.WriteArray(result)
    out_band.ComputeStatistics(False)

    src_ds = None
    out_ds = None
# end func

if __name__ == "__main__":
    process()
# end if
