#!/bin/env python
"""
Description:
    Replace nans in geotiff files with user-provided nodataval
References:

CreationDate:   20/09/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     20/09/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import numpy as np
from collections import defaultdict
import sys

import gdal
from gdalconst import *

def process(fn, ofn, nodataval):
    src_ds = gdal.Open(fn, gdal.GA_ReadOnly)
    
    driver = gdal.GetDriverByName('GTiff')

    sb = src_ds.GetRasterBand(1)
    od = sb.ReadAsArray()

    od[np.isnan(od)] = nodataval

    of = driver.CreateCopy(ofn, src_ds, strict=0)
    rb = of.GetRasterBand(1)
    rb.WriteArray(od)
    rb.ComputeStatistics(0)
    rb.SetNoDataValue(np.int_(nodataval))
    of = None
# end func

# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    if len(sys.argv) < 4:
        sys.exit('Usage: python %s in.tif out.tif NO_DATA_VALUE' % sys.argv[0])
    # end if
    process(sys.argv[1], sys.argv[2], sys.argv[3])
