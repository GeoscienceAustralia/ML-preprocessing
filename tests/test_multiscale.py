import pytest
from preprocessing.multiscale import Multiscale
import numpy as np
from osgeo import gdal
import os
import pywt

import tempfile

@pytest.fixture(params=np.random.randint(100, 2000, 10))
def nx(request):
    return request.param

@pytest.fixture(params=np.random.randint(100, 2000, 10))
def ny(request):
    return request.param

def create_tif(npx, npy, fn):
    image_size = (npx, npy)
    lat_range = [-90, 90]
    lon_range = [-180, 180]

    x = np.linspace(0, 1, npx)
    y = np.linspace(0, 1, npy)
    x, y = np.meshgrid(x, y)
    vals = np.float32(np.sin((x**2+y**2)*np.pi))

    # create 1-band raster file
    dst_ds = gdal.GetDriverByName('GTiff').Create(fn, int(npx), int(npy), 1, gdal.GDT_Float32)
    dst_ds.GetRasterBand(1).WriteArray(vals)
    dst_ds = None
# end func

def test_multiscale(nx, ny):
    path = tempfile.mkdtemp()
    fn = os.path.join(path, 'test.%d.%d.tif'%(nx,ny))
    create_tif(nx, ny, fn)
    
    flist_fn = os.path.join(path, 'flist.txt')
    flist = open(flist_fn, 'w+')
    flist.write(fn)
    flist.close()

    w = pywt.Wavelet('coif6')
    ml = int(np.min(np.array([pywt.dwt_max_level(int(nx), w.dec_len),
                              pywt.dwt_max_level(int(ny), w.dec_len)])))
    print ('Testing Multiscale with nx:%d, ny:%d, max_level:%d'%(nx,ny,ml))
    ms = Multiscale(flist_fn, path, level=ml, file_extension='.tif',
                    mother_wavelet_name='coif6',
                    extension_mode='symmetric',
                    extrapolate=True,
                    max_search_dist=5,
                    smoothing_iterations=10)
    ms.process()

    # os.system() is an expensive call, replace it by os.remove 
    #os.system('rm -f %s' % fn)
    #os.system('rm -f /tmp/flist.txt')
    os.remove(fn)
    os.remove(flist_fn)

    # cleanup output files produced
    for l in range(1,ml+1):
        fn = os.path.join(path, 'test.%d.%d.level_%03d.tif'%(nx,ny, l))
        print("File to be removed ", fn)
        # os.system('rm -f %s' % fn)
        os.remove(fn)
