#!/bin/env python
"""
Description:
    Class for generating multiscale covariates based on 2D wavelet
    decomposition and reconstruction.
References:

CreationDate:   04/12/17
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     04/12/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os

import numpy as np
from mpi4py import MPI
import glob
from collections import defaultdict
import pywt
import click

try:
    import gdal
    import gdalconst
except ImportError:
    from osgeo import gdal
    from osgeo import gdalconst
# end try

from netCDF4 import Dataset

import logging

log = logging.getLogger('multiscale')

class Coefficients():

    def __init__(self, input_file=None):
        self._input_file = input_file
        self._array = None
        self._level = None
        self._nd = None
        self._indices_dict = None

        if(self._input_file):
            try:
                self._nd = Dataset(self._input_file)
            except Exception as e:
                print('Failed to load file: %s'%(self._input_file))
            # end try

            self._level = self._nd['Parameters'].level

            # load indices
            self._indices_dict = defaultdict(lambda: defaultdict(list))
            count = 1
            for i in np.arange(self._level, 0, -1):
                for k in ['v', 'h', 'd']:
                    self._indices_dict[i][k] = self._nd['indices'][count]

                    # approximation coefficient indices at highest level
                    if(count==1):
                        self._indices_dict[i]['a'] = self._nd['indices'][count-1]
                    # end if

                    count += 1
                # end for
            # end for
        # end if
    # end func

    @staticmethod
    def _save_nc(slices, coeff_array, level, output_file_name):
        root_grp = Dataset(output_file_name, 'w', clobber=True, format='NETCDF4')
        root_grp.description = 'Decomposition coefficients for level %d' % level

        # store level
        pg = root_grp.createGroup('Parameters')
        setattr(pg, 'level', level)
        setattr(pg, 'coeff_nrows', coeff_array.shape[0])
        setattr(pg, 'coeff_ncols', coeff_array.shape[1])

        # store coefficients
        root_grp.createDimension('i', coeff_array.shape[0])
        root_grp.createDimension('j', coeff_array.shape[1])

        coeffs = root_grp.createVariable('coeffs', 'f4', ('i', 'j',))
        coeffs[:, :] = coeff_array

        # store slice indices
        root_grp.createDimension('si', level*3+1)
        root_grp.createDimension('sj', 4) # start and end indices

        indices = root_grp.createVariable('indices', 'i4', ('si', 'sj',))

        idx = 0
        for item in slices:
            if(idx==0):
                row = np.zeros(4)
                for si, ns in enumerate(item):
                    if(ns.start): row[si*2+0] = ns.start
                    if(ns.stop): row[si*2+1] = ns.stop
                # end for

                indices[idx, :] = row
                idx += 1
            else:
                keys=['ad', 'da', 'dd']
                for k in keys:
                    row = np.zeros(4)
                    for si, ns in enumerate(item[k]):
                        if(ns.start): row[si*2+0] = ns.start
                        if(ns.stop): row[si*2+1] = ns.stop
                    # end for

                    indices[idx, :] = row
                    idx += 1
                # end for
            # end if
        # end for
        root_grp.close()
    # end func

    @staticmethod
    def _save_tif(slices, coeff_array, src_ds, level, output_file_basename):
        # extract coefficient indices from slices.
        # note that with respect to the original raster, these indices start from the
        # top-left corner.

        indices = np.int_(np.zeros((level*3+1, 4)))

        idx = 0
        for item in slices:
            if(idx==0):
                row = np.zeros(4)
                for si, ns in enumerate(item):
                    if(ns.start): row[si*2+0] = ns.start
                    if(ns.stop): row[si*2+1] = ns.stop
                # end for

                indices[idx, :] = row
                idx += 1
            else:
                keys=['ad', 'da', 'dd']
                for k in keys:
                    row = np.zeros(4)
                    for si, ns in enumerate(item[k]):
                        if(ns.start): row[si*2+0] = ns.start
                        if(ns.stop): row[si*2+1] = ns.stop
                    # end for

                    indices[idx, :] = row
                    idx += 1
                # end for
            # end if
        # end for

        # get geo-transform from the raster
        gt = np.array(src_ds.GetGeoTransform())

        # output tif files for each component within each decomposition level
        for i in range(indices.shape[0]):
            driver = gdal.GetDriverByName('GTiff')

            outRaster = driver.Create(output_file_basename+'.%03d.tif'%(i),
                                      int(indices[i][3]-indices[i][2]),
                                      int(indices[i][1]-indices[i][0]), 1,
                                      src_ds.GetRasterBand(1).DataType)

            ct = np.array(gt)
            ct[0] += indices[i][2] * ct[1]
            ct[3] += indices[i][0] * ct[5]
            outRaster.SetGeoTransform(tuple(ct))
            outband = outRaster.GetRasterBand(1)
            outband.WriteArray(coeff_array[indices[i][0]:indices[i][1],
                                           indices[i][2]:indices[i][3]])
            outRaster.SetProjection(src_ds.GetProjection())
            outband.FlushCache()

            outband = None
            outRaster = None
            driver = None
        # end for
    # end func

    def get_coefficient_dimensions(self):
        """
        Returns combined coefficient matrix dimensions
        :return: numpy array [nrows, ncols]
        """

        if(self._input_file):
            return np.array([self._nd['Parameters'].coeff_nrows, self._nd['Parameters'].coeff_ncols])
        else:
            return None
        # enf if
    # end func

    def get_coefficients(self, level, component):
        """
        Returns coefficients for a given level and component
        :param level: coefficient level
        :param component: can be either: 'a', 'v', 'h' or 'd', which are the
                          approximation, vertical detail, horizontal detail and
                          diagonal detail coefficients, respectively. Note that
                          approximation coefficients are only available at the
                          highest level computed -- e.g. for a level 4 decomposition,
                          approximation coefficients are only available at that level.
        :return: 1. coefficient array
                 2. starting index (i,j) in the combined coefficient matrix
        """

        if(level < 1 or level > self._level):
            print('Invalid level')
            return None
        # end if

        if(component not in ['a', 'v', 'h', 'd']):
            print('Invalid component')
            return None
        else:
            if(component == 'a' and level < self._level):
                print('Approximation coefficients are only available for level %d'%(self._level))
                return None
            # end if
        # end if

        row_start, row_end, col_start, col_end = self._indices_dict[level][component]
        return self._nd['coeffs'][row_start:row_end, col_start:col_end], np.array([row_start, col_start])
    # end func
# end class

class Multiscale():
    def __init__(self, input, output_folder,
                 level=2, zero='detail', 
                 file_extension='.tif',
                 mother_wavelet_name='coif6',
                 extension_mode='smooth',
                 extrapolate=True,
                 max_search_dist=400,
                 smoothing_iterations=10,
                 keep_level=(), 
                 output_coeffs=False,
                 output_coeffs_format='nc'):
        """
        :param input: a file containing a list of input files (with full path) or a folder containing
                      input files
        :param output_folder: output folder
        :param level: maximum decomposition level
        :param zero: coefficients to set to zero prior to reconstructing
        :param file_extension: file extension e.g. '.tif'
        :param mother_wavelet_name: name of mother wavelet
        :param extension_mode: method of signal extrapolation during computation of wavelet transforms
        :param extrapolate: note, this is separate to the extrapolation done internally by pywavelets,
                            controlled by extention_mode. This parameter controls whether image values
                            are extrapolated into masked regions with NO_DATA_VALUE assigned to them
        :param max_search_dist: this parameter sets the search radius -- in number of pixels -- of
                                extrapolation, controlled by the previous parameter
        :param smoothing_iterations: number of smoothing iterations to be performed after the extrapolation,
                                     controlled by the previous two parameters
        :param keep_level: list of integers that specify the levels to save, while the rest are culled. By
                           default all levels are saved
        :param output_coeffs: output decomposition coefficients
        :param output_coeffs_format: format to output decomposition coefficients in
        """
        self._input = input
        self._output_folder = output_folder
        self._level = level
        self._zero = zero
        self._file_extension = file_extension
        self._mother_wavelet_name = mother_wavelet_name
        self._extension_mode = extension_mode
        self._extrapolate = extrapolate
        self._max_search_dist = max_search_dist
        self._smoothing_iterations = smoothing_iterations
        self._keep_level = keep_level
        self._output_coeffs = output_coeffs
        self._output_coeffs_format = output_coeffs_format

        self._comm = MPI.COMM_WORLD
        self._nproc = self._comm.Get_size()
        self._chunk_index = self._comm.Get_rank()
        self._proc_files = defaultdict(list)

        self.__split_work()
    # end func

    def __get_files(self):
        """
        Function to get a list of input files from a text file or a folder.

        :return: list of files
        """
        files = []
        if(os.path.isdir(self._input)):
            log.info(' Searching for input files with extension %s in folder %s'%
                     (self._file_extension, self._input))
            # prepare case insensitive glob pattern,
            # e.g. for '.pdf', this will produce '*.[Pp][Dd][Ff]'
            if (self._file_extension.count('.') != 1):
                raise RuntimeError('Invalid file extension')

            glob_pattern = '*' + ''.join(sum([['[%s%s]' % (a.upper(), a.lower()) 
                                         if a.isalpha() else a for a in list(x)] for x in self._file_extension], []))

            files = glob.glob(os.path.join(self._input, glob_pattern))
        elif(os.path.isfile(self._input)):
            try:
                fh = open(self._input)
                files_raw = fh.read().splitlines()
                fh.close()

                # filter out commented files
                for f in files_raw:
                    if(len(f) and f.strip()[0] != '#'): files.append(f)
                # end for
            except:
                raise RuntimeError('Failed to read input file')
        log.info(' Found %d files to process ' % len(files))
        return files
    # end func

    def __split_work(self):
        """
        Splits up workload and sends each processor a list of files to process.
        """

        if(self._chunk_index==0):
            files = self.__get_files()

            def split_list(lst, npartitions):
                k, m = divmod(len(lst), npartitions)
                return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
            # end func
            
            self._proc_files = split_list(files, self._nproc)
        # end if

        # broadcast workload to all procs
        if(self._chunk_index==0): log.info(' Distributing workload over %d processors'%(self._nproc))

        self._proc_files = self._comm.bcast(self._proc_files, root=0)
    # end func

    def __generate_reconstructions(self, fname):
        """
        Computes wavelet decompositions and reconstructions.

        :param fname: file name
        """

        # need all data at once
        src_ds = gdal.Open(fname, gdal.GA_ReadOnly)
        od = None
        nodataval = None
        nodataval_mask = None
        if(src_ds.GetRasterBand(1).GetMaskBand() != None):
            driver = gdal.GetDriverByName('GTiff')

            mem_driver = gdal.GetDriverByName('MEM')
            scratch = mem_driver.CreateCopy('', src_ds, strict=0)

            sb = scratch.GetRasterBand(1)
            
            # GetNoDataValue() returns a double, which can be problematic when comparing against
            # pixel values which may be stored as floating point values of lower precision, e.g.
            # float32, float16, etc. We need to cast the 'nodatavalue' to the same format as the
            # pixel values.
            nodataval = sb.GetNoDataValue()
            
            if(nodataval != None):
                od = sb.ReadAsArray()
                
                nodataval = getattr(np, str(od.dtype))(sb.GetNoDataValue()) if nodataval is not None else None
                nodataval_mask = od==nodataval
            # end if

            if(nodataval is not None and self._extrapolate==False):
                log.warning(' NO_DATA_VALUES found in raster %s, but not extrapolating values. This may'%(fname)+\
                            ' cause \'ringing\' artefacts at the edges')
            elif(nodataval is not None and self._extrapolate):
                log.info(' Extrapolating raster %s by %d pixels'%(fname, self._max_search_dist))
                result = gdal.FillNodata(targetBand=sb, maskBand=None,
                                         maxSearchDist=self._max_search_dist,
                                         smoothingIterations=self._smoothing_iterations,
                                         options=['TEMP_FILE_DRIVER=MEM'])

            od = sb.ReadAsArray()
           
            # set NO_DATA_VALUE pixels to the global mean. Note that pywavelets cannot handle
            # masked values
            od[od==nodataval] = np.mean(od[od!=nodataval])

            # clean up
            scratch = None
        else:
            od = src_ds.GetRasterBand(1).ReadAsArray()

        # Compute maximum level computable for raster size
        w = pywt.Wavelet(self._mother_wavelet_name)
        ml = int(np.min(np.array([pywt.dwt_max_level(od.shape[0], w.dec_len),
                                  pywt.dwt_max_level(od.shape[1], w.dec_len)])))

        # generate wavelet decompositions up to required level
        assert(od.ndim==2)
        #print('orig shape:', od.shape)

        # reconstruct each level, starting from the highest
        for l in np.arange(1, self._level+1)[::-1]:
            fn, ext = os.path.splitext(os.path.basename(fname))

            # Culling reconstructed levels based on highest level computable for raster size
            if(l > ml):
                log.warning('Maximum level computable for raster %s is %d; '
                            'skipping level %d'%(os.path.basename(fname), ml, l))
                continue
            # end if

            # Culling reconstructed levels based on user-selection
            if(len(self._keep_level)):
                if(l not in self._keep_level): continue

            log.info('\tReconstructing level: %d'%(l))

            coeffs = pywt.wavedec2(od, self._mother_wavelet_name,
                                   mode=self._extension_mode, level=l)

            if(self._output_coeffs):
                coeff_array, slices = pywt.coeffs_to_array(coeffs)

                if(self._output_coeffs_format == 'nc'):
                    ofn = os.path.join(self._output_folder, '%s.level_%03d.nc' % (fn, l))
                    Coefficients._save_nc(slices, coeff_array, self._level, ofn)
                elif(self._output_coeffs_format == 'tif'):
                    # create output folder
                    os.makedirs(os.path.join(self._output_folder, 'level_%03d_coeffs'%(l)), exist_ok=True)
                    ofbase = os.path.join(self._output_folder, 'level_%03d_coeffs'%(l))
                    ofbase = os.path.join(ofbase, '%s.level_%03d' % (fn, l))
                    
                    Coefficients._save_tif(slices, coeff_array, src_ds, l, ofbase)
                #end if
            # end if

            # zero coefficients
            if(self._zero == 'detail'):
                for i in range(1,len(coeffs)):coeffs[i] = tuple([np.zeros_like(c) for c in coeffs[i]])
            elif(self._zero == 'approx'):
                coeffs[0] = tuple(np.zeros_like(coeffs[0]))
            # end if

            r = pywt.waverec2(coeffs, self._mother_wavelet_name, mode=self._extension_mode)

            p = np.array(r.shape) - np.array(od.shape)
            #print(p, d.shape, od.shape)
            psx = pex = psy = pey = None
            if (p[0] % 2):
                psx = np.floor(p[0] / 2.)
                pex = np.ceil(p[0] / 2.)
            else:
                psx = np.floor(p[0] / 2.)
                pex = np.floor(p[0] / 2.)

            if (p[1] % 2):
                psy = np.floor(p[1] / 2.)
                pey = np.ceil(p[1] / 2.)
            else:
                psy = np.floor(p[1] / 2.)
                pey = np.floor(p[1] / 2.)

            psx, pex, psy, pey = np.int_([psx, pex, psy, pey])
            #print psx,pex,psy,pey

            if(psx != 0 or pex != 0):
                r = r[psx:-pex, :]
            if(psy != 0 or pey != 0):
                r = r[:, psy:-pey]

            if(r.shape != od.shape):
                print ([r.shape, od.shape])
                raise RuntimeError('Error encountered in wavelet reconstruction.')

            ofn = os.path.join(self._output_folder, '%s.level_%03d%s' % (fn, l, ext))
            of = driver.CreateCopy(ofn, src_ds, strict=0)
            rb = of.GetRasterBand(1)

            if(nodataval != None and isinstance(nodataval_mask, np.ndarray)):
                r[nodataval_mask]=nodataval
            # end if

            rb.WriteArray(r)
            rb.ComputeStatistics(False)
            of = None
        # end for

        src_ds = None
    # end func

    def process(self):
        """
        Iterates over a list of files and processes them
        """

        for f in self._proc_files[self._chunk_index]:
            log.info(' Processing %s..'%(f))
            self.__generate_reconstructions(f)
        # end for
    # end func
# end class

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.argument('max-level', required=True,
                type=np.int8)
@click.option('--file-extension', default='.tif',
              help='File extension to use (e.g. \'.tif\') to search for input files; only applicable '
                   'if the \'input\' argument is a folder.',
              type=str)
@click.option('--zero', default='detail', type=click.Choice(['approx', 'detail']), 
              help="Coefficients to set to zero (either 'approx' or 'detail'). Default is 'detail'. "
                   "For lowpass reconstructions 'detail' coefficients should be zeroed out, while "
                   "for highpass reconstructions, 'approx' coefficients should be zeroed out.")
@click.option('--mother-wavelet', default='coif6',
              help='Name of the mother wavelet',
              type=click.Choice(['bior1.1', 'bior1.3', 'bior1.5', 'bior2.2',  'bior2.4',  'bior2.6',  'bior2.8',
                               'bior3.1', 'bior3.3',  'bior3.5',  'bior3.7',  'bior3.9',  'bior4.4',  'bior5.5',
                               'bior6.8',  'coif1',  'coif2',  'coif3',  'coif4',  'coif5',  'coif6',  'coif7',
                               'coif8', 'coif9', 'coif10',  'coif11',  'coif12',  'coif13',  'coif14',  'coif15',
                               'coif16', 'coif17', 'db1',  'db2',  'db3',  'db4',  'db5',  'db6',  'db7',  'db8',
                               'db9',  'db10',  'db11',  'db12',  'db13',  'db14',  'db15',  'db16',  'db17',
                               'db18', 'db19',  'db20',  'db21',  'db22',  'db23',  'db24',  'db25',  'db26',
                               'db27', 'db28', 'db29',  'db30',  'db31',  'db32',  'db33',  'db34', 'db35', 'db36',
                               'db37', 'db38', 'dmey', 'haar', 'rbio1.1', 'rbio1.3', 'rbio1.5', 'rbio2.2', 'rbio2.4',
                               'rbio2.6', 'rbio2.8', 'rbio3.1', 'rbio3.3', 'rbio3.5', 'rbio3.7', 'rbio3.9', 'rbio4.4',
                               'rbio5.5', 'rbio6.8', 'sym2', 'sym3', 'sym4', 'sym5', 'sym6', 'sym7', 'sym8', 'sym9',
                               'sym10', 'sym11', 'sym12', 'sym13', 'sym14', 'sym15', 'sym16', 'sym17', 'sym18',
                               'sym19', 'sym20']))
@click.option('--extension-mode', default='symmetric',
              help="Signal extension mode used for padding",
              type=click.Choice(['zero', 'constant', 'symmetric', 'reflect', 'periodic', 'smooth', 'periodization']))
@click.option('--extrapolate', default=True,
              type=bool,
              help="Extrapolate raster if NO_DATA_VALUES are found. 'Ringing' artefacts can result near sharp contrasts"
                   " in image values -- especially at the edges of NO_DATA_VALUE regions. By extrapolating image values"
                   " to regions of NO_DATA_VALUE, 'ringing' artefacts can be pushed further outward, away from the region"
                   " of interest in the original image. Raster data are extrapolated by default and this parameter has no"
                   " effect when the input raster has no masked"
                   " regions")
@click.option('--max-search-dist', default=500,
              help="Maximum search distance (in pixels) for extrapolating image values to regions of NO_DATA_VALUE in "
                   "input raster; not used if raster has no masked regions")
@click.option('--smoothing-iterations', default=10,
              help="Number of smoothing iterations used for smoothing extrapolated values; see option --max-search-dist")
@click.option('--keep-level', multiple=True, type=(int),
              help="Level to keep. Note that by default all levels up to max-level are saved, which may cause disk"
                   " space issues. This option allows users to save only those levels that are of interest; e.g. to "
                   "keep only levels 5 and 6, this option must be repeated twice for the corresponding levels")
@click.option('--output-coeffs', is_flag=True, help="Output coefficients for each reconstruction level.")
@click.option('--output-coeffs-format', type=click.Choice(['nc', 'tif']),
              default='nc',
              help="The format to output the decomposition coefficients in ('nc' or 'tif'). Default is 'nc'. When "
                   "output in tif format, all constituent coefficients at each decomposition level are output as "
                   "separate files in corresponding folders in the output folder. Note that this parameter has no "
                   "effect when coefficients are not output.")
@click.option('--log-level', default='INFO',
              help="Logging verbosity",
              type=click.Choice(['DEBUG', 'INFO', 'WARN']))
def process(input, output_folder, max_level, file_extension, zero,
            mother_wavelet, extension_mode, extrapolate, max_search_dist,
            smoothing_iterations, keep_level, output_coeffs, output_coeffs_format, log_level):
    """
    INPUT: Path to raster files, or a file containing a list of raster file names (with full path)\n
    OUTPUT_FOLDER: Output folder \n
    MAX_LEVEL: Maximum level up to which wavelet reconstructions are to be computed. Each consecutive
               level halves the raster resolution (doubles the spatial wavelength); thus the spatial
               wavelength at a given level is given by:

               spatial_wavelength = orig_spatial_wavelength * 2^level

    Example usage:
    mpirun -np 2 python multiscale.py filelist.txt /tmp/output 10 --max-search-dist 500

    Running in Cluster Environments:
    This script requires all raster data to be loaded into memory for processing. Furthermore, each
    raster being processed requires 4 times as much memory, e.g. a 7 GB raster will require ~28 GB
    of RAM. Hence, for parallel runs, one must be judicious in terms of allocating processors,
    bearing in mind the amount of memory available on a given compute node on the cluster.
    """

    logMap = {'DEBUG':logging.DEBUG, 'INFO':logging.INFO, 'WARN':logging.WARNING}
    logging.basicConfig(level=logMap[log_level])

    m = Multiscale(input, output_folder, level=max_level, zero=zero, file_extension=file_extension,
                   mother_wavelet_name=mother_wavelet, extension_mode=extension_mode,
                   extrapolate=extrapolate, max_search_dist=max_search_dist,
                   smoothing_iterations=smoothing_iterations, keep_level=keep_level,
                   output_coeffs=output_coeffs, output_coeffs_format=output_coeffs_format)
    m.process()
    return
# end

# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    process()
