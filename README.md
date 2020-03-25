# geo-wavelets
2D Wavelet decomposition/reconstruction for raster data

## Dependencies:

Install the following after loading required modules (listed in the next section) on Gadi:
* `pywavelts`:
  
  Install as: `pip3.7 install pywavelets --user`

There are no major additional dependencies other than the ones defined below, which are already available on
Gadi.

## Running on Gadi:

Load modules and define paths (amend as required) as follows for interactive runs:

```
module purge
module load python3/3.7.4
module load gdal/3.0.2
module load openmpi/2.1.6-mt
module load hdf5/1.10.5p
export GDAL_DATA=/apps/gdal/3.0.2/share/gdal/
export PYTHONPATH=/apps/gdal/3.0.2/lib64:/apps/gdal/3.0.2/lib64/python3.7/site-packages
export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=/g/data/ge3/rakib/raijin/soft/geo-wavelets:$PYTHONPATH
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8
```

For `pbs` batch jobs, check example pbs scripts in the `preprocessing/pbs` folder in this repository.

## Running Tests on Gadi:

Note that if you have a local version of Python3 installed in your environment, running `pytest` in the `tests` folder will likely fail, with `gdal` being reported as missing. In such cases, run `pytest` as follows:

`python3 -m pytest`

The above will run the system version of `pytest` and the tests should all runs to completion.

## User Parameters 

`multiscale.py` takes the following input parameters:

```
Usage: multiscale.py [OPTIONS] INPUT OUTPUT_FOLDER MAX_LEVEL

  INPUT: Path to raster files, or a file containing a list of raster file
  names (with full path)

  OUTPUT_FOLDER: Output folder

  MAX_LEVEL: Maximum level up to which wavelet reconstructions are to be
  computed. Each consecutive            level halves the raster resolution
  (doubles the spatial wavelength); thus the spatial            wavelength
  at a given level is given by:

             spatial_wavelength = orig_spatial_wavelength * 2^level

  Example usage: mpirun -np 2 python multiscale.py filelist.txt /tmp/output
  10 --max-search-dist 500

  Running in Cluster Environments: This script requires all raster data to
  be loaded into memory for processing. Furthermore, each raster being
  processed requires 4 times as much memory, e.g. a 7 GB raster will require
  ~28 GB of RAM. Hence, for parallel runs, one must be judicious in terms of
  allocating processors, bearing in mind the amount of memory available on a
  given compute node on the cluster.

Options:
  --file-extension TEXT           File extension to use (e.g. '.tif') to
                                  search for input files; only applicableif
                                  the 'input' argument is a folder.

  --mother-wavelet [bior1.1|bior1.3|bior1.5|bior2.2|bior2.4|bior2.6|bior2.8|bior3.1|bior3.3|bior3.5|bior3.7|bior3.9|bior4.4|bior5.5|bior6.8|coif1|coif2|coif3|coif4|coif5|coif6|coif7|coif8|coif9|coif10|coif11|coif12|coif13|coif14|coif15|coif16|coif17|db1|db2|db3|db4|db5|db6|db7|db8|db9|db10|db11|db12|db13|db14|db15|db16|db17|db18|db19|db20|db21|db22|db23|db24|db25|db26|db27|db28|db29|db30|db31|db32|db33|db34|db35|db36|db37|db38|dmey|haar|rbio1.1|rbio1.3|rbio1.5|rbio2.2|rbio2.4|rbio2.6|rbio2.8|rbio3.1|rbio3.3|rbio3.5|rbio3.7|rbio3.9|rbio4.4|rbio5.5|rbio6.8|sym2|sym3|sym4|sym5|sym6|sym7|sym8|sym9|sym10|sym11|sym12|sym13|sym14|sym15|sym16|sym17|sym18|sym19|sym20]
                                  Name of the mother wavelet
  --extension-mode [zero|constant|symmetric|reflect|periodic|smooth|periodization]
                                  Signal extension mode used for padding
  --extrapolate BOOLEAN           Extrapolate raster if NO_DATA_VALUES are
                                  found. 'Ringing' artefacts can result near
                                  sharp contrasts in image values --
                                  especially at the edges of NO_DATA_VALUE
                                  regions. By extrapolating image values to
                                  regions of NO_DATA_VALUE, 'ringing'
                                  artefacts can be pushed further outward,
                                  away from the region of interest in the
                                  original image. Raster data are extrapolated
                                  by default and this parameter has no effect
                                  when the input raster has no masked regions

  --max-search-dist INTEGER       Maximum search distance (in pixels) for
                                  extrapolating image values to regions of
                                  NO_DATA_VALUE in input raster; not used if
                                  raster has no masked regions

  --smoothing-iterations INTEGER  Number of smoothing iterations used for
                                  smoothing extrapolated values; see option
                                  --max-search-dist

  --keep-level INTEGER            Level to keep. Note that by default all
                                  levels up to max-level are saved, which may
                                  cause disk space issues. This option allows
                                  users to save only those levels that are of
                                  interest; e.g. to keep only levels 5 and 6,
                                  this option must be repeated twice for the
                                  corresponding levels

  --log-level [DEBUG|INFO|WARN]   Logging verbosity
  -h, --help                      Show this message and exit.

```
