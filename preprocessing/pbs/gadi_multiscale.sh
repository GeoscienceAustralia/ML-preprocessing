#PBS -P ge3
#PBS -N multiscale
#PBS -q hugemem
#PBS -l walltime=24:00:00,mem=3068GB,ncpus=96,jobfs=1024GB
#PBS -l storage=scratch/rxh562+gdata/ha3
#PBS -l wd
#PBS -j oe
#PBS -M rakib.hassan@ga.gov.au
#PBS -m bae

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

mpirun -np 32 --map-by ppr:16:node python3 /g/data/ge3/rakib/raijin/soft/geo-wavelets/preprocessing/multiscale.py file_list.txt /g/data/ge3/covariates/national/temp-multiscale 5 --max-search-dist 500 --keep-level 2 --keep-level 3 --keep-level 4 --keep-level 5 --log-level DEBUG 

