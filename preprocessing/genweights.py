#!/bin/env python
"""
Description:
    Generate geochemical weights
References:

CreationDate:   23/05/21
Developer:      rakib.hassan@ga.com.au

Revision History:
    LastUpdate:     23/05/21   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import numpy as np
from collections import defaultdict
import sys
import gdal
import osr
import pyproj
import shapefile
from shapely.geometry.polygon import Polygon, Point
from tqdm import tqdm
import datetime
import geopandas
from pandas import read_csv
from rtree import index
from scipy.signal import tukey
from scipy.interpolate import interp1d
from scipy.stats import rankdata
from mpi4py import MPI
from affine import Affine
import click
from scipy.interpolate import interp1d

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func

def gen_colname(stem, ext, delim='_'):
    mapping = {'FE2O3TOT':'FE2O3T',
               'F_99999':'F',
               'H2OPLUS':'H2O_P',
               'H2OMINUS':'H2O_M'}

    if(stem in mapping.keys()): stem = mapping[stem]

    colname = '%s%s%s'%(stem, delim, ext)
    assert(len(colname) <= 10)

    #if(len(colname)>8): print([colname, len(colname)])
    return colname
# end func

def process_species(rank, sname, vals, mapsymbols, in_ocean):
    def weight_function(n):
        x = np.linspace(0, 1, n)
        a=0.15
        b=.07
        wf = (0.5*np.tanh((x-a)/b) + 0.5*np.tanh((x[::-1]-a)/b))
        wf[np.argwhere(wf==np.max(wf))[0]] = 1.

        return x, wf
    # end func

    colname_wp = gen_colname(sname, 'Wp')
    colname_whs = gen_colname(sname, 'Whs')
    
    result_wp = np.ones(len(vals))*-99999
    result_whs = np.ones(len(vals))*-99999
    
    for mapsym in tqdm(set(mapsymbols), desc='rank: %d'%rank):
        indices = np.argwhere(mapsymbols == mapsym)[:,0]
        valid_indices = indices[(vals[indices] != -99999) & \
                                (~in_ocean[indices])]
        
        valid_vals = vals[valid_indices]
        
        #print '%s, %s -> %d'%(mapsym, sname, len(valid_vals))
        
        if(len(valid_vals)>3):
            NP = len(valid_vals)

            x, wf = weight_function(NP)
            #print np.max(wf), np.argwhere(wf==np.max(wf))

            rank = rankdata(valid_vals, method='ordinal')-1
            
            # multiple instances of the median value should also have a weight of 1.
            mval = valid_vals[wf[rank] == 1.]
            wf[rank[valid_vals==mval]] = 1.
            
            wperc = wf[rank]
            result_wp[valid_indices] = wperc
        elif(len(valid_vals) == 3):
            # For 3 sites within a mapsym == median gets highest weight; the lower values get 25th percentile 
            # and the highest gets the 75th percentile weight
            mean = np.mean(valid_vals)
            order = np.argsort(valid_vals)
            x, wf = weight_function(3)

            # note we use a quadratic polynomial fit to extract the 25th and 75th percentile values
            pf = np.polyfit(x, wf, 2)
            p = np.poly1d(pf)

            result_wp[valid_indices[order[0]]] = p(0.25)
            result_wp[valid_indices[order[1]]] = 1.
            result_wp[valid_indices[order[2]]] = p(0.75)
        elif(len(valid_vals) == 2):
            # For 2 sites within a mapsym == if they are within 50% of the average they get full weight; 
            # else they get the 25th or 75th percentile weight
            mean = np.mean(valid_vals)
            order = np.argsort(valid_vals)
            
            if((np.fabs(valid_vals[0] - mean) < mean*0.5) and (np.fabs(valid_vals[1] - mean) < mean*0.5)):
                result_wp[valid_indices] = 1.
            else:
                x, wf = weight_function(3)
                # note we use a quadratic polynomial fit to extract the 25th and 75th percentile values
                pf = np.polyfit(x, wf, 2)
                p = np.poly1d(pf)

                result_wp[valid_indices[order[0]]] = p(0.25)
                result_wp[valid_indices[order[1]]] = p(0.75)
            # end if
        elif(len(valid_vals) == 1):
            # For 1 site within a mapsym == get highest weight
            result_wp[valid_indices] = 1.
        # end if

        # compute whs
        if(len(valid_vals)):
            p25 = np.percentile(valid_vals, 25)
            p50 = np.percentile(valid_vals, 50)
            p75 = np.percentile(valid_vals, 75)
            iqr = p75 - p25

            result_whs[valid_indices[(valid_vals < (p25 - 3*iqr))]] = 0.1 # far lower outliers
            result_whs[valid_indices[(valid_vals > (p75 + 3*iqr))]] = 0.1 # far upper outliers

            result_whs[valid_indices[(valid_vals >= (p25 - 3*iqr)) & (valid_vals <= (p25 - 1.5*iqr))]] = 0.25 # lower outliers
            result_whs[valid_indices[(valid_vals >= (p75 + 1.5*iqr)) & (valid_vals <= (p75 + 3*iqr))]] = 0.25 # upper outliers
            
            result_whs[valid_indices[(valid_vals >= (p25 - 1.5*iqr)) & (valid_vals <= p25)]] = 0.75 # lower peripherals 
            result_whs[valid_indices[(valid_vals >= p75) & (valid_vals <= (p75 + 1.5*iqr))]] = 0.75 # upper peripherals

            result_whs[valid_indices[(valid_vals >= p25) & (valid_vals <= p75)]] = 1 # central
        # end if
    # end for
    return {colname_wp:result_wp, colname_whs:result_whs}
# end func

def process_collocated(df, mapsym_raster_fn, mapsym_grid_fn, species):
    mapsym_csv = read_csv(mapsym_grid_fn)
    mapsym_lt = defaultdict(list)

    # create value->mapsym reverse-lookup-table
    for i in np.arange(len(mapsym_csv)):
        mapsym_lt[mapsym_csv['Value'][i]] = mapsym_csv['MAPSYMBOL'][i]
    # end for

    mapsym_r = gdal.Open(mapsym_raster_fn, gdal.GA_ReadOnly)
    band = mapsym_r.GetRasterBand(1)

    # Convert geographic co-ordinates to pixel co-ordinates
    affine_forward_transform = Affine.from_gdal(*mapsym_r.GetGeoTransform())
    affine_reverse_transform = ~(affine_forward_transform)

    # find collocated sites
    collision_map = defaultdict(list)
    mapsym_match = np.ones(len(df))
    collocated_px = np.ones(len(df))*-99999
    collocated_py = np.ones(len(df))*-99999
    for i in np.arange(len(df)):
        tx, ty = df.geometry[i].x, df.geometry[i].y

        # Convert geographic co-ordinates to pixel co-ordinates
        px, py = affine_reverse_transform * (tx, ty)
        px, py = int(px + 0.), int(py + 0.)
        
        pixel_id = py * mapsym_r.RasterXSize + px
        
        value = 0
        try:
            value = band.ReadAsArray(px, py, 1, 1)[0,0]
        except:
            print ('Site (lon:%f, lat:%f) is outside the mapsym raster..'%(df['LONGITUDE'][i], df['LATITUDE'][i]))
            continue
        # end try 
        
        mapsym_match[i] = (mapsym_lt[int(value)] == df['MAPSYMBOL'][i])
        collision_map[(px,py)].append((i, pixel_id)) # append index in df and pixel id
    # end for
    
    for k,v in collision_map.iteritems():
        collision_map[k] = np.array(v)
    # end for

    #populate collocated arrays
    for k,v in collision_map.iteritems():
        if(len(v) > 1): 
            for i in np.arange(len(v)):
                collocated_px[v[i, 0]] = k[0]
                collocated_py[v[i, 0]] = k[1]
            # end for
        # end if
    # end for

    # create collocated-wperc columns for all species and assign max(wperc) to each collocated entry.
    # note that collocated entries where raster mapsys does not match df[MAPSYM] are ingored.
    for s in tqdm(species, desc='Processing collocated sites'):
        wpc_array = np.zeros(len(df))
        wpc_colname = gen_colname(s, 'Wpc') # new collocated wperc for current specie
        wp_colname = gen_colname(s, 'Wp') # wperc for current specie (already created in stage one)
        
        for k,v in collision_map.iteritems():
            if(len(v) > 1):

                indices = []
                for i in np.arange(len(v)):
                    idx = v[i, 0]
                    if(mapsym_match[idx]):
                        indices.append(idx)
                    # end if
                # end for
                indices = np.array(indices)
                
                if(len(indices)):
                    wpc_array[indices] = np.max(df[wp_colname][indices])
                # end if

            # end if
        # end for

        df[wpc_colname] = wpc_array
    # end for
    
    df['MAPSYM_EQ'] = mapsym_match
    df['COLLOC_PX'] = collocated_px
    df['COLLOC_PY'] = collocated_py
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('master-shapefile', required=True,
                type=click.Path(exists=True))
@click.argument('geology-lines-buffer-shapefile', required=True,
                type=click.Path(exists=True))
@click.argument('mapsym-raster', required=True,
                type=click.Path(exists=True))
@click.argument('mapsym-grid', required=True,
                type=click.Path(exists=True))
@click.argument('output-file-basename', required=True,
                type=click.Path())
def process(master_shapefile, geology_lines_buffer_shapefile, mapsym_raster, mapsym_grid, output_file_basename):
    """
    MASTER_SHAPEFILE: master shapefile containing all the geochemical species

    GEOLOGY_LINES_BUFFER_SHAPEFILE: shapefile with buffers around line geological boundaries

    MAPSYM_RASTER: A rasterized version of MAPSYMBOL in geoTIF format

    MAPSYM_GRID: A text file containing the Raster Attribute Table of MAPSYM_RASTER

    OUTPUT_FILE_BASENAME: output file basename
    """

    species = ['SIO2', 'TIO2', 'AL2O3', 'FE2O3TOT', 'FE2O3', 'FEO', 'MNO', 'MGO', 'CAO', 'NA2O', 'K2O', 'P2O5', 'H2OPLUS', 'H2OMINUS', 'CO2', 'SO3', 'MLOI', 'LOITOT', 'TOT_OX', 'AG', 'AL', 'AS_', 'AU', 'B', 'BA', 'BE', 'BI', 'BR', 'C', 'CA', 'CD', 'CE', 'CL', 'CO', 'CR', 'CS', 'CU', 'DY', 'ER', 'EU', 'F', 'FE', 'GA', 'GD', 'GE', 'HF', 'HG', 'HO', 'IR', 'K', 'LA', 'LI', 'LU', 'MG', 'MN', 'MO', 'F_99999', 'NB', 'ND', 'NI', 'OS', 'P', 'PB', 'PD', 'PR', 'PT', 'RB', 'RE', 'RH', 'RU', 'S', 'SB', 'SC', 'SE', 'SM', 'SN', 'SR', 'TA', 'TB', 'TE', 'TH', 'TI', 'TL', 'TM', 'U', 'V', 'W', 'Y', 'YB', 'ZN', 'ZR']

    oxides = ['SIO2', 'TIO2', 'AL2O3', 'FE2O3', 'FEO', 'MNO', 'MGO', 'CAO', 'NA2O', 'K2O', 'P2O5', 'SO3']
    trace_elems = ['AG', 'AL', 'AS_', 'AU', 'B', 'BA', 'BE', 'BI', 'BR', 'C', 'CA', 'CD', 'CE', 'CL', 'CO', 'CR', 'CS', 'CU', 'DY', 'ER', 'EU', 'F', 'FE', 'GA', 'GD', 'GE', 'HF', 'HG', 'HO', 'IR', 'K', 'LA', 'LI', 'LU', 'MG', 'MN', 'MO', 'F_99999', 'NB', 'ND', 'NI', 'OS', 'P', 'PB', 'PD', 'PR', 'PT', 'RB', 'RE', 'RH', 'RU', 'S', 'SB', 'SC', 'SE', 'SM', 'SN', 'SR', 'TA', 'TB', 'TE', 'TH', 'TI', 'TL', 'TM', 'U', 'V', 'W', 'Y', 'YB', 'ZN', 'ZR']

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_workload = None

    in_buff = None
    in_ocean = None
    mapsymbols = None
    df = None

    if(nproc==1):
        print('At least two MPI processors are needed. Please launch the program with MPI on two on more processors. Aborting..')
        exit(0)
    # end if

    if(rank == 0):
        bbox=None
        #bbox = (1121457,-3250284,1327617,-3040209)
        
        print('Reading: %s..'%(master_shapefile))
        df = geopandas.read_file(master_shapefile, bbox=bbox)
        
        print('Reading: %s..'%(geology_lines_buffer_shapefile))
        buff = geopandas.read_file(geology_lines_buffer_shapefile, bbox=bbox)
        
        tree = index.Index()
        for i, geom in tqdm(enumerate(buff.geometry), desc='Building spatial index for %s..'%(geology_lines_buffer_shapefile)):
            tree.insert(i, geom.bounds)
        # end for
        
        print('Locating sites that fall within buffer regions..')
        in_buff = np.bool_(np.zeros(len(df.geometry)))
        inbuff_count = 0
        for i, site in enumerate(df.geometry):
            site_bounds = [site.coords[0][0]-1e-5, site.coords[0][1]-1e-5, 
                           site.coords[0][0]+1e-5, site.coords[0][1]+1e-5]
            igeoms = tree.intersection(site_bounds)
            
            for igeom in igeoms:
                if(buff.geometry[igeom].contains(site)):
                    in_buff[i] = True
                    inbuff_count += 1                            
                    break
                # end if
            # end for
        # end for
        print ('Found %d sites within buffer regions..'%(inbuff_count))

        in_ocean = np.bool_(df['Mask2019_j'] == -9999)
        mapsymbols = np.array(df['MAPSYMBOL'])
    # end if

    in_buff = comm.bcast(in_buff, root=0)
    in_ocean = comm.bcast(in_ocean, root=0)
    mapsymbols = comm.bcast(mapsymbols, root=0)
    comm.barrier()

    payload = None
    if(rank==0):
        work_items = split_list(species, nproc-1)
        
        for i in np.arange(1, nproc):
            payload = defaultdict(list)
            for work_item in work_items[i-1]:
                payload[work_item] = np.array(df[work_item])
            # end for
            else: comm.send(payload, dest=i)
        # end for
    else:
        payload = comm.recv(source=0)
    # end if
    
    comm.barrier()
    print('Received payload on all worker ranks..')

    if(rank>0):
        # start processing
        local_results = defaultdict(list)
        for sname, vals in payload.iteritems():
            result_dict = process_species(rank, sname, vals, mapsymbols, in_ocean)

            for (k, v) in result_dict.iteritems():
                local_results[k] = v
            # end for
        # end for
        comm.send(local_results, dest=0)
    # end if
    
    if(rank == 0):
        for i in np.arange(1, nproc):
            received_results = comm.recv(source=i)
            for colname, vals in received_results.iteritems():
                df[colname] = vals
            # end for
        # end for
        
        # set weights for in_buffer regions
        df['Wbuff'] = np.ones(len(df))
        df.loc[in_ocean, 'Wbuff'] = -99999
        df.loc[in_buff, 'Wbuff'] = 0.7
        
        # set Wg1 for oxides
        for ox in oxides:
            colname = gen_colname(ox, 'Wg1')
            
            df[colname] = np.zeros(len(df))    
            df.loc[(df[ox]>=0) & (df[ox]<=100), colname] = 1
        # end for        
        
        # set Wg1 for trace elements
        for telem in trace_elems:
            colname = gen_colname(telem, 'Wg1')
            
            df[colname] = np.zeros(len(df))    
            df.loc[(df[telem]>=0) & (df[telem]<=1000000), colname] = 1
        # end for        
        
        # set Wg2 for tot_ox
        colname = 'Wg2'
        df[colname] = np.zeros(len(df))

        df.loc[(df['TOT_OX']<90) | (df['TOT_OX']>110), colname] = 0.5
        df.loc[(df['TOT_OX']>=90) & (df['TOT_OX']<=110), colname] = 0.8
        df.loc[(df['TOT_OX']>=95) & (df['TOT_OX']<=105), colname] = 0.9
        df.loc[(df['TOT_OX']>=97) & (df['TOT_OX']<=103), colname] = 1

        df.loc[(df['TOT_OX']==-99999), colname] = 1.0        
        
        # stage two
        process_collocated(df, mapsym_raster, mapsym_grid, species)
        
        # output file
        ofn = output_file_basename + '.shp'
        print ('Writing output shapefile..')
        df.to_file(ofn)
    # end if
# end func

if __name__ == "__main__":
    process()
# end if

