import numpy as np
import os
import operator
import sys
import glob
import re
import csv
import pandas as pd
np.set_printoptions(threshold=np.inf)
from astropy.coordinates import SkyCoord,EarthLocation,AltAz
from astropy.wcs import WCS 
from astropy import units as u 
from astropy.table import Table 
from astropy.io import ascii
from astropy.time import Time

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import argparse

def err_doppler(delta_v,c):
    return delta_v/c

def DM_delay(p,flow,fhigh):
    return 0.1*p/(4.15*((1/flow)**2-(1/fhigh)**2))


### speed of light (km/s)
c = 300_000

#### Earth's orbit and rotation (km/s)
### 2 comes from the highest
v_orb = 2 * 30
v_rot = 2 * 0.460

### Doppler from redback (km/s)
v_red = 2 * 300

##### Pdot (<0.01, 0.01-1, 1<)
pdot_msp = 1E-20
pdot_fast = 1E-16
pdot_norm = 1E-14


def parseOptions(parser):
    parser.add_argument('-f', dest='file', 
				nargs = 2, 
				type = str, 
				help="Input for two csv files",
				required=True)
    parser.add_argument('--s', dest='separation', 
				nargs = 1, 
				type = float, 
				help="A factor for sky separation",
                                default = 1.,
				required=False)
    parser.add_argument('--d', dest='dm', 
				nargs = 1, 
				type = float, 
				help="A factor for DM",
                                default = 1.,
				required=False)
    parser.add_argument('--p', dest='period', 
				nargs = 1, 
				type = float, 
				help="A factor for period",
                                default = 1.)
    

    options= parser.parse_args()

    return options




#pointing_id, ra, dec, f0, dm = np.loadtxt('candidates.csv',str,unpack=True,skiprows=1,usecols=(0,4,5,11,20),delimiter=',')
#print(pointing_id[0])
#print(ra[0])
#print(dec[0])
#print(f0[0])
#print(dm[1])


def matching():
    
    pointing_id1, ra1, dec1, time1, f01, dm1 = np.loadtxt(options.file[0],str,unpack=True,skiprows=1,usecols=(0,4,5,8,11,20),delimiter=',')
    f01 = f01.astype(float)
    dm1 = dm1.astype(float)
    time1 = time1.astype(float)

    pointing_id2, ra2, dec2, time2, f02, dm2 = np.loadtxt(options.file[1],str,unpack=True,skiprows=1,usecols=(0,4,5,8,11,20),delimiter=',')
    f02 = f02.astype(float)
    dm2 = dm2.astype(float)
    time2 = time2.astype(float)
    #print(f02)
    
    #source_name, beamshape_y, beamshape_x, obs_id, cen_freq, bw, obs_time = np.loadtxt(options.beam[0],str,unpack=True,skiprows=1,usecols=(1,6,7,12,9,10,4),delimiter=',')
    #cen_freq = cen_freq.astype(float)
    #bw = bw.astype(float)
    #print(beam_id1[0])
    '''
    for i in range(0,len(pointing_id1)):
        #print(options.file[0][0:6].replace('_',''))
        if options.file[0][0:6].replace('_','') == obs_id[i][0:7].replace('"','').replace('_',''):
            beam1x = float(beamshape_x[i].replace('"','')[6:])
            beam1y = float(beamshape_y[i].replace('"','').replace('{\y\: ',''))
            freq_min1 = cen_freq[i] - bw[i]/2
            freq_max1 = cen_freq[i] + bw[i]/2
            ### from Hz to GHz
            freq_min1 = freq_min1/10**9
            freq_max1 = freq_max1/10**9
            ### obs time
            time1 = Time(obs_time[i].replace('"',''), scale='utc')
            break
    
        
    for i in range(0,len(pointing_id1)):
        #print(options.file[1][0:6].replace('_',''))
        if options.file[1][0:6].replace('_','') == obs_id[i][0:7].replace('"','').replace('_',''):
            beam2x = float(beamshape_x[i].replace('"','')[6:])
            beam2y = float(beamshape_y[i].replace('"','').replace('{\y\: ',''))
            freq_min2 = cen_freq[i] - bw[i]/2
            freq_max2 = cen_freq[i] + bw[i]/2
            ### from Hz to GHz
            freq_min2 = freq_min2/10**9
            freq_max2 = freq_max2/10**9
            ### obs time
            time2 = Time(obs_time[i].replace('"',''), scale='utc')
            break
   
    dt = time1-time2
    time_diff = np.abs(dt.value)
    print('time_diff = ' + str(time_diff))
   '''

    time1 = np.array(time1)
    time2 = np.array(time2)
    time1 = np.array(time1[:,None])
    time2 = np.array(time2[None,:])
       
    time_diff = np.abs(time1 - time2) * 24 * 3600
    print(time_diff)
    #time_diff = np.abs(dt.value)
    #print('time_diff = ' + str(time_diff))
    
    ra1 = np.array(ra1)
    ra2 = np.array(ra2)
    dec1 = np.array(dec1)
    dec2 = np.array(dec2)
    ra1 = np.array(ra1[:,None])
    dec1 = np.array(dec1[:,None])
    ra2 = np.array(ra2[None,:])
    dec2 = np.array(dec2[None,:])

    ### Sky separation ###
    c1 = SkyCoord(ra1,dec1,unit=(u.deg, u.deg))
    c2 = SkyCoord(ra2,dec2,unit=(u.deg, u.deg))
    print(c1.shape)
    print(c2.shape)

    sep = c1.separation(c2).deg
    #sep_err = np.sqrt((beam1x**2 + beam1y**2)/2)**2 + np.sqrt((beam2y**2 + beam2x**2)/2)**2
    beam1x = 0.007
    beam1y = 0.007
    beam2x = 0.007
    beam2y = 0.007
    sep_err = np.sqrt((beam1x**2 + beam2y**2)/2)**2 + np.sqrt((beam2y**2 + beam2x**2)/2)**2
    chi2_sep = sep**2/sep_err
    print(sep_err)
    boo_sep = chi2_sep <= options.separation



    ### DM ###

    freq_min1 = 0.544 * 10**6
    freq_max1 = 1.088 * 10**6

    # l =856$-$1712
    freq_min2 = 544 * 10**6
    freq_max2 = 1088 * 10**6


    dm1 = np.array(dm1)
    dm2 = np.array(dm2)
    dm1 = np.r_[dm1[:,None]]
    dm2 = np.r_[dm2[None,:]]

    DM_off = np.abs(dm1-dm2)
    #print(freq_min1)
    dm_err1 = DM_delay(1/f01*1000,freq_min1/(10**(9)),freq_max1/(10**(9)))**2
    dm_err1 = np.r_[dm_err1[:,None]]
    dm_err2 = DM_delay(1/f02*1000,freq_min2/(10**(9)),freq_max2/(10**(9)))**2
    dm_err2 = np.r_[dm_err2[None,:]]
    #print(dm_err1)
    chi2_dm = DM_off**2/(dm_err1+dm_err2)
    #print(chi2_dm)
    boo_dm = chi2_dm <= options.dm


    ### period ###

    f01 = np.array(f01)
    f02 = np.array(f02)
    p01 = 1/f01
    p02 = 1/f02
    p01 = np.r_[p01[:,None]]
    p02 = np.r_[p02[None,:]]
    f01 = np.r_[f01[:,None]]
    f02 = np.r_[f02[None,:]]

    #time_diff = 600
    
    diff_p = np.abs(p01 - p02)
    var_p = (time_diff*pdot_fast)**2 + (err_doppler(v_orb+v_rot+v_red,c)*p01)**2 + (err_doppler(v_orb+v_rot+v_red,c)*p02)**2
    chi2_p = diff_p**2/var_p
    #print(chi2_p)
    boo_p = chi2_p <= options.period

    ###
    #id1 = np.array(id1)
    #id2 = np.array(id2)
    #id1 = np.r_[id1[:,None]]
    #id2 = np.r_[id2[None,:]]
    #beam_id1 = np.array(beam_id1)
    #beam_id2 = np.array(beam_id2)
    #beam_id1 = np.r_[beam_id1[:,None]]
    #beam_id2 = np.r_[beam_id2[None,:]]
    

    BOO = boo_sep * boo_dm * boo_p
    INDEX_BOO = BOO.nonzero()
    print(len(INDEX_BOO[0]))
    BOO_sep = boo_sep.nonzero()
    BOO_p = boo_p.nonzero()
    BOO_dm = boo_dm.nonzero()
    print(INDEX_BOO)
    #print(len(BOO_sep[0]))
    #print(len(BOO_p[0]))
    #print(len(BOO_dm[0]))
    #print(INDEX_BOO[0][8])
    #print(INDEX_BOO[1][8])
    #print(p02[0][INDEX_BOO[1][8]])
    #print(p01[INDEX_BOO[0][8]][0])

    ### print result ###
    #snr1 = np.array(snr1[:,None])
    #snr2 = np.array(snr2[None,:])
    #cand_id1 = np.array(cand_id1[:,None])
    #cand_id2 = np.array(cand_id2[None,:])
    #file_id1 = np.array(file_id1[:,None])
    #file_id2 = np.array(file_id2[None,:])
    
    #s_name = source_name[i].replace('"','')
    #oid1 =  options.file[0][0:6].replace('_','')
    #oid2 = options.file[1][0:6].replace('_','')

    #B1 = np.loadtxt(options.file[0],float,unpack=True,skiprows=1,usecols=(26),delimiter=',')
    #B2 = np.loadtxt(options.file[1],int,unpack=True,skiprows=1,usecols=(26),delimiter=',')
        
    txt1=open('r_cross.txt','w')
    print('ra1 ra2 dec1 dec2 f01 f02 dm1 dm2', file=txt1)

    ##(n_sep n_dm n_p = ' + str(len(boo_sep.nonzero()[0])) + ' ' + str(len(boo_dm.nonzero()[0])) + ' ' + str(len(boo_p.nonzero()[0])) + ') ' + '#' +  str(len(ra1)) + '_' + str(len(ra2[0])) + '(' + str(len(ra1) * len(ra2[0])) + ')'
    #print(f01[INDEX_BOO[0][1]][0])
    
    for i in range(0,len(INDEX_BOO[0])):
        print(str(ra1[INDEX_BOO[0][i]][0]) + ' ' + str(ra2[0][INDEX_BOO[1][i]]) + ' ' + str(dec1[INDEX_BOO[0][i]][0]) + ' ' + str(dec2[0][INDEX_BOO[1][i]]) + ' ' + str(f01[INDEX_BOO[0][i]][0]) + ' ' + str(f02[0][INDEX_BOO[1][i]]) + ' ' + str(dm1[INDEX_BOO[0][i]][0]) + ' ' + str(dm2[0][INDEX_BOO[1][i]]), file= txt1)

    #print(len(boo_sep.nonzero()))
    #print(len(boo_p.nonzero())) 
    #print(len(boo_dm.nonzero()))
                       
    
if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    options = parseOptions(parser)
    print('offset_factor = ' + str(options.separation))
    print('dm_factor = ' + str(options.dm))
    print('period_factor = ' + str(options.period))
    matching()


