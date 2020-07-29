from __future__ import print_function, division
import sys, os
import astropy.io.fits as fits
import numpy as np
from  martxclib import martxcfun, martevtlister
from astropy.table import Table as tab
from astropy import wcs
from astropy.time import Time, TimeDelta
from astropy.io import fits as pyfits
import tqdm
from matplotlib import pyplot as plt
from time import gmtime, strftime
import glob


tilenumber = sys.argv[1]

test_dir='/mnt/artxc/data/swartz/level2/'

toglob = test_dir + "/%s" %(tilenumber)
print(toglob)

evtlists = glob.glob(test_dir + "%s/*_urd.fits" %(tilenumber))

#print(evtlists)

ratio_sciband_list = []
ratio_bgband_list = []

for this_evtfile in tqdm.tqdm(evtlists):
    
    this_evtlister = martevtlister.martEventLister(this_evtfile)
    
    this_flag0data = this_evtlister.FilterData_Flag(this_evtlister.evt_data_original,0)
    this_flag2data = this_evtlister.FilterData_Flag(this_evtlister.evt_data_original,2)
    
    flag0_sciband = this_evtlister.FilterData_Energy(this_flag0data,4,30)
    flag2_sciband = this_evtlister.FilterData_Energy(this_flag2data,4,30)  
    
    flag0_bgband = this_evtlister.FilterData_Energy(this_flag0data,30,85)
    flag2_bgband = this_evtlister.FilterData_Energy(this_flag2data,30,85)  
    
    flag0_scicts = np.shape(flag0_sciband)[0] / this_evtlister.flag0_area
    flag0_bgcts = np.shape(flag0_bgband)[0] / this_evtlister.flag0_area
    
    flag2_scicts = np.shape(flag2_sciband)[0] / this_evtlister.flag2_area
    flag2_bgcts = np.shape(flag2_bgband)[0] / this_evtlister.flag2_area
    
    ratio_sciband = flag0_scicts / flag0_bgcts
    ratio_bgband = flag2_scicts / flag2_bgcts
    
    ratio_sciband_list.append(ratio_sciband)
    ratio_bgband_list.append(ratio_bgband)
    
ratio_sciband_arr = np.array(ratio_sciband_list)
ratio_bgband_arr = np.array(ratio_bgband_list)

plt.scatter(ratio_bgband_arr,ratio_sciband_arr)

sorted_ratios = np.sort(ratio_sciband_arr)
plt.plot(sorted_ratios,sorted_ratios)
plt.xlabel("Ratio ((4-30 keV) / (30-85 keV), Flag2")
plt.ylabel("Ratio ((4-30 keV) / (30-85 keV), Flag0")
plt.show()
