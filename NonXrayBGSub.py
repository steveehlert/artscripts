from martxclib.martimagetools import martImager, martExpmapper
from martxclib import martxcbgdtools as bgtools
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.io import fits
import sys
from scipy.optimize import minimize
from scipy.special import factorial
from scipy import stats
# this can be changed for machines other than camelot
#datapath = '/mnt/artxc/data/swartz'
datapath = '/home/ctchen/art-xc'


telstring="T0%s" %(sys.argv[1])

all_evtnames = glob.glob(datapath + '/level2/273024/artl2_'+telstring+'*.fits')
# can also try one of the large files : artl2_T04_20191211_183844.fits
print(len(all_evtnames))
print('Initilizing...')

j=0




def poisson(k, lamb):
    """poisson pdf, parameter lamb is the fit parameter"""
    return (lamb**k/factorial(k)) * np.exp(-lamb)



def negative_log_likelihood(params, data, expmap):
    ''' better alternative using scipy '''
    #print(params[0]*expmap)
    return -stats.poisson.logpmf(data, np.exp(params[0])*expmap).sum()



sumphotons= None
sum_exp = None
final_header = None
for evtname in all_evtnames:
    print(evtname)
    #Loading up martImager
    imager = martImager(evtname,datapath=datapath)
    #THis will get updated with better L2 times
    flaretool = bgtools.martFlareTool(evtname)
    #print(flaretool.timebins_all,flaretool.timebins_gti)
    #print(flaretool.gti_data)
    #flaretool.PlotResults(save_fig=False)

    clean_data = flaretool.FilterData_FlareTimes(flaretool.evt_data_original) 
    if len(clean_data) == 0:
        clean_data = flaretool.evt_data_original
        #print('Flare filtering results in an empty list, assuming no flares...')
    #Photons are 4-30 keV, gather corners and fov
    photon_data = imager.FilterData_Energy(clean_data,4,30) 
    idata_fov = imager.FilterData_Flag(photon_data, 0)
    idata_corner = imager.FilterData_Flag(photon_data, 2)

    #Particle data is 30-85 keV, gather corners and fov
    particle_data = imager.FilterData_Energy(clean_data,30,85)
    pdata_fov = imager.FilterData_Flag(particle_data,0)
    pdata_corner = imager.FilterData_Flag(particle_data,2)

    #Ratio of cts/raw_pixel  in corners and fov
    inout_ratio = float(len(pdata_fov)) / float(imager.evtlister.flag0_area) / ( float(len(pdata_corner)) / float(imager.evtlister.flag2_area))     


    # Make some images
    ihdu_fov = imager.Make_Image_SkyXY(idata_fov)

    ihdu_corner = imager.Make_Image_SkyXY(idata_corner)

    try:
        exp_fov = imager.Make_Expmap(imagehdu=ihdu_fov,flag=0,vig=True)

        exp_corner = imager.Make_Expmap(imagehdu=ihdu_corner,flag=2,vig=False)
    except: 
        continue

    #print("FOV Exposure Map")
    #plt.imshow(exp_fov.data, norm=LogNorm(vmin=np.max(exp_fov.data)*1e-6))
    #plt.show()

    #print("Corner Exposure Map")
    #plt.imshow(exp_corner.data, norm=LogNorm(vmin=np.max(exp_corner.data)*1e-6))
    #plt.show()
    
    #plt.hist(ihdu_corner.data.flatten(), bins=10)
    ##plt.semilogy()
    #plt.show()

    #expcor_corner = ihdu_corner.data.astype('float') / exp_corner.data.astype('float')

    expected_bg = sumexp_fov * mean_bg_ctrate

prob_bg = gammaincc(np.ones_like(sumphotons_fov), expected_bg)


prob_sky = 1. - prob_bg

non_zero_mask = sumphotons_fov > 0
zero_mask = sumphotons_fov ==0

clean_image = np.copy(sumphotons_fov)

clean_image[zero_mask] = 0
clean_image[non_zero_mask] = stats.binom.rvs(sumphotons_fov[non_zero_mask].astype('int'),prob_sky[non_zero_mask])


    expmask = np.where(exp_corner.data > 0)
    
    #mean_ctrate = np.mean(
    

    #expcor_corner[np.where(exp_corner.data < 1e-6) ] =0.
    #plt.hist(expcor_corner.flatten(),bins=100)
    #plt.semilogy()
    #plt.show()

    #print("Exposure Corrected Corner Image")
    #plt.imshow(expcor_corner)
    #plt.show()
    
    #This is the expected particle background in the 4-30 keV band in the FOV
    #plt.hist(expcor_corner[np.where(exp_corner.data >= 0)],bins=100)
    #plt.semilogy()
    #plt.show()
    #print(np.min(ihdu_corner.data[expmask],exp_corner.data[expmask])
    #result = minimize(negative_log_likelihood,  # function to minimize
    #              x0=0,            # start value
    #              args=(ihdu_corner.data[expmask],exp_corner.data[expmask]),             # additional arguments for function
    #              method='Powell',          # minimization method, see docs
    #              )
    #print("MLE bg count rate: ", 10**(result['x'][0]))
    #ctrate_nxb = inout_ratio * np.median(expcor_corner[np.where(exp_corner.data >= 0)])
    #ctrate_nxb = 0.
    
    
    
    ctrate_nxb = np.mean(ihdu_corner.data[expmask]) #np.exp(result['x'])[0]
    #print("Simple Average Count Rate: ", ctrate_nxb)
    nxb_fov = ctrate_nxb * exp_fov.data
    
    #plt.imshow(nxb_fov)
    #plt.show()
    #plt.imshow(ihdu_fov.data)
    #plt.show()

    fov_nxbprob = poisson(ihdu_fov.data, nxb_fov)
    #plt.imshow(fov_nxbprob)
    #plt.show()

    fov_rand = stats.uniform.rvs(size=np.shape(fov_nxbprob))
    
    
    bgmask =  fov_nxbprob >= fov_rand
    
    

    fov_nxbsub = np.copy(ihdu_fov.data)
    
    fov_nxbsub[np.where(bgmask)] =0.
    
    #plt.imshow(fov_nxbsub)
    #plt.show()
    
    
    if j==0: 
        sumphotons = fov_nxbsub
        sum_exp = exp_fov.data.astype('float')
        final_header = ihdu_fov.header
    else:
        sumphotons += fov_nxbsub
        sum_exp += exp_fov.data.astype('float')
    j+=1
    
    
final_expcor = sumphotons / sum_exp
#plt.imshow(final_expcor)
#plt.show()
#plt.hist(final_expcor.flatten(),bins=100)
#plt.semilogy()
#plt.show()

final_hdu = fits.ImageHDU(data = final_expcor, header = final_header)
final_hdu.writeto("Summed_ExpCorr_273024_%s.fits" %(telstring),overwrite=True)
    #fov_bgsubexpcor = fov_nxbsub / exp_fov.data.astype('float')

    #fov_bgsubexpcor[np.where(exp_fov.data < 1e-6)] =0.

    #print("Fully BG subbed exposure corrected image")
    #plt.imshow(fov_bgsubexpcor)
    #plt.show()







