from martxclib.martimagetools import martImager, martExpmapper
from martxclib import martxcbgdtools as bgtools
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.io import fits
import sys
from scipy.optimize import minimize
from scipy.special import factorial, gammaincc
from scipy import stats
# this can be changed for machines other than camelot
#datapath = '/mnt/artxc/data/swartz'
datapath = '/home/ctchen/art-xc'
import tqdm

telstring="T0%s" %(sys.argv[1])

all_fov_images = glob.glob("/mnt/artxc/data/swartz/workspace/imaging/273024/img_2*%s*.fits" %(telstring))
all_corner_images = glob.glob("/mnt/artxc/data/swartz/workspace/imaging/273024/img_corner_2*%s*.fits" %(telstring))

all_fov_expmaps = glob.glob("/mnt/artxc/data/swartz/workspace/expmap/273024/expfov*%s*.fits" %(telstring))
all_corner_expmaps = glob.glob("/mnt/artxc/data/swartz/workspace/expmap/273024/expcorner*%s*.fits" %(telstring))


all_fov_images.sort()
all_corner_images.sort()
all_fov_expmaps.sort()
all_corner_expmaps.sort()

j=0




def poisson(k, lamb):
    """poisson pdf, parameter lamb is the fit parameter"""
    return (lamb**k/factorial(k)) * np.exp(-lamb)



def negative_log_likelihood(params, data, expmap):
    ''' better alternative using scipy '''
    #print(params[0]*expmap)
    return -stats.poisson.logpmf(data, np.exp(params[0])*expmap).sum()



sumphotons_fov= None
sumexp_fov = None
sumphotons_corner = None
sumexp_corner = None
final_header = None
j=0

for img_fov, img_corner, exp_fov, exp_corner in tqdm.tqdm(zip(all_fov_images, all_corner_images, all_fov_expmaps, all_corner_expmaps)):
    #print(img_fov, img_corner)
    #print(exp_fov, exp_corner)
    
    img_fov_header = fits.open(img_fov)[0].header
    img_fov_data = fits.open(img_fov)[0].data
    img_corner_data = fits.open(img_corner)[0].data
    
    exp_fov_data = fits.open(exp_fov)[0].data
    exp_corner_data = fits.open(exp_corner)[0].data
    
    
    
    #Loading up martImager
    #imager = martImager(evtname,datapath=datapath)
    #THis will get updated with better L2 times
    #flaretool = bgtools.martFlareTool(evtname)
    #print(flaretool.timebins_all,flaretool.timebins_gti)
    #print(flaretool.gti_data)
    #flaretool.PlotResults(save_fig=False)

    #clean_data = flaretool.FilterData_FlareTimes(flaretool.evt_data_original) 
    #if len(clean_data) == 0:
    #    clean_data = flaretool.evt_data_original
        #print('Flare filtering results in an empty list, assuming no flares...')
    #Photons are 4-30 keV, gather corners and fov
    #photon_data = imager.FilterData_Energy(clean_data,4,30) 
    #idata_fov = imager.FilterData_Flag(photon_data, 0)
    #idata_corner = imager.FilterData_Flag(photon_data, 2)

    #Particle data is 30-85 keV, gather corners and fov
    #particle_data = imager.FilterData_Energy(clean_data,30,85)
    #pdata_fov = imager.FilterData_Flag(particle_data,0)
    #pdata_corner = imager.FilterData_Flag(particle_data,2)

    #Ratio of cts/raw_pixel  in corners and fov
    #inout_ratio = float(len(pdata_fov)) / float(imager.evtlister.flag0_area) / ( float(len(pdata_corner)) / float(imager.evtlister.flag2_area))     


    # Make some images
    #ihdu_fov = imager.Make_Image_SkyXY(idata_fov)

    #ihdu_corner = imager.Make_Image_SkyXY(idata_corner)

    #try:
    #    exp_fov = imager.Make_Expmap(imagehdu=ihdu_fov,flag=0,vig=True)

    #    exp_corner = imager.Make_Expmap(imagehdu=ihdu_corner,flag=2,vig=False)
    #except: 
    #    continue


    #expmask = np.where(exp_corner.data > 0)
    

    
    
    #ctrate_nxb = np.mean(ihdu_corner.data[expmask]) #np.exp(result['x'])[0]
    #print("Simple Average Count Rate: ", ctrate_nxb)
    #nxb_fov = ctrate_nxb * exp_fov.data


    #fov_nxbprob = poisson(ihdu_fov.data, nxb_fov)
    #plt.imshow(fov_nxbprob)
    #plt.show()

    #fov_rand = stats.uniform.rvs(size=np.shape(fov_nxbprob))
    
    
    #bgmask =  fov_nxbprob >= fov_rand
    
    

    #fov_nxbsub = np.copy(ihdu_fov.data)
    
    #fov_nxbsub[np.where(bgmask)] =0.

    
    if j==0: 
        sumphotons_fov = img_fov_data
        sumphotons_corner = img_corner_data
        sumexp_fov  = exp_fov_data.astype('float')
        sumexp_corner = exp_corner_data.astype('float')
        final_header = img_fov_header
    else:
        sumphotons_fov += img_fov_data
        sumphotons_corner += img_corner_data
        sumexp_fov  += exp_fov_data.astype('float')
        sumexp_corner += exp_corner_data.astype('float')
    j+=1
    
    
#final_expcor = sumphotons / sum_exp
expcor_corner = sumphotons_corner/sumexp_corner

expcor_corner[np.where(sumexp_corner < 1.e-3)] =0.


expmask = sumexp_corner > 0
result = minimize(negative_log_likelihood,  # function to minimize
                  x0=0,            # start value
                  args=(sumphotons_corner[expmask],sumexp_corner[expmask]),             # additional arguments for function
                  method='Powell',          # minimization method, see docs
                  )
print("MLE bg count rate: ", 10**(result['x'][0]))

mean_bg_ctrate = np.mean(sumphotons_corner[expmask]/sumexp_corner[expmask])

print("Simple Mean BG Count Rate: ",mean_bg_ctrate )


fov_rand = stats.uniform.rvs(size=np.shape(sumphotons_fov))

#plt.hist(expcor_corner.flatten(),bins=100)
#plt.semilogy()
#plt.show()

#plt.imshow(sumphotons_corner)
#plt.show()

#plt.imshow(sumexp_corner)
#plt.show()


#plt.imshow(expcor_corner)
#plt.show()


#plt.imshow(sumphotons_fov)
#plt.show()

#plt.imshow(sumexp_fov)
#plt.show()

expected_bg = sumexp_fov * mean_bg_ctrate

prob_bg = gammaincc(sumphotons_fov, expected_bg)


prob_sky = 1. - prob_bg

non_zero_mask = sumphotons_fov > 0
zero_mask = sumphotons_fov ==0

clean_image = np.copy(sumphotons_fov)

clean_image[zero_mask] = 0
clean_image[non_zero_mask] = stats.binom.rvs(sumphotons_fov[non_zero_mask].astype('int'),prob_sky[non_zero_mask])



#clean_image = np.copy(sumphotons_fov)

#clean_image[np.where(fov_rand < prob_bg)] =0.




expcor_fov = clean_image /sumexp_fov
expcor_fov[np.where(sumexp_fov < 1.e-3)] = 0.

#plt.imshow(expcor_fov)
#plt.show()

final_hdu = fits.ImageHDU(data = expcor_fov, header = final_header)
final_hdu.writeto("Summed_ExpCorr_273024_%s.fits" %(telstring),overwrite=True)








