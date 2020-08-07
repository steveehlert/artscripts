from martxclib.martimagetools import martImager, martExpmapper
from martxclib import martxcbgdtools as bgtools
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.io import fits
import sys
from scipy.optimize import minimize
from scipy.special import factorial, gammaincc, gamma
from scipy import stats
# this can be changed for machines other than camelot
#datapath = '/mnt/artxc/data/swartz'
datapath = '/home/ctchen/art-xc'
import tqdm

telstring="T0%s" %(sys.argv[1])

all_fov_images = glob.glob("/mnt/artxc/data/swartz/workspace/imaging/273024/img_2*%s*.fits" %(telstring))
all_corner_images = glob.glob("/mnt/artxc/data/swartz/workspace/imaging/273024/img_corner_2*%s*.fits" %(telstring))

all_fov_expmaps = glob.glob("/mnt/artxc/data/swartz/workspace/expmap/273024/expfov*%s.fits" %(telstring))
all_corner_expmaps = glob.glob("/mnt/artxc/data/swartz/workspace/expmap/273024/expcorner*%s*.fits" %(telstring))
all_fov_nvexpmaps = glob.glob("/mnt/artxc/data/swartz/workspace/expmap/273024/expfov*%s.nv.fits" %(telstring))

all_fov_images.sort()
all_corner_images.sort()
all_fov_expmaps.sort()
all_corner_expmaps.sort()
all_fov_nvexpmaps.sort()
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
sumexp_fovnv = None
final_header_exp = None
j=0

for img_fov, img_corner, exp_fov, exp_corner, exp_fovnv in tqdm.tqdm(zip(all_fov_images, all_corner_images, all_fov_expmaps, all_corner_expmaps, all_fov_nvexpmaps)):
    print(img_fov, img_corner)
    print(exp_fov, exp_corner, exp_fovnv)
    
    img_fov_header = fits.open(img_fov)[0].header
    img_fov_data = fits.open(img_fov)[0].data
    img_corner_data = fits.open(img_corner)[0].data
    
    exp_fov_data = fits.open(exp_fov)[0].data
    exp_corner_data = fits.open(exp_corner)[0].data
    exp_fovnv_data = fits.open(exp_fovnv)[0].data
    exp_fov_header = fits.open(exp_fov)[0].header
    
    
    if j==0: 
        sumphotons_fov = img_fov_data
        sumphotons_corner = img_corner_data
        sumexp_fov  = exp_fov_data.astype('float')
        sumexp_corner = exp_corner_data.astype('float')
        sumexp_fovnv = exp_fovnv_data.astype('float')
        final_header = img_fov_header
        final_header_exp = exp_fov_header
    else:
        sumphotons_fov += img_fov_data
        sumphotons_corner += img_corner_data
        sumexp_fov  += exp_fov_data.astype('float')
        sumexp_corner += exp_corner_data.astype('float')
        sumexp_fovnv += exp_fovnv_data.astype('float')
    j+=1
    
    
#final_expcor = sumphotons / sum_exp
expcor_corner = sumphotons_corner/sumexp_corner

expcor_corner[np.where(sumexp_corner < 1.e-3)] =0.


expmask = sumexp_corner > 0
#result = minimize(negative_log_likelihood,  # function to minimize
#                  x0=0,            # start value
#                  args=(sumphotons_corner[expmask],sumexp_corner[expmask]),             # additional arguments for function
#                  method='Powell',          # minimization method, see docs
#                  )
#print("MLE bg count rate: ", 10**(result['x'][0]))

mean_bg_ctrate = np.mean(sumphotons_corner[expmask]/sumexp_corner[expmask])

print("Simple Mean BG Count Rate: ",mean_bg_ctrate )





expected_bg = sumexp_fovnv * mean_bg_ctrate

#prob_bg = gammaincc(sumphotons_fov, expected_bg) / gamma(sumphotons_fov)


#prob_sky = 1. - prob_bg

non_zero_mask = sumphotons_fov > 0
zero_mask = sumphotons_fov ==0

#max_counts = np.max(sumphotons_fov,axis=None)
clean_image = np.copy(sumphotons_fov)
err_image = np.zeros_like(clean_image)

expected_bg_nz = expected_bg[non_zero_mask]
sumphotons_fov_nz = sumphotons_fov[non_zero_mask]

mu_star_list = []
mu_err_list = []
for mu_n, k in zip(expected_bg_nz, sumphotons_fov_nz):
    
    k = int(k)
    mu_s_bar_0 = (-1. * mu_n + np.sqrt(mu_n**2 + 4.) )/2.
    prob_0 = (1. - mu_n/np.sqrt(mu_n**2 + 4.) )**k
    
    mu_s_list = [mu_s_bar_0]
    weights_list = [prob_0]
    
    for k_n in range(1,k+1):
        mu_s_bar = (k -k_n) / k_n * mu_n
        prob = (k_n**k_n * (k - k_n)**(k - k_n))/(k**k)
        mu_s_list.append(mu_s_bar)
        weights_list.append(prob)
    
    mu_s_arr = np.array(mu_s_list)
    weight_arr = np.array(weights_list)
    mu_s_star = np.average(mu_s_arr,weights=weight_arr)
    err_mu = np.sqrt(np.cov(mu_s_arr,aweights = weight_arr))
    mu_star_list.append(mu_s_star)
    mu_err_list.append(err_mu)
    
mu_star_arr = np.array(mu_star_list)
mu_err_arr = np.array(mu_err_list)

clean_image[zero_mask] = 0
clean_image[non_zero_mask] = mu_star_arr

err_image[non_zero_mask] = mu_err_arr

#clean_image = np.copy(sumphotons_fov)

#clean_image[np.where(fov_rand < prob_bg)] =0.




expcor_fov = clean_image /sumexp_fov
#expcor_fov[np.where(sumexp_fov < 1.e-3)] = 0.

expcor_fov_err = err_image / sumexp_fov


#plt.imshow(expcor_fov)
#plt.show()

expmap_vig_hdu = fits.ImageHDU(data = sumexp_fov, header = final_header_exp)
expmap_vig_hdu.writeto("Summed_ExpMap_VIG_273024_%s.fits" %(telstring),overwrite=True)

expmap_nv_hdu = fits.ImageHDU(data = sumexp_fovnv, header = final_header_exp)
expmap_nv_hdu.writeto("Summed_ExpMap_NV_273024_%s.fits" %(telstring),overwrite=True)

final_hdu = fits.ImageHDU(data = expcor_fov, header = final_header)
final_hdu.writeto("Summed_ExpCorr_273024_%s.fits" %(telstring),overwrite=True)

finalerr_hdu = fits.ImageHDU(data = expcor_fov_err, header = final_header)
finalerr_hdu.writeto("Summed_ExpCorrErr_273024_%s.fits" %(telstring),overwrite=True)






