import numpy as np
from astropy.io import fits
import glob


image_list = []
image_header = None


seven_files = glob.glob("Summed_ExpCorr_273024_T*.fits")

for this_file in seven_files:

    image_hdu = fits.open(this_file)[1]

    image_list.append(image_hdu.data)

    image_header = image_hdu.header



image_arr  = np.array(image_list)

med_image = np.median(image_arr,axis=0)

final_hdu = fits.ImageHDU(data = med_image, header = image_header)
final_hdu.writeto("ExpCor_273024_AllSeven.fits",overwrite=True)

