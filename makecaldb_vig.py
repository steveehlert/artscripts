from astropy.table import Table as tab
from astropy.io import ascii, fits
import pandas as pd
import numpy as np
import os
version = '0.0.2'


eapath = '/home/ctchen/art-xc/caldbdata/EA'
destpath = '/home/ctchen/art-xc/caldbdata'

#%%
theta_list = [-18, -15, -12, -9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9, 12, 15, 18]
hduid = [8, 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16]
theta_list7 = [-27, -24, -21, -18, -15, -12, -9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9, 12, 15, 18, 21, 24, 27]
hduid7 = [11, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19]
detnamlist = ['T01', 'T02', 'T03', 'T04', 'T05', 'T06', 'T07']
mmalist = ['3', '5', '4', '1', '2', '8', '7']
urdnlist = [28, 22, 23, 24, 25, 26, 30]
detnlist = [33, 34, 24, 25, 22, 26, 23]
detnam_d = {mmalist[i]: detnamlist[i] for i in range(len(detnamlist))}
urdn_d = {mmalist[i]: urdnlist[i] for i in range(len(urdnlist))}
detn_d = {mmalist[i]: detnlist[i] for i in range(len(detnlist))}


def func_EA(theta,thetac,sig1,a1,sig2,a2):
    return a1*np.exp(-(theta-thetac)**2/(2*sig1**2)) + a2*np.exp(-(theta-thetac)**2/(2*sig2**2))
    
energies = np.array([4,5,6,7,8,9,10,12,14,16,18,20,25,30],dtype=float)
en_lo = energies[0:-1]
en_hi = energies[1:]

thetadeg = np.array([0.        , 0.01666667, 0.03333334, 0.05      , 0.06666667,
        0.08333334, 0.1       , 0.11666667, 0.13333334, 0.15      ,
        0.16666667, 0.18333334, 0.2       , 0.21666667, 0.23333333,
        0.25      , 0.26666668, 0.28333333, 0.3       , 0.31666666,
        0.33333334, 0.35      , 0.36666667, 0.38333333, 0.4       ,
        0.41666666, 0.43333334, 0.45      , 0.46666667, 0.48333332,
        0.5       ])
thetaarcmin = thetadeg * 60.
for istr in mmalist:
    #istr = str(i)
    detnam = detnam_d[istr]
    urdn = urdn_d[istr]
    detn = detn_d[istr]
    csvname = eapath + '/ARTXC_EA_M' + istr + '.csv'
    df = pd.read_csv(csvname)
    vigarr = np.zeros((1,1,31,13))
    for i in range(13):
        thetac = df.loc[i,'thetac']
        sig1 = df.loc[i,'sigma1']
        sig2 = df.loc[i,'sigma2']
        a1 = df.loc[i,'A1']
        a2 = df.loc[i,'A2']
        vigarr[0,0,:,i] = func_EA(thetaarcmin, thetac, sig1, a1, sig2, a2)/\
        np.max(func_EA(thetaarcmin, thetac, sig1, a1, sig2, a2))
    t = tab([[en_lo], [en_hi], [thetadeg], np.array([0.0]), vigarr], names=('ENERG_LO', 'ENERG_HI','THETA', 'PHI','VIGNET'))
    col1 = fits.Column(array=t['ENERG_LO'].data,name='ENERG_LO',format='13E',unit='keV')
    col2 = fits.Column(array=t['ENERG_HI'].data,name='ENERG_HI',format='13E',unit='keV')
    col3 = fits.Column(array=t['THETA'].data,name='THETA',format='31E',unit='degree')
    col4 = fits.Column(array=t['PHI'].data,name='PHI',format='1E',unit='degree')
    col5 = fits.Column(array=t['VIGNET'].data,name='VIGNET',format='403E',dim='(13, 31, 1)')
    coldefs = fits.ColDefs([col1,col2,col3,col4,col5])
    hdu = fits.BinTableHDU.from_columns(coldefs)
    hdu.header['EXTNAME']  = 'VIGNET  '
    hdu.header['HDUCLASS'] = 'OGIP'
    hdu.header['INSTRUME'] = 'ART-XC'
    hdu.header['VIGVERSN'] = '1992a'
    hdu.header['HDUCLAS1'] = 'RESPONSE'
    hdu.header['HDUVERS1'] = '1.0.0   '
    hdu.header['HDUCLAS2'] = 'VIGNET  '
    hdu.header['HDUVERS2'] = '1.1.0   '
    hdu.header['CSYSNAME'] = 'XNA_POL'
    hdu.header['CREATOR'] = 'makecaldb_vig.py v' + version
    hdu.header['ORIGIN'] = 'NASA MSFC'
    hdu.header['TELESCOP'] = 'SRG'
    hdu.header['DETN']  = detn
    hdu.header['URDN']  = urdn
    hdu.header['DETNAM']  = detnam
    hdu.header['CCLS0001'] = 'BCF'
    hdu.header['CDTP0001'] = 'DATA'
    hdu.header['CCNM0001'] = 'VIGNET'
    hdu.header['CBD10001'] = 'THETA(0-0.5)deg'
    hdu.header['CBD20001'] = 'PHI(0-360.0)deg'
    hdu.header['CBD30001'] = 'ENERG (4-30.0)keV'
    hdu.header['CVSD0001'] = "2013-10-01"
    hdu.header["CVST0001"] = "00:00:00"
    hdu.header["CDES0001"] = "Vignetting functions based on MSFC calibration database 2013-2014/"
    outname = destpath + '/art' + detnam + '_vign_20131001v0002.fits'
    hdu.writeto(outname,overwrite=True, checksum=True)
