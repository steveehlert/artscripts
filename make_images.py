from martxclib.preprocessfunc import *
import sys
import numpy as np
import astropy.io.fits as fits
from astropy import wcs
import astropy.coordinates as cd
import astropy.units as u
from martxclib.martxcfun import *
from martxclib.martxcfits import *
from martxclib.martl1tools import *
from martxclib.martimagetools import martImager, martExpmapper
from martxclib import martxcbgdtools as bgtools

from time import gmtime, strftime


ver = '0.0.1'
today = strftime("%Y%m%d", gmtime())



parser = Parser(
    description=__doc__,
    epilog="""Chien-Ting Chen <ct.chen@nasa.gov>""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parser.add_argument(
    '-datapath', type=str, required=False, default='/home/ctchen/art-xc',
    help='ART-XC data path')

parser.add_argument(
    '-imgpath', type=str, required=False, default='/home/ctchen/art-xc/workspace/imaging',
    help='ART-XC data path')

parser.add_argument(
    '-exppath', type=str, required=False, default='/home/ctchen/art-xc/workspace/expmap',
    help='ART-XC data path')

parser.add_argument(
    '-evt', type=str, required=False,
    help='Name of the event file.')

parser.add_argument(
    '-att', type=str, required=False,
    help='Name of the attitude file (with quaternions).')

parser.add_argument(
    '-exp', type=str, required=False, 
    help='Name of the exp file.')

parser.add_argument(
    '-makeexp', type=str2bool, required=False, default=True,
    help='Set to False to skip expmap making'
    )
parser.add_argument(
    '-img', type=str, required=False, 
    help='Name of the image file.')

# optional arguments
parser.add_argument('-usexy', type=str2bool, required=False, default=True,
    help='Make an image using the X/Y column from the event list')

parser.add_argument('-raw', type=bool, required=False, default=False,
    help='Using RAW_X and RAW_Y if set as True.')

parser.add_argument('-gyrocoords', type=bool, required=False, default=False,
    help='Using gyro RA/DEC as opposed to event list original RA/DEC if set as True.')

'''
parser.add_argument(
    '-rasize', type=float, required=False, default=15.,
    help='RA pixel size in arcsec.')

parser.add_argument(
    '-decsize', type=float, required=False, default=15.,
    help='DEC pixel size in arcsec.')
'''
parser.add_argument(
    '-projection', type=str, required=False, default='TAN',
    help='DEC pixel size in arcsec.')

parser.add_argument(
    '-box', nargs='+', type=float, default=[],
    help="""list of ra dec values specifying the 4 corners of the box.
    Example: -box ra1 ra2 dec1 dec2""")

parser.add_argument('-overwrite', type=str2bool, required=False, default=True,
    help='Overwrite if set as True.')

parser.add_argument(
    '-tilelist', nargs='+', type=str, default=None,
    help="""list of tiles to be processed""")
parser.add_argument(
    '-tele', nargs='+', type=int, default=None,
    help="""list of telescopes to be used.""")
parser.add_argument('-verbose', type=str2bool, required=False, default=True,
    help='A verbosity parameter')


args = parser.parse_args()
verbose = args.verbose
vprint = verboseprint(verbose)
vprint('arguments parsed:\n')
vprint(args)

datapath = os.path.abspath(args.datapath)
imgpath = os.path.abspath(args.imgpath)
exppath = os.path.abspath(args.exppath)
tilelist = args.tilelist
tele = args.tele
evt = args.evt
att = args.att
img = args.img
exp = args.exp
makeexp = args.makeexp
overwrite=args.overwrite
att = args.att
# get the file list for making images
vprint('overwrite = ', overwrite)

vprint('getting all the file names')
filer = martFileHandler(surveypath=datapath,level=2,tile=tilelist, tele=tele, file=evt)
filer.GetTiles()
filer.GetTeles()
# check tile paths first
if evt is None:
    for t in filer.tiles:
        if not os.path.exists(imgpath + '/' + t + '/'):
            os.system('mkdir '+ imgpath + '/' + t)
        if not os.path.exists(exppath + '/' + t + '/'):
            os.system('mkdir '+ exppath + '/' + t)

for i, file in enumerate(filer.allfiles):
    imager = martImager(file,datapath=datapath)
    if evt is None:
        imgname = imgpath + '/' + filer.alltiles[i] + '/img_' + imager.evtlister.lv1date + '_' + imager.evtlister.lv1time + '.' + filer.alltiles[i] + '.' + filer.allteles[i] + '.fits' 
        imgname2 = imgpath + '/' + filer.alltiles[i] + '/img_corner_' + imager.evtlister.lv1date + '_' + imager.evtlister.lv1time + '.' + filer.alltiles[i] + '.' + filer.allteles[i] + '.fits' 
        expname = exppath + '/' + filer.alltiles[i] + '/expfov_' + imager.evtlister.lv1date + '_' + imager.evtlister.lv1time + '.' + filer.alltiles[i] + '.' + filer.allteles[i] + '.fits'
        expname2 = exppath + '/' + filer.alltiles[i] + '/expcorner_' + imager.evtlister.lv1date + '_' + imager.evtlister.lv1time + '.' + filer.alltiles[i] + '.' + filer.allteles[i] + '.fits'
    if evt is not None and img is None:
        imgname = imgpath + '/img.fits'
        imgname2 = imgpath + '/img_corner.fits'
        print('image name not set (by -img), saving them as ' + imgpath + '/img.fits and img_corner.fits')
    elif evt is not None and img is not None:
        imgname = img
        imgname2 = imgname[:-5] + '.corner.fits'
    if evt is not None and exp is None:
        expname = exppath + '/exp.fits'
        expname2 = exppath + '/exp_corner.fits'
        print('exp name not set (by -img), saving them as ' + exppath + '/exp.fits and exp_corner.fits')
    elif evt is not None and exp is not None:
        expname = exp
        expname2 = expname[:-5] + '.corner.fits'
    vprint('Flare filtering...')
    flaretool = bgtools.martFlareTool(file,time_step = 10.)
    clean_data = flaretool.FilterData_FlareTimes(flaretool.evt_data_original) 
    if len(clean_data) == 0:
        clean_data = flaretool.evt_data_original
        print('Flare filtering results in an empty list, assuming no flares...')
    photon_data = imager.FilterData_Energy(clean_data,4,30) 
    idata_fov = imager.FilterData_Flag(photon_data, 0)
    idata_corner = imager.FilterData_Flag(photon_data, 2)
    ihdu_fov = imager.Make_Image_SkyXY(idata_fov)
    ihdu_corner = imager.Make_Image_SkyXY(idata_corner)
    if (overwrite == False):
        if os.path.exists(imgname):
            vprint(imgname + ' already exists, skpping')
        else:
            ihdu_fov.writeto(imgname,overwrite=overwrite)
        if os.path.exists(imgname2):
            vprint(imgname2 + ' already exists, skpping')
        else:
            ihdu_corner.writeto(imgname2,overwrite=overwrite)
    else:
        vprint('saving image as ' + imgname)
        ihdu_fov.writeto(imgname,overwrite=overwrite)
        ihdu_corner.writeto(imgname2,overwrite=overwrite)
        vprint('saving corner image as ' + imgname2)
    if makeexp:
        if overwrite == False:
            if os.path.exists(expname):
                print(expname + ' already exists, skpping')
            else:
                vprint('saving exp as ' + expname)
                exphdu = imager.Make_Expmap(attname=att, vig=True, flag=0, imagehdu=ihdu_fov)
                exphdu.writeto(expname,overwrite=overwrite)
            if os.path.exists(expname2):
                print(expname2 + ' already exists, skpping')
            else:
                vprint('saving corner expmap as ' + expname2)
                exphdu = imager.Make_Expmap(attname=att, vig=False, flag=2, imagehdu=ihdu_corner)
                exphdu.writeto(expname2,overwrite=overwrite)
        else:
            vprint('saving exp as ' + expname)
            exphdu = imager.Make_Expmap(attname=att, vig=True, flag=0, imagehdu=ihdu_fov)
            exphdu.writeto(expname,overwrite=overwrite)
            vprint('saving corner expmap as ' + expname2)
            exphdu = imager.Make_Expmap(attname=att, vig=False, flag=2, imagehdu=ihdu_corner)
            exphdu.writeto(expname2,overwrite=overwrite)
