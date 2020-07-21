from martxclib.preprocessfunc import *
import sys
import numpy as np
import astropy.io.fits as fits
from astropy import wcs
import astropy.coordinates as cd
import astropy.units as u
from martxclib.martxcfun import *
from martxclib.martxcfits import *
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
    '-out', type=str, required=False, default='expmap.fits',
    help='Name of the output file.')

# optional arguments
parser.add_argument('-usexy', type=str2bool, required=False, default=True,
    help='Make an image using the X/Y column from the event list')

parser.add_argument('-raw', type=bool, required=False, default=False,
    help='Using RAW_X and RAW_Y if set as True.')

parser.add_argument('-gyrocoords', type=bool, required=False, default=False,
    help='Using gyro RA/DEC as opposed to event list original RA/DEC if set as True.')

parser.add_argument(
    '-rasize', type=float, required=False, default=15.,
    help='RA pixel size in arcsec.')

parser.add_argument(
    '-decsize', type=float, required=False, default=15.,
    help='DEC pixel size in arcsec.')

parser.add_argument(
    '-projection', type=str, required=False, default='TAN',
    help='DEC pixel size in arcsec.')

parser.add_argument(
    '-box', nargs='+', type=float, default=[],
    help="""list of ra dec values specifying the 4 corners of the box.
    Example: -box ra1 ra2 dec1 dec2""")

parser.add_argument('-verbose', type=bool, required=False, default=True,
    help='Add a time constraint to the attitude file. Only the periods within the specified time would be considered.')

parser.add_argument('-overwrite', type=bool, required=False, default=True,
    help='Overwrite if set as True.')

parser.add_argument(
    '-tilelist', nargs='+', type=str, default=None,
    help="""list of tiles to be processed""")
parser.add_argument(
    '-tele', nargs='+', type=int, default=None,
    help="""list of telescopes to be used.""")

args = parser.parse_args()

datapath = os.path.abspath(args.datapath)
imgpath = os.path.abspath(args.imgpath)
exppath = os.path.abspath(args.exppath)
tilelist = args.tilelist
tele = args.tele
evt = args.evt
out = args.out
overwrite=args.overwrite


# get the file list for making images
filer  = martFileHandler(surveypath=datapath,level=2,tile=tilelist, tele=tele, file=evt)
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
    flaretool = bgtools.martFlareTool(file,time_step = 0.1)
    clean_data = flaretool.FilterData_FlareTimes(flaretool.evt_data_original) 
    if len(clean_data) == 0:
        clean_data = flaretool.evt_data_original
        print('Flare filtering results in an empty list, assuming no flares...')
    photon_data = imager.FilterData_Energy(clean_data,4,30) 
    idata_fov = imager.FilterData_Flag(photon_data, 0)
    idata_corner = imager.FilterData_Flag(photon_data, 2)
    ihdu_fov = imager.Make_Image_SkyXY(idata_fov)
    ihdu_corner = imager.Make_Image_SkyXY(idata_fov)


    imgname = imgpath + '/' + filer.alltiles[i] + '/img_' + imager.evtlister.lv1date + '_' + imager.evtlister.lv1time + '.' + filer.alltiles[i] + '.' + filer.allteles[i] + '.fits' 
    ihdu_fov.writeto(imgname,overwrite=overwrite)
    imgname = imgpath + '/' + filer.alltiles[i] + '/img_corner_' + imager.evtlister.lv1date + '_' + imager.evtlister.lv1time + '.' + filer.alltiles[i] + '.' + filer.allteles[i] + '.fits' 
    ihdu_corner.writeto(imgname,overwrite=overwrite)
    exphdu = imager.Make_Expmap(vig=True, flag=0, imagehdu=ihdu_fov)
    expname = exppath + '/' + filer.alltiles[i] + '/expfov_' + imager.evtlister.lv1date + '_' + imager.evtlister.lv1time + '.' + filer.alltiles[i] + '.' + filer.allteles[i] + '.fits'
    exphdu.writeto(expname,overwrite=overwrite)
    exphdu = imager.Make_Expmap(vig=False, flag=2, imagehdu=ihdu_corner)
    expname = exppath + '/' + filer.alltiles[i] + '/expfov_' + imager.evtlister.lv1date + '_' + imager.evtlister.lv1time + '.' + filer.alltiles[i] + '.' + filer.allteles[i] + '.fits'
    exphdu.writeto(expname,overwrite=overwrite)

