#!/usr/bin/env python

#python3 MUVIT.py --RA 328.5 --DEC 17.67 --flux 50.5 --ms_files *.calibrated --re 200 --z 0.233 --input_fits TEST-image.fits

#Mock UV-data Injector Tool (MUVIT)



from astropy.io import fits
import argparse
from astropy import units as u
import casacore.tables as tables
import glob
import os
import sys
import numpy as np
from pylab import meshgrid, imshow
from astropy.wcs import WCS
from astropy.wcs import WCS as pywcs
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


"""
ARGPARSE
"""

parser = argparse.ArgumentParser(description = 'Perform injection of mock visibilities in uv-dataset')



parser.add_argument('--name', type=str, help = 'Object name; default TARGET', default = 'TARGET', required=False)

parser.add_argument('--input_fits', type=str, help = 'Full resolution image from a previous wsclean run', required=True)

parser.add_argument('--RA', type=float, help = 'Right Ascension (deg)', required=True)

parser.add_argument('--DEC', type=float, help = 'Declination (deg)', required=True)

parser.add_argument('--z', type=float, help = 'Redshift', required=True)

parser.add_argument('--re', type=float, help = 'e-folding radius (kpc)', required=True)
parser.add_argument('--flux', type=float, help = 'Total injected flux density at the reference frequency (mJy)', required=True)

parser.add_argument('--model', type=str, help = 'exponential (EXP) or Gaussian (GAUSS); default EXP', default = 'EXP', required=False)

parser.add_argument('--ms_files', help = 'Measurement sets', nargs='+', required=True)

parser.add_argument('--taper_kpc', type=float, help = 'Tapering (in kpc)', required=False)
parser.add_argument('--taper_arcsec', type=float, help = 'Tapering (in arcsec)', required=False)
parser.add_argument('--spix', type=float, help = 'Spectral index; default -1.3', default = -1.3, required=False)
parser.add_argument('--do_0inj', type=bool, help = 'Do only imaging without injection and exit; default False', default = False)



args = parser.parse_args()



"""
INPUTS
"""

#input_image = 'TEST.fits'
input_image = args.input_fits
RA = args.RA
DEC = args.DEC
re = args.re
flux = 1.e-3 * args.flux
ms_files = args.ms_files
ms_files_wsclean = ' '.join(ms_files)



root_name = args.name
root_img = root_name
z = args.z
input_model = args.model



spidx = args.spix
do_0inj = args.do_0inj
taper_kpc = args.taper_kpc
taper_arcsec = args.taper_arcsec


"""
ANGULAR/LINEAR CONVERSION
"""

#arcsec/kpc factor
arcsec_to_kpc = cosmo.arcsec_per_kpc_proper(z).value
#kpc/arcsec factor
kpc_to_arcsec=1./arcsec_to_kpc


#e-folding radius in arcsec
re_as = re*arcsec_to_kpc


  
def find_data_column(text, word_to_search):
    text_array = text.split()
    try:
        for word in text_array:
            if word == word_to_search:
                index = text_array.index(word) + 1
                datacolumn_0 = text_array[index]
        return text, datacolumn_0
    except:
        datacolumn_0 = 'DATA'
        text_array.insert(1, '-data-column')
        text_array.insert(2, datacolumn_0)
        text = ' '.join(text_array)
        return text, datacolumn_0
        



data_fits = fits.open(input_image)
header_fits = data_fits[0].header


print(header_fits['HISTORY'])


#wsclean history is output as array. this raises issues from spaces and need to be properly modified before converting in string
wsclean_array = header_fits['HISTORY']
wsclean_command = []
for i in range(len(wsclean_array)):
    if len(wsclean_array[i]) == 71:
        wsclean_command.append(wsclean_array[i] + ' ')
    else:
         wsclean_command.append(wsclean_array[i])

wsclean_command = ''.join(wsclean_command)
wsclean_command = find_data_column(wsclean_command, '-data-column')[0]




def find_clean_parameter(text, word_to_search):
    text_array = text.split()
    try:
        for word in text_array:
            if word == word_to_search:
                index = text_array.index(word) + 1
                clean_parameter = text_array[index]
        return clean_parameter
    except:
        clean_parameter = False
        return clean_parameter


datacolumn_0 = find_data_column(wsclean_command, '-data-column')[1]
imsize = int(header_fits['NAXIS1'])
ref_freq = float(header_fits['CRVAL3'])
pixscale = round(float(abs((header_fits['CDELT1']*3600))),2)



nchan = find_clean_parameter(wsclean_command, '-channels-out')
name_fullres0 = find_clean_parameter(wsclean_command, '-name')


def substitute_clean_parameter(text, word_to_search, word_new, narguments):
    text_array = text.split()
    for word in text_array:
        if word == word_to_search:
            index = text_array.index(word) + narguments
            text_array[index] = word_new
    return ' '.join(text_array)


def add_clean_parameter(text, clean_option, clean_parameter):
    text_array = text.split()
    text_array.insert(1, clean_option)
    text_array.insert(2, clean_parameter)
    
    text = ' '.join(text_array)
    return text


def add_logfile(text, logfile):
    text_array = text.split()
    text_array.append(logfile)
    text = ' '.join(text_array)
    return text

"""
TAPER
"""

taper_to_use = None
if taper_arcsec != None and taper_kpc != None:
    print('Error: please select only one taper parameter as input')
    sys.exit()
elif taper_kpc != None:
    taper_to_use = int(taper_kpc/kpc_to_arcsec)
elif taper_arcsec != None:
    taper_to_use = int(taper_arcsec)

if taper_to_use != None:
    pixscale_tapered = int(taper_to_use/5.)
    if int(pixscale_tapered) == 1:
        pixscale_tapered = 2
    imsize_tapered = int(imsize*pixscale/pixscale_tapered)
    



"""
FLATTEN
"""


def flatten(filename, channel=0, freqaxis=0):
    f = fits.open(filename)
    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Cannot make map from this')
    if naxis==2:
        pass

    w = pywcs(f[0].header)
    wn = pywcs(naxis=2)

    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]

    header = wn.to_header()
    header['NAXIS']=2
    header['NAXIS1']=f[0].header['NAXIS1']
    header['NAXIS2']=f[0].header['NAXIS2']
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    dataslice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            dataslice.append(np.s_[:],)
        elif i==freqaxis:
            dataslice.append(channel)
        else:
            dataslice.append(0)
# add freq
    header["FREQ"] = f[0].header['CRVAL3']

    # add beam if present
    try:
        header["BMAJ"]=f[0].header['BMAJ']
        header["BMIN"]=f[0].header['BMIN']
        header["BPA"]=f[0].header['BPA']
    except:
        pass
    return header, f[0].data[tuple(dataslice)]


"""
LOW RESOLUTION IMAGE WITHOUT INJECTION
"""

if do_0inj == True:
            wsclean0 = substitute_clean_parameter(wsclean_command, '-size', str(imsize_tapered), 1)
            wsclean0 = substitute_clean_parameter(wsclean0, '-size', str(imsize_tapered), 2)
            wsclean0 = substitute_clean_parameter(wsclean0, '-scale', str(pixscale_tapered)+'arcsec', 1)
            wsclean0 = substitute_clean_parameter(wsclean0, '-name', root_img+'_ORIGINAL_T'+str(taper_to_use)+'arcsec', 1)
            wsclean0 = add_clean_parameter(wsclean0, '-taper-gaussian', str(taper_to_use)+'arcsec')
            wsclean0 = add_logfile(wsclean0, '>wsclean0.log')
            print('********************************************')
            print('Running wsclean without injection')
            print('********************************************')
            print()
            print(wsclean0)
            os.system(wsclean0)
            os.system('rm -rf *ORIGINAL_T*000*.fits') #channel images
            os.system('rm -rf *ORIGINAL_T*-dirty.fits') #dirty images
            os.system('rm -rf *ORIGINAL_T*-psf.fits') #psf images
            sys.exit()
            




##############################################################
#   CREATE MODELS with niter 1
##############################################################


wsclean_niter1 = substitute_clean_parameter(wsclean_command, '-name', root_img+'_mod'+str(flux*1000)+'mJy', 1)
wsclean_niter1 = substitute_clean_parameter(wsclean_niter1, '-niter', str(1), 1)
wsclean_niter1 = add_logfile(wsclean_niter1, '>wsclean_niter1.log')
os.system(wsclean_niter1)

if nchan == False:
    models = glob.glob(root_img+'_mod'+str(flux*1000)+'mJy-model.fits')
else:
    models = glob.glob(root_img+'_mod'+str(flux*1000)+'mJy-0*-model.fits')


##############################################################
#   UPDATE MODELS: add mock halo
##############################################################


head, data = flatten(models[0])
wcs = WCS(head)
center_ra_px, center_dec_px = wcs.wcs_world2pix(RA*u.deg, DEC*u.deg, 1, ra_dec_order=True)

x = np.linspace(0, imsize, imsize)
X,Y = meshgrid(x, x)


##############################################################
#   MODEL WITH EXPONENTIAL FUNCTION
##############################################################

def exponential_2D(flux, spidx, freq, re_as, center_ra, center_dec, x, y):
    re_px = re_as/pixscale
    flux_freq = flux*(freq/ref_freq)**spidx
    I_0 = flux_freq/(2.*np.pi*re_as**2) #Jy/arcsec^2
    r = np.sqrt((x-center_ra)**2+(y-center_dec)**2)
    exponential = I_0 * np.exp(- r/re_px) * pixscale**2 #Jy/pixel
    

    return exponential


##############################################################
#   MODEL WITH GAUSSIAN FUNCTION (sigma = re)
##############################################################

def gauss_2D(flux, spidx, freq, re_as, center_ra, center_dec, x, y):
    re_px = re_as/pixscale
    flux_freq = flux*(freq/ref_freq)**spidx
    I_0 = flux_freq/(2.*np.pi*re_as**2) #Jy/arcsec^2
    r=np.sqrt((x-center_ra)**2+(y-center_dec)**2)
    gaussian = I_0 * np.exp(- 0.5*(r/re_px)**2) * pixscale**2

    return gaussian


##############################################################
#   IMAGES of MOCK EMISSION
##############################################################

print('********************************************')
print('Creating model images of mock emission')
print('********************************************')
print()

for model in models:
    hdr = fits.getheader(model)
    freq = hdr['CRVAL3']
    if input_model == 'EXP':
        mock_halo = exponential_2D(flux, spidx, freq, re_as, center_ra_px, center_dec_px, X, Y)
    elif input_model == 'GAUSS':
        mock_halo = gauss_2D(flux, spidx, freq, re_as, center_ra_px, center_dec_px, X, Y)
    model_update = 'inject_'+model
    os.system('cp '+model+' '+model_update)

    fits.update(model_update, mock_halo, hdr)
    


##############################################################
#   FOURIER TRANSFORM: obtain mock visibilities
##############################################################

print('********************************************')
print('Predicting mock visibilities')
print('********************************************')
print()

if nchan != False:
    os.system('wsclean -predict -name inject_'+root_img+'_mod'+str(flux*1000)+'mJy -channels-out ' +str(nchan)+' '+ms_files_wsclean+' >wsclean_predict.log')
if nchan == False:
    os.system('wsclean -predict -name inject_'+root_img+'_mod'+str(flux*1000)+'mJy '+ms_files_wsclean+' >wsclean_predict.log')


##############################################################
#   ADD MOCK VISIBILITIES to the data
##############################################################

print('********************************************')
print('Updating MODEL_DATA with mock+real data')
print('********************************************')
print()



for ms_file in ms_files:
    taql_command = 'update "'+ms_file+'" set MODEL_DATA = '+datacolumn_0+' + MODEL_DATA'
    tables.taql(taql_command)
    
    

##############################################################
#   IMAGE THE REAL+MOCK DATA
##############################################################


if taper_to_use == None:
    wsclean_final = substitute_clean_parameter(wsclean_command, '-name', root_img+'_MOCK_'+str(flux*1000)+'mJy', 1)
if taper_to_use != None:
    wsclean_final = substitute_clean_parameter(wsclean_command, '-name', root_img+'_MOCK_'+str(flux*1000)+'mJy_T'+str(taper_to_use)+'arcsec', 1)
    wsclean_final = substitute_clean_parameter(wsclean_final, '-size', str(imsize_tapered), 1)
    wsclean_final = substitute_clean_parameter(wsclean_final, '-size', str(imsize_tapered), 2)
    wsclean_final = substitute_clean_parameter(wsclean_final, '-scale', str(pixscale_tapered)+'arcsec', 1)
    wsclean_final = add_clean_parameter(wsclean_final, '-taper-gaussian', str(taper_to_use)+'arcsec')
   
   

wsclean_final = substitute_clean_parameter(wsclean_final, '-data-column', 'MODEL_DATA', 1)
wsclean_final = add_logfile(wsclean_final, '>wsclean_final.log')




print('********************************************')
print('Imaging mock+real data')
print('********************************************')
print()
print(wsclean_final)
os.system(wsclean_final)


##############################################################
#   REMOVE USELESS FILES
##############################################################


os.system('rm -rf '+str(root_img)+'_mod'+str(flux*1000)+'*000*') #niter 1 channel images
os.system('rm -rf '+str(root_img)+'_mod'+str(flux*1000)+'*MFS*') #niter 1 MFS images
os.system('rm -rf inject_'+str(root_img)+'_mod'+str(flux*1000)+'*000*') #injected model images
os.system('rm -rf '+str(root_img)+'_MOCK_'+str(flux*1000)+'*T*000*') #mock image channel maps
os.system('rm -rf '+str(root_img)+'_MOCK_'+str(flux*1000)+'*T*psf.fits') #mock image psf map
os.system('rm -rf '+str(root_img)+'_MOCK_'+str(flux*1000)+'*T*dirty.fits') #mock image dirty map
