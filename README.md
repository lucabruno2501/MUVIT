# MUVIT (Mock UV-data Injector Tool)

## Introduction


MUVIT is a python code that allows to inject mock visibilities into a real interferometric radio observation. MUVIT was developed to generate mock radio halos, but can be tailored to simulate point sources and many other extended radio sources. 


The current release of the code:

1. includes 2D-exponential and 2D-gaussian profiles to model the surface brightness of the mock emission. 
2. was tested on LOFAR, GMRT/uGMRT, VLA/JVLA, VLBA data

Future releases will include further models.

Questions, suggestions and issues can be reported by sending an email to: luca.bruno4@unibo.it  


## Citation 

If you make use of MUVIT, please cite the following paper: "The Planck clusters in the LOFAR sky. II. LoTSS-DR2: Recovering diffuse extended emission with LOFAR", Bruno et al. 2023, A&A, 672, A41 

https://www.aanda.org/articles/aa/full_html/2023/04/aa44552-22/aa44552-22.html



## Requirements

MUVIT is written in python3 and makes use of WSClean (https://gitlab.com/aroffringa/wsclean) for imaging and Fourier-transforms. The main required python packages are:

1. NUMPY (https://numpy.org)
2. ASTROPY (https://github.com/astropy/)
3. PYLAB (https://www.tutorialspoint.com/matplotlib/matplotlib_pylab_module.htm)
4. ARGPARSE (https://docs.python.org/3/library/argparse.html)
5. CASACORE (https://github.com/casacore/casacore)


## Usage

STEP0: before running the code, it is necessary to image your dataset in WSClean. For this step, do not apply uv-tapering and specify the parameter "-no-update-model-required". Your imaging parameters will be then automatically exploited by MUVIT to image the mock visibilities.   

WARNING: MUVIT updates the MODEL_DATA column with the mock+real visibilities. At each running, the MODEL_DATA column is overwritten!

USAGE EXAMPLE: python3 MUVIT.py --RA 328.5 --DEC 17.67 --flux 50.5 --z 0.233 --re 200 --input_fits imagefromstep0.fits --ms_files *.ms



### REQUIRED ARGUMENTS:

  --input_fits INPUT_FITS  Image from STEP0
                        
  --RA RA               Right Ascension (deg)
  
  --DEC DEC             Declination (deg)
  
  --z Z                 Redshift
  
  --re RE               e-folding radius (kpc)
  
  --flux FLUX           Total injected flux density at the reference frequency (mJy)
  
  --ms_files MS_FILES [MS_FILES ...]   Measurement sets


### OPTIONAL ARGUMENTS:

  -h, --help            show this help message and exit
  
  --model MODEL         exponential (EXP) or Gaussian (GAUSS); default EXP
  
  --name NAME           Object name; default TARGET
  
  --taper_kpc TAPER_KPC
                        Tapering (in kpc)
                        
  --taper_arcsec TAPER_ARCSEC
                        Tapering (in arcsec)
                        
  --spix SPIX           Spectral index; default -1.3
  
  --do_0inj DO_0INJ     Do only imaging without injection and exit; default False



### Notes on arguments

1. By default, MUVIT will produce an image of the mock emission with the same imaging parameters of STEP0. Tapering of the baselines can be activated with --taper_kpc or --taper_arcsec to enhance faint extended emission
2. --do_0inj True allows to produce a (tapered) image of the data without performing any injections. 
3. For --model GAUSS, the e-folding radius is assumed to be the standard deviation of the Gaussian profile (see Bruno et al. 2023 for details). 
4. --spix allows to change the default spectral index of the mock source (for multi-frequency models)





