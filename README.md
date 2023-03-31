# MUVIT (Mock UV-data Injector Tool)

## Introduction


MUVIT is a python code that allows to inject mock visibilities into a real interferometric radio observation. MUVIT was developed to generate mock radio halos, but can be tailored to simulate point sources and many other extended radio sources. 


The current release of the code:

1. includes 2D-exponential and 2D-gaussian profiles to model the surface brightness of the mock emission. 
2. was tested on LOFAR, GMRT/uGMRT, VLA/JVLA, VLBA data


Questions, suggestions and issues can be reported by sending an email to: luca.bruno4@unibo.it  


### Citation 

If you make use of MUVIT, please cite the following paper: "The Planck clusters in the LOFAR sky. II. LoTSS-DR2: Recovering diffuse extended emission with LOFAR", Bruno et al. 2023, A&A, 672, A41 

https://ui.adsabs.harvard.edu/search/fq=%7B!type%3Daqp%20v%3D%24fq_database%7D&fq_database=(database%3Aastronomy)&q=%20%20author%3A%22%5Ebruno%22%20year%3A2023&sort=date%20desc%2C%20bibcode%20desc&p_=0


### Requirments

MUVIT is written in python3 and makes use of WSClean (https://gitlab.com/aroffringa/wsclean) for imaging and Fourier-transforms. The main required python packages are:

1. NUMPY (https://numpy.org)
2. ASTROPY (https://github.com/astropy/)
3. PYLAB (https://www.tutorialspoint.com/matplotlib/matplotlib_pylab_module.htm)
4. ARGPARSE (https://docs.python.org/3/library/argparse.html)
5. CASACORE (https://github.com/casacore/casacore)


### Usage

Before running the code (STEP 0): image your dataset in WSClean at FULL RESOLUTION (no uvtaper in this step) (please specify "-no-update-model-required"). The imaging parameters will be automatically exploited by MUVIT.  


Example of how to run the code (see details with python3 MUVIT.py --help):
python3 MUVIT.py --RA 328.5 --DEC 17.67 --flux 50.5 --z 0.233 --re 200 --input_fits TEST-image.fits --ms_files *.calibrated

Required parameters are:
1. coordinates (in deg) of the centre of the injection
2. total injected flux density (in mJy)
3. redshift
4. e-folding radius (exponential model) or standard deviation (Gaussian model) (in kpc)
5. image created at STEP 0
6. measurment set(s) 

By default, MUVIT will produce an image of the mock emission with the same imaging parameters of STEP 0. Tapering of the baselines can be done (see help). 

WARNING: MUVIT updates the MODEL_DATA column with the mock+real visibilities. At each running, the MODEL_DATA column is overwritten!


