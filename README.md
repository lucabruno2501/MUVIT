# MUVIT
Mock UV-data Injector Tool

MUVIT is a python code that allows to inject mock visibilities into a real dataset of a radio observation.
Requirments: python3 and WSClean


MUVIT was developed to produce images of mock radio halos, but can be tailored to many other radio sources. The current release of the code:
a. includes 2D-exponential and 2D-gaussian profiles to model the surface brightness of the mock emission. 
b. was tested on LOFAR, GMRT/uGMRT, VLA/JVLA data


USAGE:

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


