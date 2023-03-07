# MUVIT
Mock UV-data Injector Tool

MUVIT is a python code that allows to inject mock visibilities into a real dataset of a radio observation. 

Before running the code (STEP 0): image your dataset in WSClean at FULL RESOLUTION (no uvtaper in this step) (please specify "-no-update-model-required" --> to check if this is default) 

Example of how to run the code:
python3 MUVIT.py --RA 328.5 --DEC 17.67 --flux 50.5 --ms_files *.calibrated --re 200 --z 0.233 --input_fits TEST-image.fits

RA, DEC coordinate of the centre of the injection (in deg) 

flux total injected flux

*.calibrated all the measurement sets 

re e-folding radius in kpc

z redshift

input_fits is the image created at STEP 0


Check for further option with python3 MUVIT.py --help


