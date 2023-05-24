import numpy as np
import argparse


"""
ARGPARSE
"""

parser = argparse.ArgumentParser(description = 'Estimate flux density of upper limit for LOFAR HBA dataset (see Eq. 5 in Bruno et al. 2023)')

parser.add_argument('--beam', type=float, nargs=2, help = 'Restoring beam axes (arcsec, arcsec)', required=True)
parser.add_argument('--D', type=float, help = 'diameter of mock halo (arcsec)', required=True)
parser.add_argument('--rms', type=float, help = 'Noise of the image (mJy/beam)', required=True)
parser.add_argument('--spix', type=float, help = 'Spectral index to re-scale flux density at 150 MHz', required=False, default=-1.3)

args = parser.parse_args()

D = args.D
beam = args.beam
sigma = args.rms
alpha = args.spix

m = 0.789
q = -0.091


Nbeam = D**2 / (beam[0]*beam[1])
print(Nbeam)
S_UL = 10**(q) * sigma * Nbeam**(m)

S_UL = round(S_UL,2)
print('Number of beams: ', Nbeam)
print('Expected upper limit at 144 MHz: ', S_UL, ' mJy')


S_UL_150 = S_UL * (150/144)**alpha
S_UL_150 = round(S_UL_150,2)
print('Expected upper limit at 150 MHz: ', S_UL_150, ' mJy')
    



