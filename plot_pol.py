#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from astropy.io import fits

PRO_CATG_I = 'S2D_POL_I'
PRO_CATG_STOKES = 'S2D_POL_STOKES'
PRO_CATG_NULL = 'S2D_POL_NULL'
PRO_CATG_IN = 'S2D_A'

def main(prefix):
    W = fits.open(f'{prefix}_{PRO_CATG_I}.fits')['WAVEDATA_VAC_BARY'].data
    I = fits.open(f'{prefix}_{PRO_CATG_I}.fits')['SCIDATA'].data
    Ie = fits.open(f'{prefix}_{PRO_CATG_I}.fits')['ERRDATA'].data
    S = fits.open(f'{prefix}_{PRO_CATG_STOKES}.fits')['SCIDATA'].data
    Se = fits.open(f'{prefix}_{PRO_CATG_STOKES}.fits')['ERRDATA'].data
    if os.path.exists(f'{prefix}_{PRO_CATG_NULL}.fits'):
        N = fits.open(f'{prefix}_{PRO_CATG_NULL}.fits')['SCIDATA'].data
        Ne = fits.open(f'{prefix}_{PRO_CATG_NULL}.fits')['ERRDATA'].data
    else:
        N = np.zeros_like(I)
        Ne = np.zeros_like(N)
    
    fig,axs = plt.subplots(3,1, sharex=True)
    plt.suptitle(f'{prefix}')
    for order in range(I.shape[0]):
        axs[0].plot(W[order],I[order],'-')
        axs[0].fill_between(W[order], I[order]-Ie[order], I[order]+Ie[order], alpha=0.1)
        axs[0].set_ylabel('I')
        axs[1].plot(W[order],S[order],'-')
        axs[1].fill_between(W[order], S[order]-Se[order], S[order]+Se[order], alpha=0.1)
        axs[1].set_ylabel('Stokes/I')
        axs[2].plot(W[order],N[order],'-')
        axs[2].fill_between(W[order], N[order]-Ne[order], N[order]+Ne[order], alpha=0.1)
        axs[2].set_ylabel('Null')
        axs[2].set_xlabel('Wavelength')
    plt.show()
    
if __name__ == "__main__":
    prefix = sys.argv[1]
    main(prefix)