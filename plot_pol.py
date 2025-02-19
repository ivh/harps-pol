#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import os
import sys
from astropy.io import fits
from pathlib import Path

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
        
    return locals()

def plot_oldreduc(dir, prefix, W, axs, alpha=0.5):
    oldred = dir/'../oldreduc'
    oI=fits.open(oldred/(pref+f'_I.fits'))[0].data
    if oI.shape[0]==71:
        oI=np.delete(oI, 44, axis=0)
    for order in range(oI.shape[0]):
        axs[0].plot(W[order],oI[order],'k',alpha=alpha)
    oN=fits.open(oldred/(pref+f'_Null.fits'))[0].data
    if oN.shape[0]==71:
        oN=np.delete(oN, 44, axis=0)
    for order in range(oN.shape[0]):
        axs[2].plot(W[order],oN[order],'k',alpha=alpha)
    oS=fits.open(oldred/(pref+f'_v.fits'))[0].data
    if oS.shape[0]==71:
        oS=np.delete(oS, 44, axis=0)
    for order in range(oS.shape[0]):
        axs[1].plot(W[order],oS[order],'k',alpha=alpha)

def plot_alexis(dir, prefix, W, axs, alpha=0.5,scale=0.28):
    d=fits.open(dir/'../alexisred'/(pref+f'_demodulated.fits'))[1].data
    axs[0].plot(d['WAVE']*10,d['I']*scale,'k',alpha=alpha)
    axs[1].plot(d['WAVE']*10,d['Stokes'],'k',alpha=alpha)
    axs[2].plot(d['WAVE']*10,d['Null'],'k',alpha=alpha)
    return locals()
    
if __name__ == "__main__":
    prefix = sys.argv[1]
    locals().update(main(prefix))
    
    dir = Path(prefix).parent
    pref = Path(prefix).name
    pref=pref.removeprefix('r.')

    #locals().update(plot_oldreduc(dir, pref, W,axs, alpha=0.5))
        
    locals().update(plot_alexis(dir, pref, W,axs, alpha=0.5))

    axs[0].set(xlim=(5459.5,5481))
    axs[0].set(ylim=(30000,50000))
    axs[1].set(ylim=(-0.025, 0.025))
    axs[2].set(ylim=(-0.02622, 0.01409))
