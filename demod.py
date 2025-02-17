#!/usr/bin/env python3

import numpy as np
import sys
from astropy.io import fits
from scipy.interpolate import interp1d

PRO_CATG_I = "S2D_POL_I"
PRO_CATG_STOKES = "S2D_POL_STOKES"
PRO_CATG_NULL = "S2D_POL_NULL"
PRO_CATG_INA = "S2D_A"
PRO_CATG_INAB = "S2D_BLAZE_A"
PRO_CATG_INB = "S2D_B"
PRO_CATG_INBB = "S2D_BLAZE_B"


def readspec(filename):
    s = fits.open(filename)["SCIDATA"].data
    w = fits.open(filename)["WAVEDATA_VAC_BARY"].data
    e = fits.open(filename)["ERRDATA"].data

    # Fiber A has an extra order, pop it out.
    nord = s.shape[0]
    if nord == 71:
        s = np.delete(s, 44, axis=0)
        w = np.delete(w, 44, axis=0)
        e = np.delete(e, 44, axis=0)
    elif nord != 70:
        print("Wrong number of orders")
        exit(1)

    # Mask zeros and NaNs
    mask = (s == 0) | np.isnan(s)
    s = np.ma.array(s, mask=mask)
    e = np.ma.array(e, mask=mask)

    return w, s, e


def save(fname, X, Xe, I, Ie, N=None, Ne=None):
    # Convert masked arrays to regular arrays, filling masked values with NaN
    X = X.filled(np.nan) if isinstance(X, np.ma.MaskedArray) else X
    Xe = Xe.filled(np.nan) if isinstance(Xe, np.ma.MaskedArray) else Xe
    I = I.filled(np.nan) if isinstance(I, np.ma.MaskedArray) else I
    Ie = Ie.filled(np.nan) if isinstance(Ie, np.ma.MaskedArray) else Ie
    if N is not None:
        N = N.filled(np.nan) if isinstance(N, np.ma.MaskedArray) else N
        Ne = Ne.filled(np.nan) if isinstance(Ne, np.ma.MaskedArray) else Ne

    if fname.endswith(f"{PRO_CATG_INB}.fits"):
        insuff = PRO_CATG_INB
    elif fname.endswith(f"{PRO_CATG_INBB}.fits"):
        insuff = PRO_CATG_INBB
    else:
        print("Something fishy with filenames or their ordering!")

    outname = fname.replace(f"{insuff}.fits", f"{PRO_CATG_STOKES}.fits")
    hdul = fits.open(fname)
    hdul["SCIDATA"].data = X
    hdul["ERRDATA"].data = Xe
    hdul[0].header["PRO_CATG"] = PRO_CATG_STOKES
    hdul.writeto(outname, overwrite=True)
    hdul.close()

    outname = fname.replace(f"{insuff}.fits", f"{PRO_CATG_I}.fits")
    hdul = fits.open(fname)
    hdul["SCIDATA"].data = I
    hdul["ERRDATA"].data = Ie
    hdul[0].header["PRO_CATG"] = PRO_CATG_I
    hdul.writeto(outname, overwrite=True)
    hdul.close()

    if N is not None:
        outname = fname.replace(f"{insuff}.fits", f"{PRO_CATG_NULL}.fits")
        hdul = fits.open(fname)
        hdul["SCIDATA"].data = N
        hdul["ERRDATA"].data = Ne
        hdul[0].header["PRO_CATG"] = PRO_CATG_NULL
        hdul.writeto(outname, overwrite=True)
        hdul.close()


def check_wl(wls, delta=0.1):
    """Check that the starting wavelengths of all orders in all frames and fibers
    are the same within delta"""
    wlref = wls[0]
    for order in range(wlref.shape[0]):
        try:
            assert np.allclose(
                [wl[order][0] for wl in wls], wlref[order][0], atol=delta
            )
        except AssertionError as e:
            print(f"Order {order} failed:", [wl[order][0] for wl in wls])
            raise e


def demod2(fnames):
    w1a, s1a, e1a = readspec(fnames[0])
    w1b, s1b, e1b = readspec(fnames[1])
    w2a, s2a, e2a = readspec(fnames[2])
    w2b, s2b, e2b = readspec(fnames[3])
    check_wl([w1a, w1b, w2a, w2b])

    I = 0.25 * (s1a + s1b + s2a + s2b)

    Ra = s2a / s1a
    Rb = s1b / s2b
    for i in range(Rb.shape[0]):
        Rb[i] = interp1d(
            w1b[i], Rb[i].filled(np.nan), fill_value="extrapolate", bounds_error=False
        )(w1a[i])
    R = Ra * Rb
    X = np.sqrt(R) - 1
    X /= np.sqrt(R) + 1

    # Error propagation
    # First calculate Ra and Rb errors using quotient rule
    Rae = Ra * np.sqrt((e1b / s1b) ** 2 + (e2b / s2b) ** 2)
    Rbe = Rb * np.sqrt((e2a / s2a) ** 2 + (e1a / s1a) ** 2)

    # Error in R using product rule
    Re = R * np.sqrt((Rae / Ra) ** 2 + (Rbe / Rb) ** 2)

    # Error in X using chain rule for (sqrt(R)-1)/(sqrt(R)+1)
    sqrtR = np.sqrt(R)
    Xe = 2 * Re / (2 * sqrtR * (sqrtR + 1) ** 2)

    # Error in I using sum rule for mean
    Ie = 0.25 * np.sqrt(e1a**2 + e1b**2 + e2a**2 + e2b**2)

    save(fnames[1], X, Xe, I, Ie) # OBS, using *second* file as template,
                            # since fiber B has one less order and the waves
                            # are not re-written and need to match


def error_helper(Re, R):
    dX_dR = 1 / (2 * R**0.75) * 2 / ((R**0.25 + 1) ** 2)
    return abs(dX_dR) * Re


def demod4(fnames):
    w1a, s1a, e1a = readspec(fnames[0])
    w1b, s1b, e1b = readspec(fnames[1])
    w2a, s2a, e2a = readspec(fnames[2])
    w2b, s2b, e2b = readspec(fnames[3])
    w3a, s3a, e3a = readspec(fnames[4])
    w3b, s3b, e3b = readspec(fnames[5])
    w4a, s4a, e4a = readspec(fnames[6])
    w4b, s4b, e4b = readspec(fnames[7])
    check_wl([w1a, w1b, w2a, w2b, w3a, w3b, w4a, w4b])

    I = 0.125 * (s1a + s1b + s2a + s2b + s3a + s3b + s4a + s4b)
    Ie = 0.125 * np.sqrt(
        e1a**2 + e1b**2 + e2a**2 + e2b**2 + e3a**2 + e3b**2 + e4a**2 + e4b**2
    )

    # Demodulate for Stokes parameter
    Ra = s2a / s1a * s4a / s3a
    Rb = s1b / s2b * s3b / s4b
    for i in range(Rb.shape[0]):
        Rb[i] = interp1d(
            w1b[i], Rb[i].filled(np.nan), fill_value="extrapolate", bounds_error=False
        )(w1a[i])
    R = Ra * Rb
    X = np.sqrt(np.sqrt(R)) - 1
    X /= np.sqrt(np.sqrt(R)) + 1

    # Error propagation
    # First calculate Ra and Rb errors using quotient rule
    Rae = Ra * np.sqrt(
        (e1b / s1b) ** 2 + (e2b / s2b) ** 2 + (e4a / s4a) ** 2 + (e3a / s3a) ** 2
    )
    Rbe = Rb * np.sqrt(
        (e2a / s2a) ** 2 + (e1a / s1a) ** 2 + (e3b / s3b) ** 2 + (e4b / s4b) ** 2
    )
    Re = R * np.sqrt((Rae / Ra) ** 2 + (Rbe / Rb) ** 2)
    Xe = error_helper(Re, R)

    # Demodulate for null
    Ra = s2a / s1a * s3a / s4a
    Rb = s1b / s2b * s4b / s3b
    for i in range(Rb.shape[0]):
        Rb[i] = interp1d(
            w1b[i], Rb[i].filled(np.nan), fill_value="extrapolate", bounds_error=False
        )(w1a[i])
    R = Ra * Rb
    N = np.sqrt(np.sqrt(R)) - 1
    N /= np.sqrt(np.sqrt(R)) + 1

    Rae = Ra * np.sqrt(
        (e1b / s1b) ** 2 + (e2b / s2b) ** 2 + (e4a / s4a) ** 2 + (e3a / s3a) ** 2
    )
    Rbe = Rb * np.sqrt(
        (e2a / s2a) ** 2 + (e1a / s1a) ** 2 + (e3b / s3b) ** 2 + (e4b / s4b) ** 2
    )
    Re = R * np.sqrt((Rae / Ra) ** 2 + (Rbe / Rb) ** 2)
    Ne = error_helper(Re, R)

    save(fnames[1], X, Xe, I, Ie, N, Ne) # see comment at end of demod2()


if __name__ == "__main__":
    
    # This assumes filename sorting results in order 1a, 1b, 2a, 2b, 3a, 3b, 4a, 4b
    fnames = sorted(sys.argv[1:])
    print(fnames)

    # 2 or 4 raw frames means 4 or 8 input files, because fibers A and B are
    # separate files
    if len(fnames) == 4:
        demod2(fnames)
    elif len(fnames) == 8:
        demod4(fnames)
    else:
        print("Wrong number of arguments")
        exit(1)
