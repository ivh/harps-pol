This is an extension to the reduction pipeline for the HARPS 
spectrograph, carrying out the demodulation of the spectra from
fibers A and B in multiple frames (either 2 or 4).

The input to demod.py is a list of filenames. These need to be the S2D
files that are made by the ESPRESSO pipeline, applied to the HARPSPOL raw
framed. The fibers A and B are saved in separate files, therefore the
number of input files to the script is 4 or 8.

Outputs are saved in the same S2D format, with filenames that share the prefix
and timestamp with the first of the files given. The suffixes are
* _S2D_POL_I.fits for the intensity spectrum
* _S2D_POL_STOKES.fits for the Stokes parameters
* _S2D_POL_NULL.fits for the null spectrum (only when 8 files provided)

The files in old/ are the plugin to the old, python-based, HARPS
reduction pipeline.