#!/usr/bin/env DRS.DRS_name.
####################################################################
#  obj_pol.py [night_directory] [fitsfilenames]
#
#  Polarimetry with HARPS
#
###################################################################


import hadmrFITS
import hadgtVISU
import hadmrBIAS
import hadmrDARK
import hadmrEXTOR
import hadmrTHORCA
import hadmrCDB
import hadmrFLAT
import hadmrRV
import hadrgdCONFIG

import time,os
from sys import exit
import string as s
import shutil
import fitsio
import fit
from numpy.oldnumeric import *
from extract2 import *
from newbervmain import *
from e2dstos1d import *
from hadmrMATH import *

execfile(os.getenv('PYTHONSTARTUP'))
execfile(os.getenv('PYTHONSTARTUP2'))

DEBUG=False

# set variables with additional keyword names
execfile(os.path.join(dir_drs_config,'kw_pol.py'))


Vangles=[45.0, 135.0, 225.0, 315.0]
Qangles=[0.0, 45.0, 90.0, 135.0, 180, 225.0, 270.0, 315.0]
Uangles=map(add,Qangles,[22.5]*8)

Vp1=Vangles[:2]
Vp2=Vangles[2:]
Vpairs=[Vp1,Vp2]

Qpairs=[]
Upairs=[]
for i in arange(4):
    Qpairs.append(Qangles[i*2:(i+1)*2])
    Upairs.append(Uangles[i*2:(i+1)*2])

WLOG('info',log_opt,'Starting to classify the given files: %s'%str(arg_file_names))

insmodes=[]
ret25=[]
ret50=[]

firstfile=arg_file_names[0]

for filename in arg_file_names:
    fitsfilename=os.path.join(dir_data_raw,arg_night_name,filename)
    if DEBUG: print "fitsfilename: %s"%fitsfilename
    
    insmodes.append(s.strip(hadmrFITS.read_key(fitsfilename,kw_insmode)))
    r25=hadmrFITS.read_key(fitsfilename,kw_pol_ret25)
    if isinstance(r25,float): ret25.append(r25)

    r50=hadmrFITS.read_key(fitsfilename,kw_pol_ret50)
    if isinstance(r50,float): 
        if r50!=360.0:
            ret50.append(r50)
        else:
            ret50.append(0.0)
    

# check if all files are polarisation files
if insmodes.count(kw_insmode_pol) != len(arg_file_names):
    WLOG('error',log_opt,'Not all given files are polarization measurements.')
    exit(1)

if len(ret25) and (not len(ret50)):
    CPOL=True
    angles=ret25
elif len(ret50) and (not len(ret25)):
    CPOL=False
    angles=ret50
else:
    WLOG('error',log_opt,'Recieved a mixture of RET25 and RET50, or empty angle list')
    exit(1)

# double pairs or single
if len(angles)==4: DOUB=True
elif len(angles)==2: DOUB=False
else:
    WLOG('error',log_opt,'We should either have two or four files given.')
    exit(1)


## SORT THE LISTS 
tmp=zip(angles,arg_file_names)
tmp.sort()
angles, arg_file_names = map(lambda t: list(t), zip(*tmp))
if DEBUG: print angles, arg_file_names

def issubset(a,b):
    """ tests if all elemets of a are in b"""
    for i in a:
        if i not in b:
            return False
    return True

# V, Q or U
if issubset(angles,Vangles) and CPOL:
    STOKES='V'
    pairs=Vpairs
elif issubset(angles,Uangles):
    STOKES='U'
    pairs=Upairs
elif issubset(angles,Qangles):
    STOKES='Q'
    pairs=Qpairs
else:
    WLOG('error',log_opt,'Retarder angles do not belong to either V, Q or U: %s'%str(angles)) 
    exit(1)

# check that we really have angle pairs
if angles[:2] not in pairs:
    WLOG('error',log_opt,'Angles are not in pairs: %s'%str(angles)) 
    exit(1)  
if DOUB:
    if angles[2:] not in pairs:
        WLOG('error',log_opt,'Angles are not in pairs: %s'%str(angles)) 
        exit(1)  



WLOG('info',log_opt,'Measuring %s from angles: %s'%(STOKES,str(angles)))

# So far so good, if we have not exited so far, everything looks nice
# and we can start doing the real work, starting with the extraction
# on the raw files.
#

WLOG('info',log_opt,'Starting the loop over the files and running the extraction obj_sub_pol.py on each of them.')
data=[]
for i,filename in enumerate(arg_file_names):
    e2dsffb_fitsfilename={\
        'A': os.path.join(dir_data_reduc,arg_night_name,string.replace(filename,'.fits','_e2dsffb_A.fits')),\
            'B': os.path.join(dir_data_reduc,arg_night_name,string.replace(filename,'.fits','_e2dsffb_B.fits'))}
    
    if not (os.path.exists(e2dsffb_fitsfilename['A']) and os.path.exists(e2dsffb_fitsfilename['B']) ):
        cmd='obj_sub_pol.py %s %s'%(arg_night_name,filename)
        WLOG('info',log_opt,'Running sub-script: %s'%cmd)
        err=os.system('obj_sub_pol.py %s %s'%(arg_night_name,filename))
        if DEBUG: print 'Error code: %d'%err
        if err !=0:
            WLOG('error',log_opt,'Running obj_sub_pol.py on %s failed!'%filename)
            exit(1)

    for fiber in ['A','B']:
        e2dsffb,fx,fy=hadmrFITS.read_data(e2dsffb_fitsfilename[fiber])
        if fiber=='A':
            # exclude order 46 from A, because it is not in B
            e2dsffb=concatenate((e2dsffb[:45,:],e2dsffb[46:,:]))
        data.append(e2dsffb.astype('double'))
        WLOG('info',log_opt,'re-read data from %s.'%e2dsffb_fitsfilename[fiber])

###
### Finished extraction
###
WLOG('info',log_opt,'Finished running obj_sub_pol.py on all files.')

def ratio_single(data):
    """ IMPORTANT: This expects the list "data" to be sorted by angle, then fiber"""
    if len(data) != 4:
        WLOG('error',log_opt,'Function ratio_single expects 4 measurements')
        exit(1)
    
    ang1a,ang1b,ang2a,ang2b=tuple(data)
    
    R = ang1b*ang2a
    R /= ang1a*ang2b
    X = sqrt(R) - 1
    X /= sqrt(R) + 1
    I = 0.25 * (ang1a+ang1b+ang2a+ang2b)
    
    return X,I

def R2X(R):
    """sub-function for the double ratio calculation """
    return (sqrt(sqrt(R)) - 1) / (sqrt(sqrt(R)) + 1)
    
def ratio_double(data):
    """ IMPORTANT: This expects the list "data" to be sorted by angle, then fiber"""
    if len(data) != 8:
        WLOG('error',log_opt,'Function diff_double expects 8 measurements')
        exit(1)

    ang1a,ang1b,ang2a,ang2b,ang3a,ang3b,ang4a,ang4b=tuple(data)
    
    R = ang1b*ang2a*ang3b*ang4a
    R /= ang1a*ang2b*ang3a*ang4b
    RN = ang1b*ang2a*ang3a*ang4b
    RN /= ang1a*ang2b*ang3b*ang4a

    X = R2X(R)
    N = R2X(RN)
    I = 0.125 * (ang1a+ang1b+ang2a+ang2b+ang3a+ang3b+ang4a+ang4b)

    return X,I,N

def diff_single(data):
    """ IMPORTANT: This expects the list "data" to be sorted by angle, then fiber"""
    if len(data) != 4:
        WLOG('error',log_opt,'Function diff_single expects 4 measurements')
        exit(1)
   
    ang1a,ang1b,ang2a,ang2b=tuple(data)
    
    X = (ang1b-ang1a) / (ang1b+ang1a)
    X -= (ang2b-ang2a) / (ang2b+ang2a)
    X *= 0.5

    I = 0.25 * (ang1a+ang1b+ang2a+ang2b)
    
    return X,I

def diff_double(data):
    """ IMPORTANT: This expects the list "data" to be sorted by angle, then fiber"""
    if len(data) != 8:
        WLOG('error',log_opt,'Function diff_double expects 8 measurements')
        exit(1)

    ang1a,ang1b,ang2a,ang2b,ang3a,ang3b,ang4a,ang4b=tuple(data)
    
    X1 = (ang1b-ang1a) / (ang1b+ang1a)
    X2 = (ang2b-ang2a) / (ang2b+ang2a)
    X3 = (ang3b-ang3a) / (ang3b+ang3a)
    X4 = (ang4b-ang4a) / (ang4b+ang4a)
    
    X = X1 - X2 + X3 - X4
    X *= 0.25

    N = X1 - X2 - X3 + X4
    N *= 0.25
    
    I = 0.125 * (ang1a+ang1b+ang2a+ang2b+ang3a+ang3b+ang4a+ang4b)

    return X,I,N

# Set all filenames
filename_base=os.path.join(dir_data_reduc,arg_night_name,firstfile).replace('.fits','')
Q_filename= filename_base+'_q.fits'
V_filename= filename_base+'_v.fits'
U_filename= filename_base+'_u.fits'
N_filename= filename_base+'_Null.fits'
I_filename= filename_base+'_I.fits'

if STOKES=='Q':
    X_filename=Q_filename
elif STOKES=='U':
    X_filename=U_filename
elif STOKES=='V':    
    X_filename=V_filename
else:
    WLOG('error',log_opt,'This should never have happened.')
    exit(1)

header_base=filename_base+'_e2dsffb_A.fits'

def fix_headers(filename,what):
    hadmrFITS.copy_key(header_base,filename)
    
    hadmrFITS.delete_key(filename,kw_pol_ret25) 
    hadmrFITS.delete_key(filename,kw_pol_ret50)

    # which emodulation?
    if 'diff' in filename:
        demod=kw_pol_demod_diff
    else:
        demod=kw_pol_demod_ratio
    hadmrFITS.write_newkey(filename,[kw_pol_demod, demod, kw_pol_demod_com])

    stokes=[kw_pol_stokes, kw_pol_stokes_v, kw_pol_stokes_com]
    if what=='V':
        hadmrFITS.write_newkey(filename,stokes)
    elif what=='Q':
        stokes[1]=kw_pol_stokes_q
        hadmrFITS.write_newkey(filename,stokes)
    elif what=='U':
        stokes[1]=kw_pol_stokes_u
        hadmrFITS.write_newkey(filename,stokes)
    elif what=='N':
        stokes[1]=kw_pol_stokes_n
        hadmrFITS.write_newkey(filename,stokes)
    elif what=='I':
        stokes[1]=kw_pol_stokes_i
        hadmrFITS.write_newkey(filename,stokes)
       


## RUN THE DEMODULATION
WLOG('info',log_opt,'Starting the demodulation calculation.')
  
if DOUB: # We have double set of files

    X,I,N=ratio_double(data)
    
    hadmrFITS.write_data(X_filename,X)
    fix_headers(X_filename,STOKES)
    hadmrFITS.write_data(I_filename,I)
    fix_headers(I_filename,'I')
    hadmrFITS.write_data(N_filename,N)
    fix_headers(N_filename,'N')
 
    X,I,N=diff_double(data)
    
    hadmrFITS.write_data(X_filename.replace('.fits','_diff.fits'),X)
    fix_headers(X_filename.replace('.fits','_diff.fits'),STOKES)
    hadmrFITS.write_data(I_filename.replace('.fits','_diff.fits'),I)
    fix_headers(I_filename.replace('.fits','_diff.fits'),'I')
    hadmrFITS.write_data(N_filename.replace('.fits','_diff.fits'),N)
    fix_headers(N_filename.replace('.fits','_diff.fits'),'N')

else: # We have a single set of files
    
    X,I=ratio_single(data)
    
    hadmrFITS.write_data(X_filename,X)
    fix_headers(X_filename,STOKES)
    hadmrFITS.write_data(I_filename,I)
    fix_headers(I_filename,'I')

    X,I=diff_single(data)
    
    hadmrFITS.write_data(X_filename.replace('.fits','_diff.fits'),X)
    fix_headers(X_filename.replace('.fits','_diff.fits'),STOKES)
    hadmrFITS.write_data(I_filename.replace('.fits','_diff.fits'),I)
    fix_headers(I_filename.replace('.fits','_diff.fits'),'I')

WLOG('info',log_opt,'Finished.')
