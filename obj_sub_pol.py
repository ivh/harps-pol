#!/usr/bin/env DRS.DRS_name.
####################################################################
#  obj_sub_pol.py [night_directory] [fitsfilename]
#
#  Polarization. This file is a variant of obj_TWO_harps.py
#  and extracts the spectra from individual raw frames for
#  use in obj_pol.py.
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

import sys,string,time,os
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

s=string
VISU=hadgtVISU
tmp_file="%s"%int(time.time())


fitsfilename=os.path.join(dir_data_raw,arg_night_name,arg_file_names[0])

dirfits=os.path.join(dir_data_reduc,arg_night_name)
insmode=string.strip(hadmrFITS.read_key(fitsfilename,kw_insmode))
WLOG('info',log_opt,'INSTRUMENT MODE: '+insmode)

if (string.find(insmode,ic_insmode)<0):
   WLOG('error',log_opt,'Wrong Image keyword INS MODE: '+insmode+ ' (should be : '+ic_insmode+ ')')
   WLOG('info',log_opt,'Recipe '+process_running+' is terminated')
   exit(1)


if os.path.exists(os.path.join(dir_data_reduc,arg_night_name))==0:
   WLOG('error',log_opt,'Directory :'+dir_data_reduc+arg_night_name+' does not exist')
   exit(1)
 
e2ds_fitsfilename={\
     'A': os.path.join(dir_data_reduc,arg_night_name,string.replace(arg_file_names[0],'.fits','_e2ds_A.fits')),\
     'B': os.path.join(dir_data_reduc,arg_night_name,string.replace(arg_file_names[0],'.fits','_e2ds_B.fits'))}
s1d_fitsfilename={\
     'A': os.path.join(dir_data_reduc,arg_night_name,string.replace(arg_file_names[0],'.fits','_s1d_A.fits')),\
     'B': os.path.join(dir_data_reduc,arg_night_name,string.replace(arg_file_names[0],'.fits','_s1d_B.fits'))}

#
# Read and load CalibDataBase 
#
dic_db=hadmrCDB.get_master(ic_calib_db_master_file,hadmrCDB.filename2tunix(os.path.basename(fitsfilename)),log_opt)
if dic_db.has_key('NONE'): 
   WLOG('error',log_opt,'Missing CalibDB can not proceed')
   exit(0)
hadmrCDB.cp_db_files(dir_calib_db,dirfits,dic_db,log_opt)

#
# Read image keywords
#
execfile(ic_exec_kw_list)

fiber='A'

if os.path.exists(fitsfilename):
    _kw_out_=fitsio.read_keys(fitsfilename,kw,hdu=-1,verbose=0)
elif os.path.exists(e2ds_fitsfilename[fiber]):
    _kw_out_=fitsio.read_keys(e2ds_fitsfilename[fiber],kw,hdu=-1,verbose=0)
else:
     WLOG('error',log_opt,'Files '+fitsfilename+' AND '+e2ds_fitsfilename[fiber]+' missing')
     exit(0)


WLOG('info',log_opt,'Now processing Image TYPE: '+_kw_out_[kw_dprtype][0]+'  with '+process_running+' recipe')
dprtype=_kw_out_[kw_dprtype][0]

execfile(ic_exec_kw_allocate)
execfile(ic_exec_kw_display)


if not (os.path.exists(e2ds_fitsfilename['A']) and os.path.exists(e2ds_fitsfilename['B'])): 
  ccd_gain1=0
  ccd_gain2=0
  ccd_gain3=0
  ccd_gain4=0
  ccd_noise1=0
  ccd_noise2=0  
  ccd_noise3=0  
  ccd_noise4=0  
  WLOG('',log_opt,'Reading Image: '+fitsfilename)
  if hadmrFITS.read_ext(fitsfilename):
     if (locals().has_key('ic_ccdx1')): ccd_gain1=hadmrFITS.read_key(fitsfilename,kw_det_mout1_conad,hdu=kw_det_mout1_HDU)
     if (locals().has_key('ic_ccdx2')): ccd_gain2=hadmrFITS.read_key(fitsfilename,kw_det_mout2_conad,hdu=kw_det_mout2_HDU)
     if (locals().has_key('ic_ccdx3')): ccd_gain3=hadmrFITS.read_key(fitsfilename,kw_det_mout3_conad,hdu=kw_det_mout3_HDU)
     if (locals().has_key('ic_ccdx4')): ccd_gain4=hadmrFITS.read_key(fitsfilename,kw_det_mout4_conad,hdu=kw_det_mout4_HDU)
     if (locals().has_key('ic_ccdx1')): ccd_noise1=hadmrFITS.read_key(fitsfilename,kw_det_mout1_ron,hdu=kw_det_mout1_HDU)
     if (locals().has_key('ic_ccdx2')): ccd_noise2=hadmrFITS.read_key(fitsfilename,kw_det_mout2_ron,hdu=kw_det_mout2_HDU)
     if (locals().has_key('ic_ccdx3')): ccd_noise3=hadmrFITS.read_key(fitsfilename,kw_det_mout3_ron,hdu=kw_det_mout3_HDU)
     if (locals().has_key('ic_ccdx4')): ccd_noise4=hadmrFITS.read_key(fitsfilename,kw_det_mout4_ron,hdu=kw_det_mout4_HDU)
  else :
     if (locals().has_key('ic_ccdx1')): ccd_gain1=hadmrFITS.read_key(fitsfilename,kw_det_out1_conad)
     if (locals().has_key('ic_ccdx2')): ccd_gain2=hadmrFITS.read_key(fitsfilename,kw_det_out2_conad)
     if (locals().has_key('ic_ccdx3')): ccd_gain3=hadmrFITS.read_key(fitsfilename,kw_det_out3_conad)
     if (locals().has_key('ic_ccdx4')): ccd_gain4=hadmrFITS.read_key(fitsfilename,kw_det_out4_conad)
     if (locals().has_key('ic_ccdx1')): ccd_noise1=hadmrFITS.read_key(fitsfilename,kw_det_out1_ron)
     if (locals().has_key('ic_ccdx2')): ccd_noise2=hadmrFITS.read_key(fitsfilename,kw_det_out2_ron)
     if (locals().has_key('ic_ccdx3')): ccd_noise3=hadmrFITS.read_key(fitsfilename,kw_det_out3_ron)
     if (locals().has_key('ic_ccdx4')): ccd_noise4=hadmrFITS.read_key(fitsfilename,kw_det_out4_ron)

 
  data,nx,ny=hadmrFITS.read_data_all(fitsfilename)
  if len(arg_file_names)>1:
     for i in range(len(arg_file_names)-1):
        fitsfilename=os.path.join(dir_data_raw,arg_night_name,arg_file_names[i+1])
        WLOG('',log_opt,'Reading Image '+fitsfilename+ ' and added')
        data2,nx,ny=hadmrFITS.read_data_all(fitsfilename)
        data=add(data,data2)

  WLOG('',log_opt,'Image '+repr(nx)+'x'+repr(ny)+' loaded')

#############################
# Correction of bad columns #
#############################

  WLOG('',log_opt,'Doing Bad columns Correction')

  hadmrBIAS.cosmetic1(data,ic_ccdbadcol1)
  hadmrBIAS.cosmetic2(data,ic_ccdbadcol2)
  hadmrBIAS.cosmetic3(data,ic_ccdbadcol3)

######################
# Correction of BIAS #
######################

  WLOG('',log_opt,'Doing BIAS Measurement')
  meas_ccdsigdet=0
  if (locals().has_key('ic_ccdx1')): 
          hadmrBIAS.cor_real_bias(data,ic_ccdx1,ic_ccdxbias1,ic_smooth_size_bias)
          overscan1=hadmrBIAS.meas_bias(data,ic_ccdxbias1,ic_ccdybias)
          WLOG('',log_opt,'Overscan Zone #1 : mean = %.2f ADU - rms = %.3f ADU' %(overscan1))
          meas_ccdsigdet=meas_ccdsigdet+overscan1[1]

  if (locals().has_key('ic_ccdx2')): 
          hadmrBIAS.cor_real_bias(data,ic_ccdx2,ic_ccdxbias2,ic_smooth_size_bias)
          overscan2=hadmrBIAS.meas_bias(data,ic_ccdxbias2,ic_ccdybias)
          WLOG('',log_opt,'Overscan Zone #2 : mean = %.2f ADU - rms = %.3f ADU' %(overscan2))
          meas_ccdsigdet=meas_ccdsigdet+overscan2[1]

  if (locals().has_key('ic_ccdx3')): 
          hadmrBIAS.cor_real_bias(data,ic_ccdx3,ic_ccdxbias3,ic_smooth_size_bias)
          overscan3=hadmrBIAS.meas_bias(data,ic_ccdxbias3,ic_ccdybias)
          WLOG('',log_opt,'Overscan Zone #3 : mean = %.2f ADU - rms = %.3f ADU' %(overscan3))
          meas_ccdsigdet=meas_ccdsigdet+overscan3[1]
		
  if (locals().has_key('ic_ccdx4')): 
          hadmrBIAS.cor_real_bias(data,ic_ccdx4,ic_ccdxbias4,ic_smooth_size_bias)
          overscan4=hadmrBIAS.meas_bias(data,ic_ccdxbias4,ic_ccdybias)
          WLOG('',log_opt,'Overscan Zone #4 : mean = %.2f ADU - rms = %.3f ADU' %(overscan4))
          meas_ccdsigdet=meas_ccdsigdet+overscan4[1]
  
  meas_ccdsigdet=meas_ccdsigdet/(locals().has_key('ic_ccdx1')+locals().has_key('ic_ccdx2')+\
                               locals().has_key('ic_ccdx3')+locals().has_key('ic_ccdx4'))/sqrt(len(arg_file_names))
	
  WLOG('',log_opt,'Doing BIAS Correction')
  ccdgain=max(ccd_gain1,ccd_gain2,ccd_gain3,ccd_gain4)
  try:
    ccdsigdet=max(ccd_noise1,ccd_noise2,ccd_noise3,ccd_noise4)
    diff_ccdsigdet=ccdsigdet-meas_ccdsigdet*ccdgain
    if abs(diff_ccdsigdet)>ic_diff_ccdsigdet:
       WLOG('',log_opt,'Differences between expected and measured CCD RON is unusual %.3f, using measured value'%diff_ccdsigdet)
       ccdsigdet=meas_ccdsigdet*ccdgain

  except TypeError:
        ccdsigdet=meas_ccdsigdet*ccdgain

  noise=ccdsigdet*sqrt(ic_extmeanzone)
  WLOG('',log_opt,'Using CCD RON: %.3f[e-] and spectral pixel RON %.3f[-e]'%(ccdsigdet,noise))
  kw_CCD_SIGDET[1]=ccdsigdet
  kw_CCD_CONAD[1]=ccdgain

######################
# Correction of DARK #
######################

  WLOG('',log_opt,'Doing DARK Substraction')
  meandark=exp_time*(ic_ccddark/ccdgain)/3600.
  WLOG('',log_opt,'Estimated  DARK level = %.2f ADU' %meandark+' used for correction')
  if meandark > 1:
    data=data-int(meandark+0.5)

#######################
# Resize of the image #
#######################

  data2,nx2,ny2=hadmrBIAS.redim(data,ic_ccdx,ic_ccdy)
  WLOG('',log_opt,'Image format changed to '+repr(nx2)+'x'+repr(ny2))
  if locals().has_key('ic_exec_image_geometry'): 
    WLOG('',log_opt,'Doing geometry correction')
    execfile(ic_exec_image_geometry)




################################
#  Earth Velocity calculation  #
################################


  WLOG('',log_opt,'Computing Earth RV correction')
  berv,bjd,berv_max=newbervmain(target_alpha,target_delta,target_equinox,obs_year,obs_month,obs_day,obs_hour,\
                                ic_longit_obs,ic_latit_obs,ic_altit_obs,target_rapm,target_decpm)
  WLOG('info',log_opt,'Barycentric Earth RV correction: %9.5f km/s' %(berv))
 

#####################
# Orders Extraction #
#####################
  

wave={}
e2dsff={}
blaze={}
param_ll={}
for fiber in ['A','B']:
  if os.path.exists(e2ds_fitsfilename[fiber]):
    WLOG('info',log_opt,'Reading existing E2DS file '+os.path.split(e2ds_fitsfilename[fiber])[1])
    e2dsff[fiber]=fitsio.read_data(e2ds_fitsfilename[fiber],hdu=-1,verbose=0)[:,:]
    wave[fiber],param_ll[fiber],param_x_=hadmrFITS.get_e2ds_ll(e2ds_fitsfilename[fiber],_kw_TH_CAL_)
    ccdsigdet=hadmrFITS.read_key(e2ds_fitsfilename[fiber],kw_CCD_SIGDET[0])
    ccdgain=hadmrFITS.read_key(e2ds_fitsfilename[fiber],kw_CCD_CONAD[0])
    noise=ccdsigdet*sqrt(ic_extmeanzone)
  else:  
    
#
# get localization param
#
    if dic_db.has_key('LOC_'+fiber)==0:
      WLOG('error',log_opt,'No order geometry is defined in the calibDB for fiber: '+fiber)
      exit(3)
    WLOG('',log_opt,'Reading localization parameters of Fiber '+fiber)
    loco_file=os.path.join(dirfits,dic_db['LOC_'+fiber][1])
    nbo=int(hadmrFITS.read_key(loco_file,kw_LOCO_NBO[0]))
    nbcoeff_ctr=int(hadmrFITS.read_key(loco_file,kw_LOCO_DEG_C[0]))+1
    nbcoeff_width=int(hadmrFITS.read_key(loco_file,kw_LOCO_DEG_W[0]))+1
    delta_gau=hadmrFITS.read_key(loco_file,kw_LOCO_DELTA[0])
    ac=hadmrEXTOR.readkeyloco(loco_file,kw_LOCO_CTR_COEFF[0],nbo,nbcoeff_ctr)
    as=hadmrEXTOR.readkeyloco(loco_file,kw_LOCO_FWHM_COEFF[0],nbo,nbcoeff_width)

    
    e2ds=zeros((nbo,ny2),'d')
    nbcos=zeros(nbo,'i')
    flat=zeros((nbo,ny2),'d')
    nbcos=zeros(nbo,'i')
    S_N=zeros(nbo,'d')
  
#
# extraction
#
    WLOG('',log_opt,'On fiber '+fiber+' doing  extraction of '+repr(nbo)+' orders')   
    saturation=0
    for no in range(0,nbo):
      e2ds[no],nbcos[no]=extract2(data2.ravel(),ac[no],as[no],ny2,
           ic_extopt,ic_extnbsig,ccdgain,ccdsigdet,delta_gau,ic_extseuilcosmic)
      dd,cc=extract2(data2.ravel(),ac[no],as[no],ny2,
           2,ic_extnbsig,ccdgain,ccdsigdet,delta_gau,ic_extseuilcosmic)
# cut outstanding pixel hit
      hbin=arange(min(e2ds[no]),max(e2ds[no]),(max(e2ds[no])-min(e2ds[no]))/4)
      H=histogram(e2ds[no],hbin)
      detect_lim=max(e2ds[no])     
      if ((H[2])==0)&(H[3]>0):      
          nbcos[no]=nbcos[no]+H[3]
	  detect_lim=(hbin[2])
      e2ds[no]=where(greater(e2ds[no],detect_lim),e2ds[no]*0+(hbin[1]+hbin[0])/2,e2ds[no])
      if ic_debug: 
        print H
        hadgtVISU.megaplot(range(ny2),dd,ymax=max(e2ds[no]))
        hadgtVISU.megaplot(range(ny2),e2ds[no],overplot=1,title='order %s'%(no),ymax=max(e2ds[no]))
        time.sleep(ic_disptimeout)
	  
      det_noise=sqrt(ic_extmeanzone)*ccdsigdet
      flux=sum(e2ds[no,ny2/2-ic_extfblaz:ny2/2+ic_extfblaz])/(2.*ic_extfblaz)
      S_N[no]=flux/sqrt(flux+det_noise**2)
      WLOG('',log_opt,'On fiber '+fiber+'  order %2s' %(repr(no))+\
                         ': S/N= %.1f  %d cosmics found' %(S_N[no],nbcos[no]))
# QC Saturation
      if ((flux/ccdgain)/ic_meanwidth)>qc_max_signal: 
	     saturation=1
          
      if ic_debug: hadgtVISU.plotxy(range(ny2),e2ds[no])

    if saturation: WLOG('warning',log_opt,'Saturation level reached')

    WLOG('info',log_opt,'On fiber '+fiber+\
         ': S/N[%snm]= %.1f S/N[%snm]= %.1f S/N[%snm]= %.1f' %\
	 (ic_sn_display_ll[fiber][0],S_N[ic_sn_display_o[fiber][0]],ic_sn_display_ll[fiber][1],S_N[ic_sn_display_o[fiber][1]],ic_sn_display_ll[fiber][2],S_N[ic_sn_display_o[fiber][2]])+\
	 '   (%d cosmic removed)'%(sum(nbcos)))
  
#
# Do Flat Field correction  
#
    if dic_db.has_key('FLAT_'+fiber)==0:
       WLOG('error',log_opt,'No Flat-Field is defined in the calibDB for fiber: '+fiber)
       exit(3)

    flat_file=os.path.join(dirfits,dic_db['FLAT_'+fiber][1])
    flat,fx,fy=hadmrFITS.read_data(flat_file)
    WLOG('',log_opt,'Doing Flat-Field correction for fiber: '+fiber+' using '+dic_db['FLAT_'+fiber][1])

    e2dsff[fiber]=e2ds/flat



###############################
# Insert wavelength calibration  
###############################

    if dic_db.has_key('TH_'+fiber)==0:
      WLOG('error',log_opt,'No wavelength Calibration for fiber '+fiber+'  FATAL!')
      exit(3)

    th_file=os.path.join(dirfits,dic_db['TH_'+fiber][1])
    wave[fiber],param_ll[fiber],param_x=hadmrFITS.get_e2ds_ll(th_file,_kw_TH_CAL_) 
    WLOG('',log_opt,'Wavelength calibration is set using '+os.path.split(th_file)[1])


#############################
# Save E2DSFF  on fits file #
#############################
		
    WLOG('',log_opt,'Saving E2DS spectrum of Fiber '+fiber+' in '+os.path.split(e2ds_fitsfilename[fiber])[1])

    hadmrFITS.write_data(e2ds_fitsfilename[fiber],e2dsff[fiber])
# Copy keywords of the raw frame 
    hadmrFITS.copy_key(fitsfilename,e2ds_fitsfilename[fiber])
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_version)

    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_CCD_SIGDET)
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_CCD_CONAD)

# Copy keywords of the localization frame 
    kw_LOCO_FILE[1]=dic_db['LOC_'+fiber][1]
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_LOCO_FILE)
    hadmrFITS.copy_keys_root(loco_file,e2ds_fitsfilename[fiber],kw_root_drs_loc)
# Copy keywords of the flat field frame 
    kw_FLAT_FILE[1]=dic_db['FLAT_'+fiber][1]
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_FLAT_FILE)
    hadmrFITS.copy_keys_root(flat_file,e2ds_fitsfilename[fiber],kw_root_drs_flat)
# Copy keywords of the blaze frame 
    kw_BLAZE_FILE[1]=dic_db['BLAZE_'+fiber][1]
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_BLAZE_FILE)
# Add keywords of the extraction process 
    kw_EXTRA_OPT[1]=ic_extopt
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_EXTRA_OPT)
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_EXTRA_SIG)
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_EXTRA_COSM)
# Add keywords of the S/N 
    for no in range(0,nbo):
      SNname=kw_EXTRA_SN[0]+repr(no)
      SNval=int(S_N[no]*10)/10.
      SNcom=kw_EXTRA_SN[2]+repr(no)
      hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],[SNname,SNval,SNcom])
    for no in range(0,nbo):
      Nbcosname=kw_EXTRA_NBCOS[0]+repr(no)
      Nbcosval=nbcos[no]
      Nbcoscom=kw_EXTRA_NBCOS[2]+repr(no)
      hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],[Nbcosname,Nbcosval,Nbcoscom])
# Add keywords of the RV earth correction process 
    kw_BERV[1]=berv
    kw_BJD[1]=bjd
    kw_BERVMX[1]=berv_max    
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_BERV)
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_BJD)
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_BERVMX)
# Copy keywords of the wavelength calibration
    kw_TH_FILE[1]=dic_db['TH_'+fiber][1]
    kw_LAMP_USED[1]=hadmrFITS.read_key(th_file,kw_lamp_id)
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_TH_FILE)
    hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_LAMP_USED)
    hadmrFITS.copy_keys_root(th_file,e2ds_fitsfilename[fiber],kw_root_drs_th)


##### NO INSTRUMENTAL DRIFT COMPUTATION  

kw_DRIFT_NOISE[1] = ic_typical_drift

for fiber in ['A','B']:
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_REF)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_VR)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_NBCOS)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_RFLUX)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_NBORDKILL)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_NOISE)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_REF_CCF)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_REF_RV)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_VR_CCF)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_VR_USED)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_ALGO)
  hadmrFITS.update_key(e2ds_fitsfilename[fiber],kw_DRIFT_QC)




for fiber in ['A','B']:

#
# Cosmic filter correction
#

  newcosmic=0
  WLOG('',log_opt,'Doing Cosmic filtering on fiber: '+fiber)

  for no in range(0,len(e2dsff[fiber])):
     e2dsff[fiber][no],nbcos2=hadmrEXTOR.cosmic_filter(e2dsff[fiber][no])
     newcosmic=newcosmic+nbcos2
    
  WLOG('info',log_opt,'Number of Cosmic corrected:  %i'%(newcosmic)+' on fiber: '+fiber)



#
# divide by blaze and save
#

 
for fiber in ['A','B']:
  if dic_db.has_key('BLAZE_'+fiber)==0:
       WLOG('error',log_opt,'No Blaze is defined in the calibDB for fiber: '+fiber)
       exit(3)
  blaze_file=os.path.join(dirfits,dic_db['BLAZE_'+fiber][1])
  blaze[fiber],fx,fy=hadmrFITS.read_data(blaze_file)
  e2dsffb=e2dsff[fiber]/blaze[fiber]

  e2dsffb_fitsfilename=e2ds_fitsfilename[fiber].replace('e2ds','e2dsffb')
  WLOG('',log_opt,'Saving E2DSFFB spectrum of Fiber '+fiber+' in '+os.path.split(e2dsffb_fitsfilename)[1])
  
  hadmrFITS.write_data(e2dsffb_fitsfilename,e2dsffb)
  hadmrFITS.copy_key(e2ds_fitsfilename[fiber],e2dsffb_fitsfilename)

# If any new keywords need to be set, do it here.
#  hadmrFITS.write_newkey(e2ds_fitsfilename[fiber],kw_version)
