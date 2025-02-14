###############
#
# Keywords for polarization
#
##############


kw_insmode_pol = 'HARPSPOL'             # This is what kw_insmode should contain for polarisation
kw_pol_ret25 = 'HIERARCH ESO INS RET25 POS' # The keyword name for the lambda/4 plate position
kw_pol_ret50 = 'HIERARCH ESO INS RET50 POS' # The keyword name for the lambda/2 plate position
kw_unknown = 'UNKNOWN'

kw_pol_stokes = 'HIERARCH ESO DRS POL STOKES'
kw_pol_stokes_com = 'Extracted Stokes parameter'
kw_pol_stokes_v = 'V/I'
kw_pol_stokes_q = 'Q/I'
kw_pol_stokes_u = 'U/I'
kw_pol_stokes_i = 'I'
kw_pol_stokes_n = 'Null'

kw_pol_demod = 'HIERARCH ESO DRS POL DEMOD'
kw_pol_demod_com = 'Demodulation algorithm'
kw_pol_demod_ratio = 'RATIO'
kw_pol_demod_diff = 'DIFFERENCE'


#HIERARCH ESO DRS POL STOKES  = 'V/I'        / Extracted Stokes parameter
#HIERARCH ESO DRS POL STOKES  = 'Q/I'        / Extracted Stokes parameter
#HIERARCH ESO DRS POL STOKES  = 'U/I'        / Extracted Stokes parameter
#HIERARCH ESO DRS POL STOKES  = 'I'          / Extracted Stokes parameter
#HIERARCH ESO DRS POL DEMOD   = 'DIFFERENCE' / Demodulation algorithm
#HIERARCH ESO DRS POL DEMOD   = 'RATIO'      / Demodulation algorithm
