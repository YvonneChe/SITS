# SITS
# SITS amber input demo
step 0: output the energy of enhanced part (sluene), 
       non-enhanced part (solene) and their interaction
        

 &cntrl
  imin = 0, irest = 1, ntx = 5,
  ntb = 2, pres0 = 1.0, ntp = 1,
  taup = 2.0, cut = 10, ntr = 0,
  ntc = 2, ntf = 2, ig = -1,
  tempi = 300.0, temp0 = 300.0,
  ntt = 3, gamma_ln = 5.0,
  nstlim = 500000, dt = 0.002,
  ntpr = 1000, ntwx = 1000, ntwr = 1000,
  ioutfm = 1, iwrap=1,
  isub =1
 /
  &subunit
  bgsluatm=1, nsluatm=179, bgsolatm=180
/


step 1: Start SITS enhanced sampling, get converged set of {nk}
        which are bias factors for different temperatures 

 &cntrl
  imin = 0, irest = 1, ntx = 5,
  ntb = 2, pres0 = 1.0, ntp = 1,
  taup = 2.0, cut = 10, ntr = 0,
  ntc = 2, ntf = 2, ig = -1,
  tempi = 300.0, temp0 = 300.0,
  ntt = 3, gamma_ln = 5.0,
  nstlim = 500000, dt = 0.002,
  ntpr = 1000, ntwx = 1000, ntwr = 1000,
  ioutfm = 1, iwrap=1,
  isub =1, doits=1
/
  &subunit
  bgsluatm=1, nsluatm=179, bgsolatm=180
/
  ntemp=200, update_step=100,
  fb_rest=0,fb_const=0,
  templ=260.0, temph=500.0,
  peshift=300.0, fb_bias=0.00
/

step 2: with fixed {nk} run SITS

 &cntrl
  imin = 0, irest = 1, ntx = 5,
  ntb = 2, pres0 = 1.0, ntp = 1,
  taup = 2.0, cut = 10, ntr = 0,
  ntc = 2, ntf = 2, ig = -1,
  tempi = 300.0, temp0 = 300.0,
  ntt = 3, gamma_ln = 5.0,
  nstlim = 500000, dt = 0.002,
  ntpr = 1000, ntwx = 1000, ntwr = 1000,
  ioutfm = 1, iwrap=1,
  isub =1, doits=1
/
  &subunit
  bgsluatm=1, nsluatm=179, bgsolatm=180
/
  ntemp=200, update_step=100,
  fb_rest=1,fb_const=1,
  templ=260.0, temph=500.0,
  peshift=300.0, fb_bias=0.00
/
