#first model is 'Tidal buldge height approach' (Cuntz et al.,2000):
#del_g/g = (M_pl/M_star)*(2.*R_star**3.)/(a-R_star)**3.


import numpy as np
from astropy.io import ascii
from astropy.table import Column

TSI = ascii.read('tidal/params.csv',format = 'csv')#system,P_rot[d],P_rot_err[d],R_star[R_sun],R_star_err[R_sun],M_star[M_sun],M_star_err[M_sun],P_orb[d],P_orb_err[d],M_pl[M_jup],M_pl_err[M_jup],a[AU],a_err[AU]

grav_pert = (TSI['M_pl'] / 
            1047.57 / 
            TSI['M_star'] * 2. *
            (TSI['R_star'])**3. / 
            (TSI['a'] * 215.032 - TSI['R_star'])**3.)

TSI.add_column(Column(grav_pert,name = 'grav_pert'))

grav_pert_up_err = ((TSI['M_pl'] + TSI['M_pl_up_err']) / 
                    1047.57 / 
                   ( TSI['M_star'] - TSI['M_star_low_err']) * 2. *
                   (TSI['R_star'] + TSI['R_star_up_err'])**3. / 
                    ((TSI['a'] - TSI['a_err']) * 215.032 - (TSI['R_star'] - TSI['R_star_low_err']))**3. - grav_pert)

grav_pert_low_err = (grav_pert - (TSI['M_pl'] - TSI['M_pl_low_err']) / 
                    1047.57 / 
                   ( TSI['M_star'] + TSI['M_star_up_err']) * 2. *
                   (TSI['R_star'] - TSI['R_star_low_err'])**3. / 
                    ((TSI['a'] + TSI['a_err']) * 215.032 - (TSI['R_star'] + TSI['R_star_up_err']))**3.)

#grav_pert_up_err = (TSI['M_pl'] / 1047.57 /TSI['M_star'] * 2. * TSI['R_star']**3. / 
#                    ( TSI['a'] * 215.032 - TSI['R_star'] )**3. * 
#                    np.sqrt((TSI['M_pl_up_err'] / TSI['M_pl'])**2. + 
#                            (TSI['M_star_up_err'] / TSI['M_star'])**2. + 
#                            (3 * TSI['a_err'] * 215.032 / (TSI['a'] * 215.032 - TSI['R_star']) )**2. + 
#                            (3 * TSI['a'] * 215.032 / (TSI['a'] * 215.032 - TSI['R_star']) * TSI['R_star_up_err'] / TSI['R_star'])**2.)
            #               )
# grav_pert_err = grav_pert_err.filled(np.nan)


# grav_pert_low_err = TSI['M_pl']/1047.57/TSI['M_star']*2.*TSI['R_star']**3./(TSI['a']*215.032-TSI['R_star'])**3.*np.sqrt((TSI['M_pl_low_err']/TSI['M_pl'])**2.+(TSI['M_star_low_err']/TSI['M_star'])**2.+(3*TSI['a_err']*215.032/(TSI['a']*215.032-TSI['R_star']))**2.+(3*TSI['a']*215.032/(TSI['a']*215.032-TSI['R_star'])*TSI['R_star_low_err']/TSI['R_star'])**2.)
# grav_pert_err = grav_pert_err.filled(np.nan)
TSI.add_column(Column(grav_pert_low_err,name = 'grav_pert_low_err'))
TSI.add_column(Column(grav_pert_up_err,name = 'grav_pert_up_err'))



#second model is 'Tidal time skale approach' (Albrecht et al.,2012):
q = TSI['M_pl'] / 1047.57 / TSI['M_star']

q_up = (TSI['M_pl'] + TSI['M_pl_up_err']) / 1047.57 / (TSI['M_star'] - TSI['M_star_low_err'])
q_low = (TSI['M_pl'] - TSI['M_pl_low_err']) / 1047.57 / (TSI['M_star'] + TSI['M_star_up_err'])
#q is the planetary to stellar mass ratio
tau_inv = []
tau_inv_up = []
tau_inv_low = []
for i in range(len(TSI)):
    if TSI['Teff'][i] < 6250.:
        #convective envelope
        tau_inv_el = ( 1 / (10 * 10**9.) * q[i]**2. * (TSI['a'][i] * 215.032 / TSI['R_star'][i] / 40.)**(-6.)) * (5 * 10**9.)
        tau_inv_el_up = (( 1 / 
                         (10 * 10**9.) * 
                         q_up[i]**2. * 
                         ((TSI['a'][i] + TSI['a_err'][i]) *
                         215.032 / 
                         ( TSI['R_star'][i] - TSI['R_star_low_err'][i]) / 40.)**(-6.)) *
                          (5 * 10**9.))
        tau_inv_el_low = (( 1 / 
                         (10 * 10**9.) * 
                         q_low[i]**2. * 
                         ((TSI['a'][i] - TSI['a_err'][i]) *
                         215.032 / 
                         ( TSI['R_star'][i] + TSI['R_star_up_err'][i]) / 40.)**(-6.)) *
                          (5 * 10**9.))
    else:
        #radiative envelope
        raise ValueError
    tau_inv.append(tau_inv_el)
    tau_inv_up.append(tau_inv_el_up)
    tau_inv_low.append(tau_inv_el_low)

tau_inv = np.asarray(tau_inv)


tau = 1./tau_inv

tau_up =  1 / np.asarray(tau_inv_low)
tau_low = 1 / np.asarray(tau_inv_up) 

tau_up_err = np.amax(np.asarray([tau_up, tau_low]), axis=0) - tau
tau_low_err = tau - np.amin(np.asarray([tau_up, tau_low]), axis=0) 


TSI.add_column(Column(tau,name = 'tidal_disip_timescale'))
TSI.add_column(Column(tau_up_err,name = 'tidal_disip_timescale_up_err'))
TSI.add_column(Column(tau_low_err,name = 'tidal_disip_timescale_low_err'))




G = 1.90509*10.**5. #Rsun/Msun*(km/s)**2.
Q_star = 10.**7.
a = TSI['a']*215.032 #in solar radii
a_up = (TSI['a'] + TSI['a_err']) * 215.032
a_low = (TSI['a'] - TSI['a_err']) * 215.032
#thirs model is 'Modelling of tidal interactions and stellar wind' (Penev et al.,2012)
#torque_conv = -0.5*M_pl*TSI['M_star']*np.sqrt(G/(a*(TSI['M_star']+M_pl)))*np.sign((1./TSI['P_rot'])-(1./TSI['P_orb']))*9./2.*np.sqrt(G/(a*TSI['M_star']))*(TSI['R_star']/a)**5.*M_pl/Q_star -- can be rewritten as below:
torque_conv = (-9. / 4. * 
               np.sign((1./TSI['P_rot'])-(1./TSI['P_orb'])) *
               G / 
               Q_star * 
               TSI['R_star']**5./ 
               a**6. * 
               np.sqrt(1. / (1. + q)) * (TSI['M_pl'] / 1047.57)**2.)

torque_conv_up_err = (- 9. / 4. * 
                      np.sign((1./TSI['P_rot'])-(1./TSI['P_orb'])) *
                      G / 
                      Q_star * 
                      (TSI['R_star'] + TSI['R_star_up_err'])**5./ 
                      a_low**6. * 
                      np.sqrt(1. / (1. + q)) *
                      ((TSI['M_pl'] + TSI['M_pl_up_err']) / 1047.57)**2. - torque_conv)


torque_conv_low_err = (torque_conv + 9. / 4. * 
                      np.sign((1./TSI['P_rot'])-(1./TSI['P_orb'])) *
                      G / 
                      Q_star * 
                      (TSI['R_star'] - TSI['R_star_low_err'])**5. / 
                      a_up**6. * 
                      np.sqrt(1. / (1. + q)) *
                      ((TSI['M_pl'] - TSI['M_pl_low_err']) / 1047.57)**2.)

#G is the gravitational constant
#the omegas are the rotational angular frequency of the convectiv envelope and the planetary orbit, respectively
#the Q_star is the stellar tidal quality factor =~ 10**7.

#torque_conv = torque_conv.filled(np.nan)
#torque_conv_err = torque_conv_err.filled(np.nan)

TSI.add_column(Column(torque_conv,name = 'torque_conv'))
TSI.add_column(Column(torque_conv_up_err,name = 'torque_conv_up_err'))
TSI.add_column(Column(torque_conv_low_err,name = 'torque_conv_low_err'))

ascii.write(TSI,'tidal/TIS_all.csv',format='csv',overwrite=True)
