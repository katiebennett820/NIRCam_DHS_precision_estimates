import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import matplotlib as mpl
import spectres
import sys
from scipy import interpolate
from astropy.io import ascii, fits

'''
Calculates precision for NIRCam DHS observations with using different subarrays/readout patterns/groups per integration. 
Useful for comparing different observing strategies. 

Takes input S/N from ETC and calculates (unbinned) photometric precision (PP):

PP = sqrt(2) * (1/SNR) / sqrt(n_ints_in)

where "n_ints_in" is the number of in-transit integrations. The sqrt(2) factor comes from the fact that we must account for the in and
out of transit data. This assumes a 1:1 ratio and is a conservative estimate even if you are planning for a longer out-of-transit baseline!

Another approach is to take: 

PP = sqrt(1/n_ints_in + 1/n_ints_out) * (1/SNR)

but this assumes your observation will execute exactly as you planned.. 

STEPS:

1. Run ETC for desired subarray, readout pattern, and number of groups.
2. Download .tar files for short-wavelength (SW) and long-wavelength (LW) data. 
3. Run APT file to check number of integrations total (n_ints = length of observation / exposure time per integration). Estimate length of 
    observation via 2 x T_dur + 1 hr. Exposure time per integration is calculated from the exposure time given for 1 integration for a 
    given readout pattern and groups/int. (Can also estimate number of integrations yourself from JDox, but I think this is easier.)
4. Check data excess in APT. 5-15 GB is the "yellow zone". 
5. Add inputs below and run .py file!
6. Outputs: plot of precision, median error bar size in SW and LW, and ascii file of error bars.  

'''

#----------------INPUTS-----------------------------

planet = 'LTT1445Ab'
mode = 'NIRCam'
substripe = 'SUB41S1' #see https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-observing-modes/nircam-time-series-observations/nircam-short-wavelength-grism-time-series#gsc.tab=0
dhs = 'DHS3' #which readout mode?

obs_mode_sw = 'F150W2' #what observing mode are you using for the short-wavelength component? 
obs_mode_lw = 'F444W' #what observing mode are you using for the short-wavelength component? 

#Paths to untarred data from ETC (SW and LW) and where to save error outputs and figures
data_path_sw = '/Users/katiebennett/Documents/Research/Projects/LTT1445Ab_JWST_NIRCam/Observation_Planning/error_estimates/ETC_Calculations/F150W2_F444W/SUB41_DHS3_9grps/SW/'
data_path_lw = '/Users/katiebennett/Documents/Research/Projects/LTT1445Ab_JWST_NIRCam/Observation_Planning/error_estimates/ETC_Calculations/F150W2_F444W/SUB41_DHS3_9grps/LW/'
save_path = '/Users/katiebennett/Documents/Research/Projects/LTT1445Ab_JWST_NIRCam/Observation_Planning/error_estimates/error_estimates/F150W2_F444W/SUB41_DHS3_9grps/'
fig_path = '/Users/katiebennett/Documents/Research/Projects/LTT1445Ab_JWST_NIRCam/Observation_Planning/error_estimates/figures/'

length_obs = 3.73456 # length of observation in hours. 
length_transit = 1.36728 #transit duration (hours)
ngrps = 9 #Calculated in ETC
time_first_last = 5.50 #total science time per int ("Time Between First Reset and Last Measurement, per integration" in ETC)
nint = 2341 #Total number of integrations (estimate from APT)

#----------------DEFINITIONS-----------------------------

def snr_to_ppm(flux, noise, nints):
    '''
    Convert SNR to precision
    '''
    snr  = flux / noise 
    precision = 1 / snr # get relative precision in single integration
    precision_ints = np.sqrt(2) * precision / np.sqrt(nints) #how many integrations do you have? 
    noise_ppm = precision_ints * 1e6
    return noise_ppm

#----------------PULL IN ETC DATA-----------------------------

#Pull flux (in e-/s) and SNR calculated in ETC (per integration)
flx_data_sw = fits.open(data_path_sw+'lineplot/lineplot_extracted_flux.fits')
snr_data_sw = fits.open(data_path_sw+'lineplot/lineplot_sn.fits')
noise_data_sw = fits.open(data_path_sw+'lineplot/lineplot_extracted_noise.fits')

flx_data_lw = fits.open(data_path_lw+'lineplot/lineplot_extracted_flux.fits')
snr_data_lw = fits.open(data_path_lw+'lineplot/lineplot_sn.fits')
noise_data_lw = fits.open(data_path_lw+'lineplot/lineplot_extracted_noise.fits')

#Label wavelength, flux, and SNR variables
wvl_sw = np.array([flx_data_sw[1].data[i][0] for i in range(len(flx_data_sw[1].data))])
flx_sw = np.array([flx_data_sw[1].data[i][1] for i in range(len(flx_data_sw[1].data))])
sn_sw = np.array([snr_data_sw[1].data[i][1] for i in range(len(snr_data_sw[1].data))])
noise_sw = np.array([noise_data_sw[1].data[i][1] for i in range(len(noise_data_sw[1].data))])

wvl_lw = np.array([flx_data_lw[1].data[i][0] for i in range(len(flx_data_lw[1].data))])
flx_lw = np.array([flx_data_lw[1].data[i][1] for i in range(len(flx_data_lw[1].data))])
sn_lw = np.array([snr_data_lw[1].data[i][1] for i in range(len(snr_data_lw[1].data))])
noise_lw = np.array([noise_data_lw[1].data[i][1] for i in range(len(noise_data_lw[1].data))])


#-----PLOT TO MAKE SURE S/N MATCHES ETC PLOT----------
'''
#UNCOMMENT TO PLOT HERE
plt.title('S/N comparison')
plt.plot(wvl_sw, flx_sw/noise_sw, label = 'SNR from given flux and noise')
plt.xlabel('Wavelength (microns)')
plt.ylabel('S/N')
plt.legend()
plt.show()

plt.title('S/N comparison')
plt.plot(wvl_lw, flx_lw/noise_lw, label = 'SNR from given flux and noise')
plt.xlabel('Wavelength (microns)')
plt.ylabel('S/N')
plt.legend()
plt.show()
'''

#-------------CALCULATIONS----------------------------
print('------------------------------')
print('Estimating precision for:', str(planet)+', '+str(mode)+', '+str(obs_mode_sw)+'+'+str(obs_mode_lw)+', '+str(substripe)+', '+ \
    str(dhs)+', '+str(ngrps), 'grps')
print('------------------------------')


#How many integrations are in transit? 
nint_in = int(nint * (length_transit / length_obs))
print('Number of in-transit integrations:', nint_in)
print('------------------------------')

#Calculate spectrophotometric precision from SNR -- THIS IS MORE CONSERVATIVE WAY AND WHAT WE SHOULD USE  
err_from_snr_sw = snr_to_ppm(flx_sw, noise_sw, nint_in)
err_from_snr_lw = snr_to_ppm(flx_lw, noise_lw, nint_in)

#Print median errors FROM ETC
print('Median SW unbinned error (ppm):', int(np.median(err_from_snr_sw)))
print('Median LW unbinned error (ppm):', int(np.median(err_from_snr_lw)))
print('------------------------------')

#----------------WAVELENGTH RANGES-----------------------------

#Set wavelength bounds of interest (taken from JDOX NIRCam page)
if obs_mode_sw == 'F150W2' and obs_mode_lw == 'F322W2':
    wvl_min_sw_1 = 1.07
    wvl_max_sw_1 = 1.66
    wvl_mmin_sw_2 = 1.71
    wvl_max_sw_2 = 2.01
    wvl_min_lw = 2.413
    wvl_max_lw = 4.083 

if obs_mode_sw == 'F200W' and obs_mode_lw == 'F322W2':
    wvl_min_sw_1 = 1.76
    wvl_max_sw_2 = 2.23
    wvl_min_lw = 2.413
    wvl_max_lw = 4.083 

if obs_mode_sw == 'F150W2' and obs_mode_lw == 'F444W':
    wvl_min_sw_1 = 1.01
    wvl_max_sw_1 = 1.23
    wvl_mmin_sw_2 = 1.28
    wvl_max_sw_2 = 1.89
    wvl_min_lw = 3.835
    wvl_max_lw = 5.084


#--------------PLOT--------------------------------

plt.plot(wvl_sw, err_from_snr_sw, label='SW Precision from ETC SNR')
plt.plot(wvl_lw, err_from_snr_lw, label='LW Precision from ETC SNR')

plt.axhline(y=np.median(err_from_snr_sw), xmin=0, xmax=0.3, ls='--', color='grey')
plt.axhline(y=np.median(err_from_snr_lw), xmin=0.5, xmax=1, ls='--', color='grey')

plt.xlabel('Wavelength (microns)')
plt.ylabel('Predicted Error (ppm)')
plt.xlim(wvl_min_sw_1, wvl_max_lw)
plt.ylim(100, np.median(err_from_snr_sw)+200)
plt.legend()
plt.savefig(fig_path+str(planet)+'_'+str(mode)+'_'+str(obs_mode_sw)+'_'+str(obs_mode_lw)+'_'+str(substripe)+'_'+ \
    str(dhs)+'_'+str(ngrps)+'grps_error_estimates.pdf')
plt.show()


#-------------Save----------------------------
ascii.write([wvl_sw, err_from_snr_sw/1e6], names = ['wavelength', 'error_per_bin'],
output = save_path+str(planet)+'_'+str(mode)+'_SW_'+str(obs_mode_sw)+'_'+str(substripe)+'_'+ \
str(dhs)+'_'+str(ngrps)+'_grps.txt', overwrite = True)

ascii.write([wvl_lw, err_from_snr_lw/1e6], names = ['wavelength', 'error_per_bin'],
output = save_path+str(planet)+'_'+str(mode)+'_LW_'+str(obs_mode_lw)+'_'+str(substripe)+'_'+ \
str(dhs)+'_'+str(ngrps)+'_grps.txt', overwrite = True)
