# NIRCam_DHS_precision_estimates
Estimates precision for JWST NIRCam DHS observations (for both the short-wavelength and long-wavelength) using the S/N estimates from ETC.
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
