import numpy as np
import os
from astropy.io import fits
from astropy.io import ascii
import astropy.units as u

dir = os.path.abspath('input_test')
object_file = '/final_TM_test.fits'
file = fits.open(dir + object_file)

hdr_o = file[0].header
datat = file[1].data
cols = datat.columns

# Columns to process for invalid values
columns_to_process = [col for col in cols.names if col not in ['quick_object_id', 'ra', 'dec']]

for col in columns_to_process:
    # Convert to float if it's an integer array to allow NaN values
    if datat[col].dtype.kind in 'i':  # integer type
        datat[col] = datat[col].astype(float)
    
    if hasattr(datat[col], "mask"):
        datat[col][datat[col].mask] = np.nan
    w = np.where((datat[col] == -99.0) | (datat[col] == 9999.0) | (datat[col] < 0.) | (datat[col] > 50))
    if len(w[0]) > 0:
        datat[col][w] = np.nan

file.writeto(dir + f"{object_file}_nan.fits", overwrite=True)
file.close()

f = fits.open(dir + f'{object_file}_nan.fits')

hdr = f[0].header
tbdata = f[1].data
cols = tbdata.columns

id_list = tbdata['quick_object_id'][:]
ra = tbdata['ra'][:]
dec = tbdata['dec'][:]

z_1 = np.ones((len(id_list)))
z = z_1 * -1.

def mag_to_flux(mag, err):    
    flux_mJy = (mag * u.ABmag).to('mJy')
    mag_2 = mag + err
    flux_err_mJy = (mag_2 * u.ABmag).to('mJy')
    err_final = np.abs(flux_mJy - flux_err_mJy)
    return flux_mJy, err_final

# Process specific magnitudes
gflux_mJy_delve, gerr_final_delve = mag_to_flux(tbdata['mag_auto_g'][:], tbdata['magerr_auto_g'][:])
rflux_mJy_delve, rerr_final_delve = mag_to_flux(tbdata['mag_auto_r'][:], tbdata['magerr_auto_r'][:])
iflux_mJy_delve, ierr_final_delve = mag_to_flux(tbdata['mag_auto_i'][:], tbdata['magerr_auto_i'][:])
zflux_mJy_delve, zerr_final_delve = mag_to_flux(tbdata['mag_auto_z'][:], tbdata['magerr_auto_z'][:])

Y_ABmag = tbdata['ypetromag'][:] + 0.6
yflux_mJy_vhs, yerr_final_vhs = mag_to_flux(Y_ABmag, tbdata['ypetromagerr'][:])

J_ABmag = tbdata['jpetromag'][:] + 0.916
jflux_mJy_vhs, jerr_final_vhs = mag_to_flux(J_ABmag, tbdata['jpetromagerr'][:])

H_ABmag = tbdata['hpetromag'][:] + 1.366
hflux_mJy_vhs, herr_final_vhs = mag_to_flux(H_ABmag, tbdata['hpetromagerr'][:])

KS_ABmag = tbdata['kspetromag'][:] + 1.827
ksflux_mJy_vhs, kserr_final_vhs = mag_to_flux(KS_ABmag, tbdata['kspetromagerr'][:])

W1_ABmag = tbdata['w1mpro'][:] + 2.699
w1flux_mJy, w1err_final = mag_to_flux(W1_ABmag, tbdata['w1sigmpro'][:])

W2_ABmag = tbdata['w2mpro'][:] + 3.339
w2flux_mJy, w2err_final = mag_to_flux(W2_ABmag, tbdata['w2sigmpro'][:])

ascii.write([id_list, ra, dec, z, 
            gflux_mJy_delve, gerr_final_delve, rflux_mJy_delve, rerr_final_delve, iflux_mJy_delve, ierr_final_delve, zflux_mJy_delve, zerr_final_delve,
            yflux_mJy_vhs, yerr_final_vhs, jflux_mJy_vhs, jerr_final_vhs, hflux_mJy_vhs, herr_final_vhs, ksflux_mJy_vhs, kserr_final_vhs,
            w1flux_mJy, w1err_final, w2flux_mJy, w2err_final], 
            dir + f'{object_file}_fluxes.dat',
            names=['#id', 'ra', 'dec', 'redshift', 
                   'g_prime_delve', 'g_prime_err_delve', 'r_prime_delve', 'r_prime_err_delve', 'i_prime_delve', 'i_prime_err_delve', 'z_prime_delve', 'z_prime_err_delve',
                   'vista.vircam.Y_vhs', 'vista.vircam.Y_err_vhs', 'vista.vircam.J_vhs', 'vista.vircam.J_err_vhs', 'vista.vircam.H_vhs', 'vista.vircam.H_err_vhs', 'vista.vircam.Ks_vhs', 'vista.vircam.Ks_err_vhs',
                   'WISE1', 'WISE1_err', 'WISE2', 'WISE2_err'],  
                   overwrite=True)
