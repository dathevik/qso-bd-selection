import shutil

import numpy as np
from decimal import Decimal
import matplotlib.pyplot as plt
import os
from astropy.io import ascii
from astropy.table import Table, QTable

from synphot import SourceSpectrum, etau_madau
from synphot.models import Empirical1D


# ------------------------------- Define functions ---------------------------------

# function for a_scale which is a scaling factor and takes as an input of observed flux, observed flux error of the object and the template model
def a_scale(vec_flux_obs, vec_fluxe_obs, vec_flux_model):
    # ---- Obtain scaling factor for Chi2
    a = np.sum((vec_flux_obs * vec_flux_model) / (vec_fluxe_obs) ** 2) / np.sum(
        (vec_flux_model) ** 2 / (vec_fluxe_obs) ** 2)
    return a


# function for calculation of chi2 statistical parameter and takes as an input of observed flux, observed flux error of the object and the template model and a scaling factor
def chi2_calc(vec_flux_obs, vec_fluxe_obs, vec_flux_model, a):
    # ---- Obtain scaling factor for Chi2
    chi2 = np.sum((vec_flux_obs - a * vec_flux_model) ** 2 / (vec_fluxe_obs) ** 2)
    return chi2


# function for converting cgs parameters to mJy as the whole script works with mJy units of fluxes
def flux_from_cgs_to_mJy(flux_cgs, ll):
    # ---- Obtain flux in mJy from flux in erg/s/cm2/A and considering effective wavelength of the filter in AA
    flux_mJy = (3.33564095E+04 * flux_cgs * ll ** 2) * 1E-3
    return flux_mJy


# ------------------------------- Main code ---------------------------------
# I: Define the path and print the data in the created directory
BD_temp_path = os.path.abspath('input_test/BDRA_fluxes_mJy_202301.dat')
data_bd_temp = ascii.read(BD_temp_path)
mdRA_path = os.path.abspath('input_test/lrt_kc04_MS.dat')
bdRA_path = os.path.abspath('input_test/lrt_a07_BDs.dat')
QSO_temp_path = os.path.abspath('input_test/Selsing+Matt_temp.dat')
data_qso_temp = ascii.read(QSO_temp_path)
QSO_spec_path = os.path.abspath('input_test/Selsing2015.dat')
data_qso_spec = ascii.read(QSO_spec_path)
output_folder = os.path.join(os.path.abspath(''), 'output_test')
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)
output_file = os.path.abspath('output_test/results.csv')


# Input the information
ls_id = input("Input the name of object").strip()
ls_id = str(ls_id)
input_file_name = input("Input the name of catalog ")
input_objects_path = os.path.abspath('input_test/'+input_file_name)
list_type = input("Input the type of list (1 for QSO, 2 for BD, 0 for test): ").strip()
list_type = int(list_type)

print("READING INPUT CATALOG:", input_objects_path, "FILE")
# Read the file with latin-1 encoding to handle special characters
data_fil = ascii.read(input_objects_path, encoding='latin-1')

id = data_fil.columns[0]
print(id)

# Initialize redshift variable
object_redshift = None

# Only read redshift column for QSO list type
if list_type == 1:
    if 'redshift' in data_fil.colnames:
        redshift = data_fil['redshift']
    else:
        print("Warning: 'redshift' column not found in the input file")
        redshift = None

# Only read optical_type column for BD list type
if list_type == 2:
    if 'optical_type' in data_fil.colnames:
        optical_type = data_fil['optical_type']
    else:
        print("Warning: 'optical_type' column not found in the input file")
        optical_type = None

# IIa : Make the vector arrays for Brown Dwarf templates
BD_g_decam = data_bd_temp.columns[1]
BD_r_decam = data_bd_temp.columns[2]
BD_i_decam = data_bd_temp.columns[3]
BD_z_decam = data_bd_temp.columns[4]
BD_Y_vhs = data_bd_temp.columns[5]
BD_J_vhs = data_bd_temp.columns[6]
BD_H_vhs = data_bd_temp.columns[7]
BD_K_vhs = data_bd_temp.columns[8]
BD_W1 = data_bd_temp.columns[9]
BD_W2 = data_bd_temp.columns[10]
#
# # Read BD template rows in a loop
BD_vec_flux = []
for i in range(len(data_bd_temp)):
    vec_flux = [BD_g_decam[i], BD_r_decam[i], BD_i_decam[i], BD_z_decam[i], BD_Y_vhs[i], BD_J_vhs[i], BD_H_vhs[i],
                BD_K_vhs[i], BD_W1[i], BD_W2[i]]
    BD_vec_flux.append(vec_flux)

# # Make list with BDs type and print first template as an example
BD_all_vec_flux = np.array(BD_vec_flux)
vec_flux_model_BD = BD_all_vec_flux.astype(float)
BD_type_vec = ["BD1", "BD2", "BD3", "BD4", "MD1", "MD2", "MD3"]

 # -------- Read RAssef BD and MD templates for plotting and comparing the results
mdRA_data = Table(ascii.read(mdRA_path))
bdRA_data = Table(ascii.read(bdRA_path))
bdRA_l1 = bdRA_data.columns[0] * 10 ** 4  # --AA
bdRA_l2 = bdRA_data.columns[1] * 10 ** 4  # --AA
bdRA_l = bdRA_l2 + (bdRA_l1 - bdRA_l2) / 2.  # ---bin center in the  BD template
mdRA_l1 = mdRA_data.columns[0] * 10 ** 4  # --AA
mdRA_l2 = mdRA_data.columns[1] * 10 ** 4  # --AA
mdRA_l = mdRA_l2 + (mdRA_l1 - mdRA_l2) / 2.  # ---bin center in the MD template

# # ----- Fluxes vectors in nuFnu (Hz Jy)
bdRA_nf1 = bdRA_data.columns[2]
bdRA_nf2 = bdRA_data.columns[3]
bdRA_nf3 = bdRA_data.columns[4]
bdRA_nf4 = bdRA_data.columns[5]
mdRA_nf1 = mdRA_data.columns[33]
mdRA_nf2 = mdRA_data.columns[34]
mdRA_nf3 = mdRA_data.columns[35]
#
# # ----- Convert fluxes from nuFnu (Hz Jy) in (erg/s/cm2/A)
bdRA_f1 = bdRA_nf1
bdRA_f2 = bdRA_nf2
bdRA_f3 = bdRA_nf3
bdRA_f4 = bdRA_nf4
mdRA_f1 = mdRA_nf1
mdRA_f2 = mdRA_nf2
mdRA_f3 = mdRA_nf3
bd_md_array = QTable([bdRA_f1, bdRA_f2, bdRA_f3, bdRA_f4, mdRA_f1, mdRA_f2, mdRA_f3], names=(BD_type_vec))

# # -----Define central wavelength for plot it later for check
ll_g = 4798.3527009231575  # Angstrom
ll_r = 6407.493598028656  # Angstrom
ll_i = 7802.488114833454  # Angstrom
ll_z = 9144.625340022629  # Angstrom
ll_J = 12325.125694338809  # Angstrom
ll_Y = 10201.359507821942  # Angstrom
ll_H = 16473.95843628733  # Angstrom
ll_K = 22045.772662096875  # Angstrom
ll_W1 = 33791.878497259444  # Angstrom
ll_W2 = 46292.93969033106  # Angstrom
ll_vec = np.array([ll_g, ll_r, ll_i, ll_z, ll_Y, ll_J, ll_H, ll_K, ll_W1, ll_W2])
#
# # IIb : Make vector arrays for Quasars templates
QSO_z_vec = data_qso_temp.columns[1]
QSO_g_decam = data_qso_temp.columns[2]
QSO_r_decam = data_qso_temp.columns[3]
QSO_i_decam = data_qso_temp.columns[4]
QSO_z_decam = data_qso_temp.columns[5]
QSO_Y_vhs = data_qso_temp.columns[6]
QSO_J_vhs = data_qso_temp.columns[7]
QSO_H_vhs = data_qso_temp.columns[8]
QSO_K_vhs = data_qso_temp.columns[9]
QSO_W1 = data_qso_temp.columns[10]
QSO_W2 = data_qso_temp.columns[11]
QSO_temp_type = data_qso_temp.columns[0]
QSO_temp_emline = data_qso_temp.columns[12]
QSO_temp_ebv = data_qso_temp.columns[13]
#
# # ----- Read Column of Redshift in a loop and print the first template as an example
QSO_vec_flux = []

for i in range(len(data_qso_temp)):
    vec_flux = [QSO_g_decam[i], QSO_r_decam[i], QSO_i_decam[i], QSO_z_decam[i], QSO_Y_vhs[i], QSO_J_vhs[i],
                QSO_H_vhs[i], QSO_K_vhs[i], QSO_W1[i], QSO_W2[i]]
    QSO_vec_flux.append(vec_flux)
QSO_all_vec_flux = np.array(QSO_vec_flux)
vec_flux_model_QSO = QSO_all_vec_flux.astype(float)

#
# IVa : Make array for central fluxes of the bands(each of them)
vec_flux = data_fil[
    ["g_prime_delve", "r_prime_delve", "i_prime_delve", "z_prime_delve", "vista.vircam.Y_vhs", "vista.vircam.J_vhs",
     "vista.vircam.H_vhs", "vista.vircam.Ks_vhs", "WISE1", "WISE2"]].copy()
# print(vec_flux)

# IVb : Make array for errors of the bands(each of them)
vec_fluxe = data_fil[
    ["g_prime_err_delve", "r_prime_err_delve", "i_prime_err_delve", "z_prime_err_delve", "vista.vircam.Y_err_vhs",
     "vista.vircam.J_err_vhs", "vista.vircam.H_err_vhs", "vista.vircam.Ks_err_vhs", "WISE1_err", "WISE2_err"]].copy()
# print(vec_fluxe)

object_name = []
data = []
object_redshift = []
vec_flux_obs = []
vec_fluxe_obs = []

print(f"CHECKING {ls_id} OBJECT IN {input_file_name} CATALOG ")

for i in range(len(id)):
    if ls_id == str(id[i]):
        print(f"{ls_id} OBJECT is FOUND in {input_file_name} CATALOG")
        object_name.append(id[i])
        data.append(data_fil[i])
        if list_type == 1 and redshift is not None:  # Only store redshift for QSO list
            object_redshift = float(redshift[i])  # Convert to float to ensure proper formatting
            print(f"Object redshift: {object_redshift}")
        x = vec_flux[i]
        y = vec_fluxe[i]
        for j in range(len(x)):
                if y[j] < 0.1 * x[j]:
                    y[j] = 0.1 * x[j]
                if x[j] <= 10 ** -30:
                    x[j] = np.nan
                vec_flux_obs.append(x[j])
                vec_fluxe_obs.append(y[j])

        BD_Chi2_array = []
        QSO_Chi2_array = []
        R_Chi = []
        a_BD_array = []
        a_QSO_array = []

        vec_flux_row = vec_flux[i]
        vec_fluxe_row = vec_fluxe[i]
        vec_flux_obs0 = np.array(
            [vec_flux_row[0], vec_flux_row[1], vec_flux_row[2], vec_flux_row[3], vec_flux_row[4], vec_flux_row[5],
             vec_flux_row[6], vec_flux_row[7], vec_flux_row[8], vec_flux_row[9]])
        vec_fluxe_obs0 = np.array(
            [vec_fluxe_row[0], vec_fluxe_row[1], vec_fluxe_row[2], vec_fluxe_row[3], vec_fluxe_row[4], vec_fluxe_row[5],
             vec_fluxe_row[6], vec_fluxe_row[7], vec_fluxe_row[8], vec_fluxe_row[9]])

        # ---- mask the nan value, which means it will not be used for the calculation
        mask_nan = ~np.isnan(vec_flux_obs0)
        vec_flux_obs = vec_flux_obs0[mask_nan]
        vec_fluxe_obs = vec_fluxe_obs0[mask_nan]
        ll_vec_best = ll_vec[mask_nan]
        # print(vec_flux_obs)

        # Va : Calculate scaling factor and Chi2 for BDs
        for j in range(len(data_bd_temp)):
            vec_flux_model_BD_j0 = vec_flux_model_BD[j]
            vec_flux_model_BD_j = vec_flux_model_BD_j0[mask_nan]
            a_BD = a_scale(vec_flux_obs, vec_fluxe_obs, vec_flux_model_BD_j)
            a_BD = a_BD.astype(float)
            Chi2_BD = chi2_calc(vec_flux_obs, vec_fluxe_obs, vec_flux_model_BD_j, a_BD)
            BD_Chi2_array.append(Chi2_BD)
            a_BD_array.append(a_BD)

        # Vb : Calculate scaling factor and Chi2 for QSOs
        for k in range(len(data_qso_temp)):
            vec_flux_model_QSO_k0 = vec_flux_model_QSO[k]
            vec_flux_model_QSO_k = vec_flux_model_QSO_k0[mask_nan]
            a_QSO = a_scale(vec_flux_obs, vec_fluxe_obs, vec_flux_model_QSO_k)
            a_QSO = a_QSO.astype(float)
            Chi2_QSO = chi2_calc(vec_flux_obs, vec_fluxe_obs, vec_flux_model_QSO_k, a_QSO)
            QSO_Chi2_array.append(Chi2_QSO)
            a_QSO_array.append(a_QSO)

        # VI : Calculate and print out the template with the best (lowest) Chi2 for QSOs and BDs
        # VIa : Calculate the lowest value for the Chi2 for the BDs
        BD_Chi2_min = np.min(BD_Chi2_array)  # -- minimum value
        BD_Chi2_min_ind = np.argmin(BD_Chi2_array)  # -- position in array
        a_BD_best = a_BD_array[BD_Chi2_min_ind]  # -- corresponding scaling factor of the best chi2
        BD_Chi2_min_temp = BD_type_vec[BD_Chi2_min_ind]  # -- corresponding best template
        vec_flux_model_BD_best = BD_all_vec_flux[BD_Chi2_min_ind][mask_nan] * a_BD_best


        # VIb. Calculate the lowest value for the Chi2 for the QSOs
        QSO_Chi2_min = np.min(QSO_Chi2_array)  # -- minimum value
        QSO_Chi2_min_ind = np.argmin(QSO_Chi2_array)  # -- po11036800232067sition in array
        a_QSO_best = a_QSO_array[QSO_Chi2_min_ind]  # -- corresponding scaling factor of the best chi2
        QSO_Chi2_min_z = QSO_z_vec[QSO_Chi2_min_ind]  # -- corresponding best template -> best redshift
        vec_flux_model_QSO_best = QSO_all_vec_flux[QSO_Chi2_min_ind][mask_nan] * a_QSO_best
        temp_ebv_best = QSO_temp_ebv[QSO_Chi2_min_ind]
        temp_type_best = QSO_temp_type[QSO_Chi2_min_ind]
        temp_emline_best = QSO_temp_emline[QSO_Chi2_min_ind]
    # flux_stripe82_test_03.dat
        # ----Calculate the ratio of chi2 of each QSO and BD templates
        R = QSO_Chi2_min / BD_Chi2_min
        R_Chi.append(R)
        print("Best Chi2 QSO:", QSO_Chi2_min)
        print("Best Chi2 BD:", BD_Chi2_min)
        print("Chi2 Ratio (Chi2_QSO/Chi2_BD) is:", R)

        # ---- parameters for plotting in a loop using Matt or S15 templates
        z = QSO_Chi2_min_z
        if temp_type_best == "S15":
            w_rest = np.arange(900., 75000, 10)
            spec_f = np.array(data_qso_spec)
            wave = data_qso_spec.columns[0]
            flux = data_qso_spec.columns[1]
            sp = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux, keep_neg=True)
            sp_z_chi2 = SourceSpectrum(sp.model, z=z)
            extcurve_best = etau_madau(w_rest, z)
            sp_ext_best = sp_z_chi2 * extcurve_best
            flux_mJY_best = flux_from_cgs_to_mJy(sp_ext_best(w_rest), w_rest) * a_QSO_best
            if BD_Chi2_min_ind > 3:
                bdRA_l_fin = mdRA_l
            else:
                bdRA_l_fin = bdRA_l

            bdRA_f_best = (bd_md_array.columns[BD_Chi2_min_ind] / bdRA_l_fin) * a_BD_best
            bdRA_fmJy_best = flux_from_cgs_to_mJy(bdRA_f_best, bdRA_l_fin)
            plt.figure(figsize=(12, 6), dpi=300)
            plt.rcParams.update({'font.size': 18})
            plt.rcParams['axes.linewidth'] = 1.5
            plt.rcParams['xtick.major.width'] = 1.5
            plt.rcParams['ytick.major.width'] = 1.5
            plt.rcParams['xtick.minor.width'] = 1.0
            plt.rcParams['ytick.minor.width'] = 1.0
            plt.rcParams['xtick.major.size'] = 8
            plt.rcParams['ytick.major.size'] = 8
            plt.rcParams['xtick.minor.size'] = 4
            plt.rcParams['ytick.minor.size'] = 4

            # Set y-axis to log scale
            plt.yscale('log')

            plt.scatter(ll_vec_best, vec_flux_model_BD_best, color="#008000", alpha=0.6, marker='^',
                        label=f"Best BD template ({BD_Chi2_min_temp}, χ² = {np.around(BD_Chi2_min, decimals=2)})",
                        zorder=5, facecolor='none', s=100, linewidth=2)

            # Only show photo_z for QSO lists
            if list_type == 1:  # QSO
                plt.scatter(ll_vec_best, vec_flux_model_QSO_best, color="#ff8c69",
                            label=f"Best QSO template (χ² = {np.around(QSO_Chi2_min, decimals=2)}, photo_z = {np.around(QSO_Chi2_min_z, decimals=2)})",
                            zorder=5, facecolor='none', s=100, linewidth=2)
            else:  # BD or test
                plt.scatter(ll_vec_best, vec_flux_model_QSO_best, color="#ff8c69",
                            label=f"Best QSO template (χ² = {np.around(QSO_Chi2_min, decimals=2)})",
                            zorder=5, facecolor='none', s=100, linewidth=2)
            
            # Modify the label based on list_type
            if list_type == 1:  # QSO
                plt.errorbar(ll_vec_best, vec_flux_obs, yerr=vec_fluxe_obs, fmt='o', color='#ff6b6b',
                             label=f'Initial flux of the quasar (spec_z = {object_redshift})',
                             markersize=8, capsize=4, capthick=2, elinewidth=2)
            elif list_type == 2:  # BD
                plt.errorbar(ll_vec_best, vec_flux_obs, yerr=vec_fluxe_obs, fmt='o', color='#008000',
                             label=f'Initial flux of the BD (SpecT = {data_fil["optical_type"][i]})',
                             markersize=8, capsize=4, capthick=2, elinewidth=2)
            else:  # test
                plt.errorbar(ll_vec_best, vec_flux_obs, yerr=vec_fluxe_obs, fmt='o', color='#FCD12A',
                             label='Initial flux of the object',
                             markersize=8, capsize=4, capthick=2, elinewidth=2)
            
            plt.plot(w_rest, flux_mJY_best, color="#ff6b6b", alpha=0.3,
                     label='QSO model spectrum', 
                     zorder=0, linewidth=2)
            plt.plot(bdRA_l_fin, bdRA_fmJy_best, color="#008000", alpha=0.3,
                     label='BD model spectrum',
                     zorder=0, linewidth=2)
 
            plt.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=18)
            plt.xlabel("Central Wavelength (Å)", fontsize=18)
            plt.ylabel("Log Flux (mJy)", fontsize=18)
            plt.xlim(0, 55000)
            # Set y-axis limits with upper limit at 1 mJy
            ymin = min(np.min(vec_flux_obs), np.min(vec_flux_model_BD_best), np.min(vec_flux_model_QSO_best)) * 0.8
            plt.ylim(ymin, 1)
            plt.grid(True, linestyle='--', alpha=0.3)
            plt.tight_layout()
            plt.savefig(f"output_test/{ls_id}_sed.png", dpi=300, bbox_inches='tight')
            plt.show()
            plt.close()

        else:
            w_rest = np.arange(900., 18900, 10)
            z = round(QSO_Chi2_min_z * 100)
            ebv_dec = round((Decimal(temp_ebv_best)), 2)
            spec_ebv = str(ebv_dec).replace(".", "")
            spec_file = (f"input_test/seds2/qsogen_sed_emlines{temp_emline_best}_ebv+{spec_ebv}_z+{z}.dat")
            spec_open = os.path.abspath(spec_file)
            data_spec_Matt = ascii.read(spec_file)
            wave = data_spec_Matt.columns[0]
            flux = data_spec_Matt.columns[1]
            flux_mJY_best = flux_from_cgs_to_mJy(flux, wave) * a_QSO_best
            if BD_Chi2_min_ind > 3:
                bdRA_l_fin = mdRA_l
            else:
                bdRA_l_fin = bdRA_l

            bdRA_f_best = (bd_md_array.columns[BD_Chi2_min_ind] / bdRA_l_fin) * a_BD_best
            bdRA_fmJy_best = flux_from_cgs_to_mJy(bdRA_f_best, bdRA_l_fin)
            plt.figure(figsize=(12, 6), dpi=300)
            plt.rcParams.update({'font.size': 18})
            plt.rcParams['axes.linewidth'] = 1.5
            plt.rcParams['xtick.major.width'] = 1.5
            plt.rcParams['ytick.major.width'] = 1.5
            plt.rcParams['xtick.minor.width'] = 1.0
            plt.rcParams['ytick.minor.width'] = 1.0
            plt.rcParams['xtick.major.size'] = 8
            plt.rcParams['ytick.major.size'] = 8
            plt.rcParams['xtick.minor.size'] = 4
            plt.rcParams['ytick.minor.size'] = 4

            # Set y-axis to log scale
            plt.yscale('log')

            plt.scatter(ll_vec_best, vec_flux_model_BD_best, color="#008000", alpha=0.6, marker='^',
                        label=f"Best BD template ({BD_Chi2_min_temp}, χ² = {np.around(BD_Chi2_min, decimals=2)})",
                        zorder=5, facecolor='none', s=100, linewidth=2)

            # Only show photo_z for QSO lists
            if list_type == 1:  # QSO
                plt.scatter(ll_vec_best, vec_flux_model_QSO_best, color="#ff8c69",
                            label=f"Best QSO template (χ² = {np.around(QSO_Chi2_min, decimals=2)}, photo_z = {QSO_Chi2_min_z})",
                            zorder=5, facecolor='none', s=100, linewidth=2)
            else:  # BD or test
                plt.scatter(ll_vec_best, vec_flux_model_QSO_best, color="#ff8c69",
                            label=f"Best QSO template (χ² = {np.around(QSO_Chi2_min, decimals=2)})",
                            zorder=5, facecolor='none', s=100, linewidth=2)
            
            # Modify the label based on list_type
            if list_type == 1:  # QSO
                plt.errorbar(ll_vec_best, vec_flux_obs, yerr=vec_fluxe_obs, fmt='o', color='#FCD12A',
                             label=f'Initial flux of the quasar (spec_z = {object_redshift})',
                             markersize=8, capsize=4, capthick=2, elinewidth=2)
            elif list_type == 2:  # BD
                plt.errorbar(ll_vec_best, vec_flux_obs, yerr=vec_fluxe_obs, fmt='o', color='#FCD12A',
                             label=f'Initial flux of the BD (SpecT = {data_fil["optical_type"][i]})',
                             markersize=8, capsize=4, capthick=2, elinewidth=2)
            else:  # test
                plt.errorbar(ll_vec_best, vec_flux_obs, yerr=vec_fluxe_obs, fmt='o', color='#FCD12A',
                             label='Initial flux of the object',
                             markersize=8, capsize=4, capthick=2, elinewidth=2)
            
            plt.plot(wave, flux_mJY_best*1e6, color="#ff8c69", alpha=0.3,
                     label='QSO model spectrum', 
                     zorder=0, linewidth=2)
            plt.plot(bdRA_l_fin, bdRA_fmJy_best, color="#008000", alpha=0.3,
                     label='BD model spectrum',
                     zorder=0, linewidth=2)

            plt.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=14)
            plt.xlabel("Central Wavelength (Å)", fontsize=18)
            plt.ylabel("Log Flux (mJy)", fontsize=18)
            plt.xlim(0, 55000)
            # Set y-axis limits with upper limit at 1 mJy
            ymin = min(np.min(vec_flux_obs), np.min(vec_flux_model_BD_best), np.min(vec_flux_model_QSO_best)) * 0.8
            plt.ylim(ymin, 1)
            plt.grid(True, linestyle='--', alpha=0.3)
            plt.tight_layout()
            plt.savefig(f"output_test/{ls_id}_sed.png", dpi=300, bbox_inches='tight')
            plt.show()
            plt.close()

