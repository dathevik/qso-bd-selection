import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.units as u

# Set the style for publication-quality plots
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['figure.figsize'] = [12, 10]
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linestyle'] = '--'

#Convert fluxes to magnitudes
def flux_to_mag(flux, flux_zero):
    return (-2.5 * np.log10(flux / flux_zero)) * u.ABmag

# Zero-point flux values in mJy
flux_zero_r = 3.631
flux_zero_I = 3.631
flux_zero_z = 3.631

# Define colorblind-friendly colors
template_colors = ['#0072B2', '#D55E00', '#009E73']  # Blue, Orange, Green for templates
selsing_color = '#CC79A7'  # Pink
qso_color = '#ff6b6b'      # Coral Red for F23 sample
bd_color = '#008000'       # Green for BD sample
fan_color = '#009E73'      # Green
cut_color = 'grey'      # Pink
yang_color = '#4169E1'     # Royal Blue for Yang sample
milliquas_color = '#FF4500'  # Orange Red for Milliquas sample

# Create the main plot
fig, ax = plt.subplots(figsize=(12, 10))

# Plot QSO templates
matt_filelist = ['qsogen_mags_emlines2.0_ebv+025.dat', 'qsogen_mags_emlines-2.0_ebv+010.dat', 'qsogen_mags_emlines0.0_ebv+000.dat']
label_temp = ['QSO template, emlines: 2, E(B-V): 0.25 (Temple+2021)', 
              'QSO template, emlines: -2, E(B-V): 0.1 (Temple+2021)', 
              'QSO template, emlines: 0, E(B-V): 0 (Temple+2021)']

for file, label, color in zip(matt_filelist, label_temp, template_colors):
    template_path = os.path.abspath(f'{file}')
    data_temp = ascii.read(template_path)
    r_decam_flux = data_temp.columns[2]
    i_decam_flux = data_temp.columns[3]
    redshift = data_temp.columns[0]
    ri = r_decam_flux - i_decam_flux
    plt.plot(redshift, ri, label=label, color=color, linewidth=2, alpha=0.7)

# Plot Selsing templates
sels_path = os.path.abspath('Selsing2015_fluxes_mJy_202301.dat')
data_sels = ascii.read(sels_path)
redshift_sels = data_sels.columns[1]
r_decam_flux_sels = data_sels.columns[3]
i_decam_flux_sels = data_sels.columns[4]
r_decam_mag_sels = flux_to_mag(r_decam_flux_sels, flux_zero_r)
i_decam_mag_sels = flux_to_mag(i_decam_flux_sels, flux_zero_I)
ri_sels = r_decam_mag_sels - i_decam_mag_sels
plt.plot(redshift_sels, ri_sels, label='QSO template (Selsing+2016)', color=selsing_color, linewidth=2, alpha=0.7)

# Plot known quasars (Fan sample)
araa_path = os.path.abspath('araa_delve_query.dat')
data_araa = ascii.read(araa_path)
r_mag_qso = data_araa.columns[0]
i_mag_qso = data_araa.columns[1]
redshift_qso = data_araa.columns[9]
ri_mag_qso = r_mag_qso - i_mag_qso
plt.scatter(redshift_qso, ri_mag_qso, s=250, color=qso_color, alpha=0.6, label='F23 Sample', marker='*', edgecolor='black')

# Plot Yang 2023 quasar sample
yang_path = os.path.abspath('yang_results.dat')
try:
    data_yang = ascii.read(yang_path)
    r_mag_yang = data_yang['mag_auto_r']
    i_mag_yang = data_yang['mag_auto_i']
    redshift_yang = data_yang['QSO_z']
    ri_mag_yang = r_mag_yang - i_mag_yang
#    plt.scatter(redshift_yang, ri_mag_yang, s=60, color=yang_color, alpha=0.5, label='Y23 Sample', marker='s', edgecolor='black')
except Exception as e:
    print(f"Warning: Could not load Yang 2023 data. Error: {str(e)}")
    print("Available columns:", data_yang.colnames if 'data_yang' in locals() else "No data loaded")

# Plot Milliquas quasar sample
milliquas_path = os.path.abspath('milliquas_results.dat')
try:
    data_milliquas = ascii.read(milliquas_path)
    r_mag_milliquas = data_milliquas['mag_auto_r']
    i_mag_milliquas = data_milliquas['mag_auto_i']
    redshift_milliquas = data_milliquas['QSO_z']
    ri_mag_milliquas = r_mag_milliquas - i_mag_milliquas
#    plt.scatter(redshift_milliquas, ri_mag_milliquas, s=4, color=milliquas_color, alpha=0.4, label='Fl23 Sample', marker='o', edgecolor='black')
except Exception as e:
    print(f"Warning: Could not load Milliquas data. Error: {str(e)}")
    print("Available columns:", data_milliquas.colnames if 'data_milliquas' in locals() else "No data loaded")

# Plot Fan et al., 2023 quasar sample
# Assuming you have a file with Fan et al., 2023 data
fan_path = os.path.abspath('fan_et_al_2023_quasars.dat')
try:
    data_fan = ascii.read(fan_path)
    r_mag_fan = data_fan.columns[0]  # Adjust column indices as needed
    i_mag_fan = data_fan.columns[1]
    redshift_fan = data_fan.columns[2]
    ri_mag_fan = r_mag_fan - i_mag_fan
    plt.scatter(redshift_fan, ri_mag_fan, s=15, color=fan_color, alpha=0.6, label='F23 Sample', marker='o')
except:
    print("Warning: Could not load Fan et al., 2023 data. File may not exist or have different format.")

# Plot spectral confirmation data (green stars, similar to statistics_qso_bd.py)
current_dir = os.path.dirname(os.path.abspath(__file__))
spectral_file = None
possible_paths = [
    os.path.join(current_dir, "final_TM_spectral_confirmation.dat"),
    os.path.join(current_dir, "..", "..", "output_catalogs", "results_catalogs", "final_TM_spectral_confirmation.dat"),
    '/Users/dathev/PhDProjects/high-z-qso-selection/output_catalogs/results_catalogs/final_TM_spectral_confirmation.dat',
    os.path.join(current_dir, "..", "results", "final_TM_spectral_confirmation.dat"),
]

for path in possible_paths:
    abs_path = os.path.abspath(path)
    if os.path.exists(abs_path):
        spectral_file = abs_path
        print(f"Found spectral confirmation file at: {abs_path}")
        break

if spectral_file:
    try:
        data_spectral = ascii.read(spectral_file)
        print(f"Successfully loaded spectral confirmation file with {len(data_spectral)} sources")
        
        # Check for column names (could be QSO_z or redshift, mag_r_delve or mag_auto_r, etc.)
        if 'QSO_z' in data_spectral.colnames:
            redshift_spectral = np.array(data_spectral['QSO_z'], dtype=float)
        elif 'redshift' in data_spectral.colnames:
            redshift_spectral = np.array(data_spectral['redshift'], dtype=float)
        else:
            print("Warning: Could not find redshift column in spectral confirmation file.")
            print("Available columns:", data_spectral.colnames)
            redshift_spectral = None
        
        # Check for r and i magnitude columns
        if 'mag_r_delve' in data_spectral.colnames and 'mag_i_delve' in data_spectral.colnames:
            r_mag_spectral = np.array(data_spectral['mag_r_delve'], dtype=float)
            i_mag_spectral = np.array(data_spectral['mag_i_delve'], dtype=float)
        elif 'mag_auto_r' in data_spectral.colnames and 'mag_auto_i' in data_spectral.colnames:
            r_mag_spectral = np.array(data_spectral['mag_auto_r'], dtype=float)
            i_mag_spectral = np.array(data_spectral['mag_auto_i'], dtype=float)
        else:
            print("Warning: Could not find r and i magnitude columns in spectral confirmation file.")
            print("Available columns:", data_spectral.colnames)
            r_mag_spectral = None
            i_mag_spectral = None
        
        if redshift_spectral is not None and r_mag_spectral is not None and i_mag_spectral is not None:
            # Filter out invalid values (99.0, NaN, negative, or out of range)
            valid_mask = (
                (redshift_spectral > 0) & (redshift_spectral < 10) &
                (r_mag_spectral > 0) & (r_mag_spectral < 30) &
                (i_mag_spectral > 0) & (i_mag_spectral < 30) &
                (r_mag_spectral != 99.0) & (i_mag_spectral != 99.0) &
                ~np.isnan(redshift_spectral) & ~np.isnan(r_mag_spectral) & ~np.isnan(i_mag_spectral)
            )
            
            redshift_spectral_valid = redshift_spectral[valid_mask]
            r_mag_spectral_valid = r_mag_spectral[valid_mask]
            i_mag_spectral_valid = i_mag_spectral[valid_mask]
            
            if len(redshift_spectral_valid) > 0:
                ri_mag_spectral = r_mag_spectral_valid - i_mag_spectral_valid
                print(f"Plotting {len(redshift_spectral_valid)} confirmed quasars")
                # Plot as hotpink stars matching statistics_qso_bd.py style
                plt.scatter(redshift_spectral_valid, ri_mag_spectral, 
                           alpha=0.9, label='confirmed quasars (this work)', color='hotpink', s=400, 
                           marker='*', edgecolor='black', linewidths=1.0, zorder=10)
            else:
                print("Warning: No valid data points found in spectral confirmation file after filtering")
    except Exception as e:
        print(f"Warning: Could not load spectral confirmation data. Error: {str(e)}")
        import traceback
        traceback.print_exc()
else:
    print("Warning: final_TM_spectral_confirmation.dat not found in any expected location")

# Add color cut lines
# Replace the two separate lines with a rectangle for the selection region
plt.axvspan(4.5, 6.0, ymin=1.3/4.5, ymax=1.0, color=cut_color, alpha=0.2, label='Sources selected')

# Add dashed black frame around the selection region
plt.plot([4.5, 6.0, 6.0, 4.5, 4.5], [1.3, 1.3, 4.5, 4.5, 1.3], 'k--', linewidth=2)

# Add grid for better readability
ax.grid(True, linestyle='--', alpha=0.3)

# Improve axis labels
plt.xlabel('Redshift', fontsize=20)
plt.ylabel('r$_{delve}$ - i$_{delve}$', fontsize=20)

# Set plot limits
plt.xlim(4.5, 7)
plt.ylim(0, 4.5)

# Improve legend - position it higher and to the left to avoid covering data points
plt.legend(fontsize=16, loc='upper left', bbox_to_anchor=(0.0, 1.01), frameon=True, framealpha=0.9, ncol=1)

# Add minor grid lines
ax.minorticks_on()
ax.grid(True, which='minor', linestyle=':', alpha=0.2)

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Show the plot without saving
#plt.show()
plt.savefig('plots/color_cut_fan.png', dpi=300)