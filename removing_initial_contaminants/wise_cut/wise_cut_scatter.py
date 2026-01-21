import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import random

# Set random seed for reproducibility
random.seed(42)

# Function to sample data
def sample_data(data, n_samples):
    if len(data) > n_samples:
        indices = random.sample(range(len(data)), n_samples)
        return data[indices]
    return data

# Read data
data_qso_araa = ascii.read('QSO_wise_spec.dat')
data_qso_milliquas = ascii.read('QSO_milliquas_wise_clean.dat')
data_qso_yang = ascii.read('QSO_yang_wise.dat')
data_bd_wise = ascii.read('BD_wise_SpT.dat', encoding='latin-1')
data_bd_temp_L = ascii.read("L_type_wise.txt")
data_bd_temp_M = ascii.read("M_type_wise.txt")
data_bd_temp_T = ascii.read("T_type_wise.txt")

# Function to convert column to float, replacing invalid entries with NaN
def convert_to_numeric(column):
    cleaned_column = [str(x).replace(',', '') for x in column]  # Remove commas
    return np.array(cleaned_column, dtype=float)

# Quasars ARAA
W1_qso_araa = data_qso_araa['w1mpro'] + 2.699
W2_qso_araa = data_qso_araa['w2mpro'] + 3.399

# Quasars Million Quasars Catalog v8
W1_qso_milliquas = data_qso_milliquas['w1mpro_1'] + 2.699
W2_qso_milliquas = data_qso_milliquas['w2mpro_1'] + 3.399

#Quasars Yang 2023
W1_qso_yang = data_qso_yang["w1mpro"] + 2.699
W2_qso_yang = data_qso_yang["w2mpro"] + 3.399

# Brown Dwarfs list (Bañados+2016)
W1_bd = data_bd_wise['w1mpro'] + 2.699
W2_bd = data_bd_wise['w2mpro'] + 3.399
Sp_T_bd = data_bd_wise['optical_type']

# Brown Dwarfs Empirical templates
Sp_T_bd_temp_L = data_bd_temp_L.columns[0]
Sp_T_bd_temp_M = data_bd_temp_M.columns[0]
Sp_T_bd_temp_T = data_bd_temp_T.columns[0]
W1_bd_temp_L = convert_to_numeric(data_bd_temp_L.columns[1]) + 2.699
W2_bd_temp_L = convert_to_numeric(data_bd_temp_L.columns[2]) + 3.399
W1_bd_temp_M = convert_to_numeric(data_bd_temp_M.columns[1]) + 2.699
W2_bd_temp_M = convert_to_numeric(data_bd_temp_M.columns[2]) + 3.399
W1_bd_temp_T = convert_to_numeric(data_bd_temp_T.columns[1]) + 2.699
W2_bd_temp_T = convert_to_numeric(data_bd_temp_T.columns[2]) + 3.399

# Filter the brown dwarfs by spectral type
M_type_mask = np.char.startswith(np.array(Sp_T_bd, dtype=str), 'M')
L_type_mask = np.char.startswith(np.array(Sp_T_bd, dtype=str), 'L')
T_type_mask = np.char.startswith(np.array(Sp_T_bd, dtype=str), 'T')
other_type_mask = ~M_type_mask & ~L_type_mask & ~T_type_mask

# Define the colors for each type
# QSO colors (distinct colors)
color_QSO_1 = "#ff6b6b"  # Darker Coral Red (from sed_plot.py QSO model spectrum, for F23)
color_QSO_2 = "#FFA500"  # Yellow/Orange (for Fl23)
color_QSO_3 = "#0000FF"  # Bright Blue (for Y23)

# BD color
color_BD = "#008000"   # Green

# Calculate number of sources in selection box
qso_araa_in_box = np.sum((W1_qso_araa - W2_qso_araa >= -0.6) & (W1_qso_araa - W2_qso_araa <= 0.6))
qso_milliquas_in_box = np.sum((W1_qso_milliquas - W2_qso_milliquas >= -0.6) & (W1_qso_milliquas - W2_qso_milliquas <= 0.6))
qso_yang_in_box = np.sum((W1_qso_yang - W2_qso_yang >= -0.6) & (W1_qso_yang - W2_qso_yang <= 0.6))
bd_in_box = np.sum((W1_bd - W2_bd >= -0.6) & (W1_bd - W2_bd <= 0.6))

# Calculate recovery percentages
qso_araa_percent = (qso_araa_in_box / len(W1_qso_araa)) * 100
qso_milliquas_percent = (qso_milliquas_in_box / len(W1_qso_milliquas)) * 100
qso_yang_percent = (qso_yang_in_box / len(W1_qso_yang)) * 100
bd_percent = (bd_in_box / len(W1_bd)) * 100

# Print the numbers
print(f"Number of sources in selection box (-0.6 to 0.6):")
print(f"Fan QSOs: {qso_araa_in_box} out of {len(W1_qso_araa)} ({qso_araa_percent:.1f}%)")
print(f"Flesch QSOs: {qso_milliquas_in_box} out of {len(W1_qso_milliquas)} ({qso_milliquas_percent:.1f}%)")
print(f"Yang QSOs: {qso_yang_in_box} out of {len(W1_qso_yang)} ({qso_yang_percent:.1f}%)")
print(f"Brown Dwarfs: {bd_in_box} out of {len(W1_bd)} ({bd_percent:.1f}%)")

# Version 1: Simple scatter plot
plt.figure(figsize=(10, 8), dpi=300)
plt.grid(True, linestyle='--', alpha=0.3)

# Sample data for plotting
W1_qso_araa_sample = sample_data(W1_qso_araa, 200)
W2_qso_araa_sample = sample_data(W2_qso_araa, 200)
W1_qso_milliquas_sample = sample_data(W1_qso_milliquas, 200)
W2_qso_milliquas_sample = sample_data(W2_qso_milliquas, 200)
W1_qso_yang_sample = sample_data(W1_qso_yang, 200)
W2_qso_yang_sample = sample_data(W2_qso_yang, 200)
# Plot BD data points (with lower opacity)
plt.scatter(W1_bd, W1_bd - W2_bd, 
           label="BD Sample", s=40, alpha=0.1, marker="o", color=color_BD, zorder=1)

# Plot QSO data points with different shades
# Y23 in the back, Fl23 in middle, F23 in front
plt.plot(W1_qso_yang_sample, W1_qso_yang_sample - W2_qso_yang_sample,'o', 
           label="Y23 Sample", alpha=0.2, marker="^", markersize=16, markerfacecolor='None', markeredgecolor=color_QSO_3, markeredgewidth=1.6, zorder=4)
plt.plot(W1_qso_milliquas_sample, W1_qso_milliquas_sample - W2_qso_milliquas_sample, 'o',
           label="Fl23 Sample", alpha=0.5, marker="D", markersize=16, markerfacecolor='None', markeredgecolor=color_QSO_2, markeredgewidth=1.6, zorder=5)
plt.plot(W1_qso_araa_sample, W1_qso_araa_sample - W2_qso_araa_sample, 'o',
           label="F23 Sample", marker="*", alpha=0.6, markersize=16, markerfacecolor='None', markeredgecolor=color_QSO_1, markeredgewidth=1.6, zorder=10)

# Add WISE color cut lines
plt.axhline(y=-0.6, color='black', linestyle='--', linewidth=1.5, label='WISE color cut')
plt.axhline(y=0.6, color='black', linestyle='--', linewidth=1.5)

# Labels and formatting
plt.xlabel("W1 (AB)", fontsize=18)
plt.ylabel("W1-W2 (AB)", fontsize=18)
plt.legend(loc="upper left", fontsize=18, ncol=2)
plt.ylim(-1, 3)
plt.xlim(14, 21)
plt.yticks([-1.0, -0.6, 0.0, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0])

# Tick parameters
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='minor', labelsize=18)

# Tight layout and save
plt.tight_layout()
plt.savefig("wise_qso_bd_scatterv1.png", dpi=300, bbox_inches='tight')
#plt.show()
plt.close()

