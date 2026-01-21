import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import random

# Set random seed for reproducibility
random.seed(42)

# Read data
data_qso_araa = ascii.read('delve_class_z_araa.dat')
data_qso_milliquas = ascii.read('delve_class_z_milliquas_clean.dat')
data_qso_yang = ascii.read('delve_class_z_yang.dat')
data_bd = ascii.read('delve_class_z_bd.dat')
data_sdss = ascii.read('delve_class_z_sdss12.dat')

ext_z_yang = data_qso_yang["extended_class_z"]
ext_z_milliquas = data_qso_milliquas["extended_class_z"]
ext_z_araa = data_qso_araa["extended_class_z"]
ext_z_bd = data_bd["extended_class_z"]
ext_z_sdss = data_sdss["extended_class_z"]

# Random selection of 10% of BD sample
bd_sample_size = int(len(ext_z_bd) * 0.2)
bd_indices = random.sample(range(len(ext_z_bd)), bd_sample_size)
ext_z_bd_sample = ext_z_bd[bd_indices]

print(f"Original BD sample size: {len(ext_z_bd)}")
print(f"Sampled BD sample size (10%): {len(ext_z_bd_sample)}")

# Custom bins
bins = [-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
plt.figure(figsize=(12, 6), dpi=300)
plt.grid(True, linestyle='--', alpha=0.3)

# Define colors consistent with wise_cut_scatter.py
color_BD = "#008000"      # Green
color_Fl23 = "#0000FF"    # Bright Blue (Milliquas/Flesch)
color_F23 = "#FF0000"     # Bright Red (Fan/ARAA)
color_Y23 = "#FFA500"     # Orange (Yang)

# Create histogram with sampled BD data
plt.hist(ext_z_bd_sample, bins=bins, edgecolor=color_BD, align='left', color=color_BD, alpha=0.1, histtype='stepfilled', stacked='True', hatch='||', label=f'BD Sample (10%: {len(ext_z_bd_sample)})')
plt.hist(ext_z_milliquas, bins=bins, align='left', color=color_Fl23, alpha=0.3, histtype='stepfilled', edgecolor=color_Fl23, stacked='True', hatch='//', label='Fl23 Sample')
plt.hist(ext_z_yang, bins=bins, edgecolor=color_Y23, align='left', color=color_Y23, alpha=0.7, histtype='stepfilled', stacked='True', hatch='\\\\', label='Y23 Sample')
plt.hist(ext_z_araa, bins=bins, edgecolor=color_F23, align='left', color=color_F23, alpha=0.4, histtype='stepfilled', stacked='True', hatch='----', label='F23 Sample')
plt.hist(ext_z_sdss, bins=bins, edgecolor='purple', align='left', color='lavender', alpha=0.7, histtype='stepfilled', stacked='True', hatch='++', label='Galaxies')

plt.xticks([-1, 0, 1, 2, 3], fontsize=20)
plt.xlim(-1, 4)
plt.ylim(0, 1200)

plt.yticks([0, 250, 500, 750, 1000], fontsize=20)

plt.xlabel('extended_class_z_delve', fontsize=20)
plt.ylabel('Number', fontsize=20)
plt.legend(loc="upper right", fontsize=20, framealpha=0.8)

plt.tight_layout()
plt.savefig("QSO_BD_extended_z.png")