import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# Read data
data_qso_araa = ascii.read('delve_araa_g_band.dat')
data_qso_milliquas = ascii.read('delve_miliquas_g_band_clean.dat')
data_qso_yang = ascii.read('delve_yang_g_band.dat')
data_bd = ascii.read('delve_bd_g_band.dat')

g_mag_yang = data_qso_yang["mag_auto_g"]
g_mag_milliquas = data_qso_milliquas["mag_auto_g_1"]
g_mag_araa = data_qso_araa["mag_auto_g"]
g_mag_bd = data_bd["mag_auto_g"]

g_magerr_yang = data_qso_yang["magerr_auto_g"]
g_magerr_milliquas = data_qso_milliquas["magerr_auto_g_1"]
g_magerr_araa = data_qso_araa["magerr_auto_g"]
g_magerr_bd = data_bd["magerr_auto_g"]

snr_yang = 1.0857/g_magerr_yang
snr_milliquas = 1.0857/g_magerr_milliquas
snr_araa = 1.0857/g_magerr_araa
snr_bd = 1.0857/g_magerr_bd

plt.figure(figsize=(10, 6), dpi=200)
plt.scatter(g_mag_bd,snr_bd, label="Brown Dwarfs Sample", s=8, alpha=0.4, marker="o", color="gray")
plt.scatter(g_mag_milliquas, snr_milliquas, label="Fl23", s=16, alpha=0.4, marker="x", color="blue")
plt.scatter(g_mag_araa, snr_araa, label="F23", s=50, alpha=0.6, marker="x", color="green")
plt.scatter(g_mag_yang, snr_yang, label="Y23", s=50, alpha=0.6, marker="x", color="red")
plt.xlabel("DELVE g_mag")
plt.ylabel("Signal to Noise")
plt.legend(loc="upper left")

# Add padding to the ranges
padding = 0.5
plt.xlim(18, 30 + padding)  # Start from 18 to ignore the very low magnitude values
plt.ylim(0, 50)  # Set upper y-axis limit to 50
# Add dashed line at y=3
plt.axhline(y=3, color='black', linestyle='--', alpha=0.5)
# Set y-axis ticks to include 3
plt.yticks([3, 10, 20, 30, 40, 50])

plt.show()
plt.savefig("BD_QSO_mag_g_snr.png")
plt.close()

# Filter out data points with magnitude >= 99
g_mag_bd = g_mag_bd[g_mag_bd < 99]
g_mag_milliquas = g_mag_milliquas[g_mag_milliquas < 99]
g_mag_yang = g_mag_yang[g_mag_yang < 99]
g_mag_araa = g_mag_araa[g_mag_araa < 99]

bins = [16, 18, 20, 22, 24, 26, 30]  # Fixed bins up to 30
plt.figure(figsize=(10,6), dpi=200 )
# Create histogram
plt.hist(g_mag_bd, bins=bins, edgecolor='gray', align='left', color='white', alpha=0.5,histtype='stepfilled',  stacked='True', hatch='||', label='Brown Dwarfs Sample')
plt.hist(g_mag_milliquas, bins=bins, align='left', color='white', alpha=0.3, histtype='stepfilled', edgecolor='green', stacked='True',hatch='//', label='Fl23')
plt.hist(g_mag_yang, bins=bins, edgecolor='blue', align='left', color='white', alpha=0.3, histtype='stepfilled',  stacked='True',hatch='\\\\', label='Y23')
plt.hist(g_mag_araa, bins=bins, edgecolor='red', align='left', color='white', alpha=0.3,histtype='stepfilled',  stacked='True', hatch='----', label='F23')
# Set x-ticks to match the bin edges
plt.xticks(bins)
plt.xlabel("DELVE g_mag")
plt.legend()
# Set fixed x-axis limits
plt.xlim(16, 30)  # Fixed range from 16 to 30
plt.savefig("BD_QSO_mag_g_hist.png")
plt.show()