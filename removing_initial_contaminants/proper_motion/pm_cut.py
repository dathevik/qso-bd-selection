import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# Define colors
color_QSO_1 = "#FF0000"  # F23
color_QSO_2 = "#0000FF"   #  Fl23
color_QSO_3 = "#FFA500"     # Y23
color_BD = "#008000"  

# Read data
data_qso_milliquas = ascii.read('proper_motion_milliquas_clean.dat')
data_qso_yang = ascii.read('proper_motion_yang.dat')
data_bd = ascii.read('proper_motion_bd.dat')

# Miliquas proper motion
pm_ra_milliquas = data_qso_milliquas["pm_ra"]
pm_ra_err_milliquas = data_qso_milliquas["pm_ra_err"]
pm_dec_milliquas = data_qso_milliquas["pm_dec"]
pm_dec_err_milliquas = data_qso_milliquas["pm_dec_err"]

# Yang proper motion
pm_ra_yang = data_qso_yang["pm_ra"]
pm_ra_err_yang = data_qso_yang["pm_ra_err"]
pm_dec_yang = data_qso_yang["pm_dec"]
pm_dec_err_yang = data_qso_yang["pm_dec_err"]

# BD proper motion
pm_ra_bd = data_bd["pm_ra"]
pm_ra_err_bd = data_bd["pm_ra_err"]
pm_dec_bd = data_bd["pm_dec"]
pm_dec_err_bd = data_bd["pm_dec_err"]


# Scatter plots
plt.figure(figsize=(12, 6), dpi=300)
plt.grid(True, linestyle='--', alpha=0.3)
plt.scatter(pm_ra_bd/pm_ra_err_bd, pm_dec_bd/pm_dec_err_bd, label="BD Sample", s=60, alpha=0.4, marker="o", c=color_BD)
plt.scatter(pm_ra_milliquas/pm_ra_err_milliquas, pm_dec_milliquas/pm_dec_err_milliquas, label="Fl23 Sample", s=40, alpha=0.2, marker="^", c=color_QSO_3)
plt.scatter(pm_ra_yang/pm_ra_err_yang, pm_dec_yang/pm_dec_err_yang, label="Y23 Sample", s=160, alpha=0.3, marker="*", c=color_QSO_2)

# Labels
plt.xlabel("μₐ/σₐ", fontsize=20)
plt.ylabel("μᵧ/σᵧ", fontsize=20)

# Drawing the square for proper motion cut
square_x = [-2, 2, 2, -2, -2]
square_y = [-2, -2, 2, 2, -2]
plt.plot(square_x, square_y, color='gray', linestyle='--', linewidth=1.5, label='Proper Motion cut')

# Legend - positioned outside the plot region at the top
plt.legend(bbox_to_anchor=(0.5, 1.1), loc='center', fontsize=16, ncol=4)

# Limits
plt.ylim(-10, 10)
plt.xlim(-10, 10)
plt.yticks([-10, -5, 0, 5, 10], fontsize=20)
plt.xticks([-10, -5, 0, 5, 10], fontsize=20)

# Count objects outside the square
out_milliquas = np.sum((pm_ra_milliquas/pm_ra_err_milliquas < -2) | (pm_ra_milliquas/pm_ra_err_milliquas > 2) |
                        (pm_dec_milliquas/pm_dec_err_milliquas < -2) | (pm_dec_milliquas/pm_dec_err_milliquas > 2))

out_yang = np.sum((pm_ra_yang/pm_ra_err_yang < -2) | (pm_ra_yang/pm_ra_err_yang > 2) |
                   (pm_dec_yang/pm_dec_err_yang < -2) | (pm_dec_yang/pm_dec_err_yang > 2))

out_bd = np.sum((pm_ra_bd/pm_ra_err_bd < -2) | (pm_ra_bd/pm_ra_err_bd > 2) |
                (pm_dec_bd/pm_dec_err_bd < -2) | (pm_dec_bd/pm_dec_err_bd > 2))

nan_count = np.sum(np.isnan(pm_ra_bd) | np.isnan(pm_ra_err_bd) |
                   np.isnan(pm_dec_bd) | np.isnan(pm_dec_err_bd))

# Print results
print(f"Objects outside the square:")
print(f"  - Quasars Flesch+2023: {out_milliquas}")
print(f"  - Quasars Yang+2023: {out_yang}")
print(f"  - Quasars BD: {out_bd}")

# Show plot
#plt.show()
plt.tight_layout()
plt.savefig("QSO_BD_pm.png")
