#!/usr/bin/env python
"""
2D and 3D Separation Plot Script

This script reads two .dat files containing astronomical coordinates,
calculates the angular separations between ra/dec and ra_wise/dec_wise,
and creates:
- 2D separation plots (RA offset vs Dec offset)
- 3D separation plots (RA offset, Dec offset, Angular separation)
- Separation distribution histograms

Based on astropy.coordinates documentation:
https://docs.astropy.org/en/stable/coordinates/matchsep.html#astropy-coordinates-separations-matching
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u

# Read the two .dat files
print("Reading wise_F23.dat...")
f23_data = ascii.read('wise_F23.dat')

print("Reading wise_TM.dat...")
tm_data = ascii.read('wise_TM.dat')

# Create SkyCoord objects for the original coordinates (ra, dec)
print("Creating SkyCoord objects for original coordinates...")
f23_coords = SkyCoord(ra=f23_data['ra']*u.degree, dec=f23_data['dec']*u.degree, frame='icrs')
tm_coords = SkyCoord(ra=tm_data['ra']*u.degree, dec=tm_data['dec']*u.degree, frame='icrs')

# Create SkyCoord objects for the WISE coordinates (ra_wise, dec_wise)
print("Creating SkyCoord objects for WISE coordinates...")
f23_wise_coords = SkyCoord(ra=f23_data['ra_wise']*u.degree, dec=f23_data['dec_wise']*u.degree, frame='icrs')
tm_wise_coords = SkyCoord(ra=tm_data['ra_wise']*u.degree, dec=tm_data['dec_wise']*u.degree, frame='icrs')

# Calculate separations using the separation() method (following astropy documentation)
print("Calculating separations for F23 data...")
f23_separations = f23_coords.separation(f23_wise_coords)

print("Calculating separations for TM data...")
tm_separations = tm_coords.separation(tm_wise_coords)

# Calculate RA and Dec differences (in arcseconds)
# Using spherical_offsets_to() method to get RA and Dec offsets
print("Calculating RA and Dec offsets...")
f23_dra, f23_ddec = f23_coords.spherical_offsets_to(f23_wise_coords)
tm_dra, tm_ddec = tm_coords.spherical_offsets_to(tm_wise_coords)

# Convert to arcseconds for plotting (extract values from Angle objects)
f23_dra_arcsec = f23_dra.to(u.arcsec).value
f23_ddec_arcsec = f23_ddec.to(u.arcsec).value
tm_dra_arcsec = tm_dra.to(u.arcsec).value
tm_ddec_arcsec = tm_ddec.to(u.arcsec).value

# Convert separations to arcseconds for plotting (extract values from Angle objects)
f23_sep_arcsec = f23_separations.to(u.arcsec).value
tm_sep_arcsec = tm_separations.to(u.arcsec).value

# Create the 2D separation plots
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: RA vs Dec separations for F23
ax1 = axes[0]
ax1.scatter(f23_dra_arcsec, f23_ddec_arcsec, alpha=0.6, s=30, label=f'F23 (N={len(f23_data)})', color='#ff6b6b')
ax1.set_xlabel('RA Offset (arcsec)', fontsize=14)
ax1.set_ylabel('Dec Offset (arcsec)', fontsize=14)
ax1.set_title('F23: WISE vs Original Coordinates', fontsize=16)
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=12)
ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
ax1.axvline(x=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)

# Plot 2: RA vs Dec separations for TM
ax2 = axes[1]
ax2.scatter(tm_dra_arcsec, tm_ddec_arcsec, alpha=0.6, s=30, label=f'TM (N={len(tm_data)})', color='hotpink')
ax2.set_xlabel('RA Offset (arcsec)', fontsize=14)
ax2.set_ylabel('Dec Offset (arcsec)', fontsize=14)
ax2.set_title('TM: WISE vs Original Coordinates', fontsize=16)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=12)
ax2.axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
ax2.axvline(x=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)

# Make axes symmetric for better visualization
for ax in axes:
    # Get current limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    # Make symmetric around zero
    max_range = max(abs(xlim[0]), abs(xlim[1]), abs(ylim[0]), abs(ylim[1]))
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)

plt.tight_layout()
plt.savefig('2d_separation_plot.png', dpi=300, bbox_inches='tight')
print("Saved: 2d_separation_plot.png")
plt.close()

# ========== 3D Separation Plots ==========
# Create 3D scatter plots showing RA offset, Dec offset, and angular separation
fig = plt.figure(figsize=(16, 6))

# 3D Plot 1: F23 data
ax1_3d = fig.add_subplot(121, projection='3d')
scatter1 = ax1_3d.scatter(f23_dra_arcsec, f23_ddec_arcsec, f23_sep_arcsec, 
                         c=f23_sep_arcsec, cmap='viridis', alpha=0.6, s=30)
ax1_3d.set_xlabel('RA Offset (arcsec)', fontsize=12)
ax1_3d.set_ylabel('Dec Offset (arcsec)', fontsize=12)
ax1_3d.set_zlabel('Angular Separation (arcsec)', fontsize=12)
ax1_3d.set_title(f'F23: 3D Separation (N={len(f23_data)})', fontsize=14)
plt.colorbar(scatter1, ax=ax1_3d, label='Separation (arcsec)', shrink=0.8)

# 3D Plot 2: TM data
ax2_3d = fig.add_subplot(122, projection='3d')
scatter2 = ax2_3d.scatter(tm_dra_arcsec, tm_ddec_arcsec, tm_sep_arcsec,
                          c=tm_sep_arcsec, cmap='plasma', alpha=0.6, s=30)
ax2_3d.set_xlabel('RA Offset (arcsec)', fontsize=12)
ax2_3d.set_ylabel('Dec Offset (arcsec)', fontsize=12)
ax2_3d.set_zlabel('Angular Separation (arcsec)', fontsize=12)
ax2_3d.set_title(f'TM: 3D Separation (N={len(tm_data)})', fontsize=14)
plt.colorbar(scatter2, ax=ax2_3d, label='Separation (arcsec)', shrink=0.8)

plt.tight_layout()
plt.savefig('3d_separation_plot.png', dpi=300, bbox_inches='tight')
print("Saved: 3d_separation_plot.png")
plt.close()

# ========== Separation Histograms ==========
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Histogram 1: F23 separations
ax1_hist = axes[0]
ax1_hist.hist(f23_sep_arcsec, bins=30, alpha=0.7, color='#ff6b6b', edgecolor='black', linewidth=1.2)
ax1_hist.axvline(np.mean(f23_sep_arcsec), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(f23_sep_arcsec):.2f} arcsec')
ax1_hist.axvline(np.median(f23_sep_arcsec), color='blue', linestyle='--', linewidth=2, label=f'Median: {np.median(f23_sep_arcsec):.2f} arcsec')
ax1_hist.set_xlabel('Angular Separation (arcsec)', fontsize=14)
ax1_hist.set_ylabel('Count', fontsize=14)
ax1_hist.set_title(f'F23: Separation Distribution (N={len(f23_data)})', fontsize=16)
ax1_hist.legend(fontsize=12)
ax1_hist.grid(True, alpha=0.3)

# Histogram 2: TM separations
ax2_hist = axes[1]
ax2_hist.hist(tm_sep_arcsec, bins=30, alpha=0.7, color='hotpink', edgecolor='black', linewidth=1.2)
ax2_hist.axvline(np.mean(tm_sep_arcsec), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(tm_sep_arcsec):.2f} arcsec')
ax2_hist.axvline(np.median(tm_sep_arcsec), color='blue', linestyle='--', linewidth=2, label=f'Median: {np.median(tm_sep_arcsec):.2f} arcsec')
ax2_hist.set_xlabel('Angular Separation (arcsec)', fontsize=14)
ax2_hist.set_ylabel('Count', fontsize=14)
ax2_hist.set_title(f'TM: Separation Distribution (N={len(tm_data)})', fontsize=16)
ax2_hist.legend(fontsize=12)
ax2_hist.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('separation_histograms.png', dpi=300, bbox_inches='tight')
print("Saved: separation_histograms.png")
plt.close()

# ========== RA/Dec Distribution with Offset Vectors ==========
# Plot showing positions with arrows indicating the offset from original to WISE coordinates
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Plot 1: F23 data
ax1_vec = axes[0]

# Plot WISE positions as points
ax1_vec.scatter(f23_data['ra_wise'], f23_data['dec_wise'], 
                alpha=0.4, s=20, color='lightblue', label='WISE positions', zorder=1)

# Plot original positions as points
ax1_vec.scatter(f23_data['ra'], f23_data['dec'], 
                alpha=0.4, s=20, color='red', marker='x', label='Original positions', zorder=2)

# Draw arrows from original to WISE positions
# Sample the data if there are too many points for clarity
if len(f23_data) > 200:
    sample_indices = np.random.choice(len(f23_data), 200, replace=False)
    sample_indices = np.sort(sample_indices)
else:
    sample_indices = np.arange(len(f23_data))

for i in sample_indices:
    ax1_vec.annotate('', 
                     xy=(f23_data['ra_wise'][i], f23_data['dec_wise'][i]),  # Arrow tip at WISE position
                     xytext=(f23_data['ra'][i], f23_data['dec'][i]),  # Arrow start at original position
                     arrowprops=dict(arrowstyle='->', color='black', alpha=0.3, lw=0.8, 
                                   connectionstyle='arc3,rad=0'))

ax1_vec.set_xlabel('Right Ascension (degrees)', fontsize=14)
ax1_vec.set_ylabel('Declination (degrees)', fontsize=14)
ax1_vec.set_title(f'F23: Position Distribution with Offsets (N={len(f23_data)})', fontsize=16)
ax1_vec.legend(fontsize=12, loc='best')
ax1_vec.grid(True, alpha=0.3)
ax1_vec.invert_xaxis()  # Astronomical convention: RA increases to the left

# Plot 2: TM data
ax2_vec = axes[1]

# Plot WISE positions as points
ax2_vec.scatter(tm_data['ra_wise'], tm_data['dec_wise'], 
                alpha=0.4, s=20, color='lightblue', label='WISE positions', zorder=1)

# Plot original positions as points
ax2_vec.scatter(tm_data['ra'], tm_data['dec'], 
                alpha=0.4, s=20, color='hotpink', marker='x', label='Original positions', zorder=2)

# Draw arrows from original to WISE positions
# Sample the data if there are too many points for clarity
if len(tm_data) > 200:
    sample_indices = np.random.choice(len(tm_data), 200, replace=False)
    sample_indices = np.sort(sample_indices)
else:
    sample_indices = np.arange(len(tm_data))

for i in sample_indices:
    ax2_vec.annotate('', 
                     xy=(tm_data['ra_wise'][i], tm_data['dec_wise'][i]),  # Arrow tip at WISE position
                     xytext=(tm_data['ra'][i], tm_data['dec'][i]),  # Arrow start at original position
                     arrowprops=dict(arrowstyle='->', color='black', alpha=0.3, lw=0.8,
                                   connectionstyle='arc3,rad=0'))

ax2_vec.set_xlabel('Right Ascension (degrees)', fontsize=14)
ax2_vec.set_ylabel('Declination (degrees)', fontsize=14)
ax2_vec.set_title(f'TM: Position Distribution with Offsets (N={len(tm_data)})', fontsize=16)
ax2_vec.legend(fontsize=12, loc='best')
ax2_vec.grid(True, alpha=0.3)
ax2_vec.invert_xaxis()  # Astronomical convention: RA increases to the left

plt.tight_layout()
plt.savefig('ra_dec_distribution_with_offsets.png', dpi=300, bbox_inches='tight')
print("Saved: ra_dec_distribution_with_offsets.png")
plt.close()

# Alternative: Quiver plot for better visualization of vector field
fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Plot 1: F23 quiver plot
ax1_q = axes[0]

# Sample data for quiver plot (too many arrows can be cluttered)
if len(f23_data) > 300:
    sample_indices = np.random.choice(len(f23_data), 300, replace=False)
    sample_indices = np.sort(sample_indices)
else:
    sample_indices = np.arange(len(f23_data))

# Calculate offsets in degrees (for quiver)
f23_dra_deg = f23_dra.to(u.degree).value
f23_ddec_deg = f23_ddec.to(u.degree).value

# Plot original positions
ax1_q.scatter(f23_data['ra'][sample_indices], f23_data['dec'][sample_indices],
              alpha=0.5, s=30, color='red', marker='o', label='Original positions', zorder=2)

# Draw quiver arrows showing offset direction and magnitude
# Scale arrows to be visible (normalize by max separation)
max_arrow_scale = np.max(f23_sep_arcsec) / 10.0  # Scale factor for visibility
q1 = ax1_q.quiver(f23_data['ra'][sample_indices], f23_data['dec'][sample_indices],
                  f23_dra_deg[sample_indices], f23_ddec_deg[sample_indices],
                  f23_sep_arcsec[sample_indices],  # Color by separation
                  angles='xy', scale_units='xy', scale=1.0/max_arrow_scale,
                  alpha=0.6, width=0.003, cmap='viridis', zorder=1, clim=(0, np.max(f23_sep_arcsec)))
plt.colorbar(q1, ax=ax1_q, label='Separation (arcsec)', shrink=0.8)

ax1_q.set_xlabel('Right Ascension (degrees)', fontsize=14)
ax1_q.set_ylabel('Declination (degrees)', fontsize=14)
ax1_q.set_title(f'F23: Offset Vector Field (N={len(sample_indices)} shown)', fontsize=16)
ax1_q.legend(fontsize=12, loc='best')
ax1_q.grid(True, alpha=0.3)
ax1_q.invert_xaxis()

# Plot 2: TM quiver plot
ax2_q = axes[1]

# Sample data for quiver plot
if len(tm_data) > 300:
    sample_indices = np.random.choice(len(tm_data), 300, replace=False)
    sample_indices = np.sort(sample_indices)
else:
    sample_indices = np.arange(len(tm_data))

# Calculate offsets in degrees
tm_dra_deg = tm_dra.to(u.degree).value
tm_ddec_deg = tm_ddec.to(u.degree).value

# Plot original positions
ax2_q.scatter(tm_data['ra'][sample_indices], tm_data['dec'][sample_indices],
              alpha=0.5, s=30, color='hotpink', marker='o', label='Original positions', zorder=2)

# Draw quiver arrows
max_arrow_scale_tm = np.max(tm_sep_arcsec) / 10.0
q2 = ax2_q.quiver(tm_data['ra'][sample_indices], tm_data['dec'][sample_indices],
                  tm_dra_deg[sample_indices], tm_ddec_deg[sample_indices],
                  tm_sep_arcsec[sample_indices],  # Color by separation
                  angles='xy', scale_units='xy', scale=1.0/max_arrow_scale_tm,
                  alpha=0.6, width=0.003, cmap='plasma', zorder=1, clim=(0, np.max(tm_sep_arcsec)))
plt.colorbar(q2, ax=ax2_q, label='Separation (arcsec)', shrink=0.8)

ax2_q.set_xlabel('Right Ascension (degrees)', fontsize=14)
ax2_q.set_ylabel('Declination (degrees)', fontsize=14)
ax2_q.set_title(f'TM: Offset Vector Field (N={len(sample_indices)} shown)', fontsize=16)
ax2_q.legend(fontsize=12, loc='best')
ax2_q.grid(True, alpha=0.3)
ax2_q.invert_xaxis()

plt.tight_layout()
plt.savefig('ra_dec_vector_field.png', dpi=300, bbox_inches='tight')
print("Saved: ra_dec_vector_field.png")
plt.close()

# Print statistics
print("\n=== F23 Statistics ===")
print(f"Total sources: {len(f23_data)}")
print(f"Mean separation: {np.mean(f23_sep_arcsec):.3f} arcsec")
print(f"Median separation: {np.median(f23_sep_arcsec):.3f} arcsec")
print(f"Max separation: {np.max(f23_sep_arcsec):.3f} arcsec")
print(f"Mean RA offset: {np.mean(f23_dra_arcsec):.3f} arcsec")
print(f"Mean Dec offset: {np.mean(f23_ddec_arcsec):.3f} arcsec")

print("\n=== TM Statistics ===")
print(f"Total sources: {len(tm_data)}")
print(f"Mean separation: {np.mean(tm_sep_arcsec):.3f} arcsec")
print(f"Median separation: {np.median(tm_sep_arcsec):.3f} arcsec")
print(f"Max separation: {np.max(tm_sep_arcsec):.3f} arcsec")
print(f"Mean RA offset: {np.mean(tm_dra_arcsec):.3f} arcsec")
print(f"Mean Dec offset: {np.mean(tm_ddec_arcsec):.3f} arcsec")

plt.show()

