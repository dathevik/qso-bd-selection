import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

def sdss_asinh_to_pogson(mag_asinh, band):
    """
    Convert SDSS asinh magnitudes to Pogson magnitudes.
    
    Parameters:
    -----------
    mag_asinh : float or array
        SDSS asinh magnitude(s)
    band : str
        Photometric band: 'u', 'g', 'r', 'i', or 'z'
    
    Returns:
    --------
    mag_pogson : float or array
        Pogson magnitude(s)
    """
    
    # SDSS asinh softening parameters (in maggies)
    softening_params = {
        'u': {'b': 1.4e-10, 'm_zero_flux': 24.63},
        'g': {'b': 0.9e-10, 'm_zero_flux': 25.11},
        'r': {'b': 1.2e-10, 'm_zero_flux': 24.80},
        'i': {'b': 1.8e-10, 'm_zero_flux': 24.36},
        'z': {'b': 7.4e-10, 'm_zero_flux': 22.83}
    }
    
    if band not in softening_params:
        raise ValueError(f"Band must be one of {list(softening_params.keys())}")
    
    b = softening_params[band]['b']
    m_zero_flux = softening_params[band]['m_zero_flux']
    
    # Step 1: Convert asinh magnitude to flux (in units of f/f₀, i.e., maggies)
    # SDSS formula: m_asinh = m₀ - (2.5/ln(10)) × [asinh((f/f₀)/(2b)) + ln(b)]
    # Inverting to get f/f₀:
    # (m₀ - m_asinh) × (ln(10)/2.5) = asinh((f/f₀)/(2b)) + ln(b)
    # (m₀ - m_asinh) × (ln(10)/2.5) - ln(b) = asinh((f/f₀)/(2b))
    # f/f₀ = 2b × sinh[(m₀ - m_asinh) × (ln(10)/2.5) - ln(b)]
    sinh_arg = (np.log(10) / 2.5) * (m_zero_flux - mag_asinh) - np.log(b)
    f_over_f0 = 2 * b * np.sinh(sinh_arg)
    
    # Step 2: Convert flux to Pogson magnitude
    # m_Pogson = m₀ - 2.5 × log₁₀(f/f₀)
    # For SDSS, m₀ is typically 22.5, but we use the band-specific zero-point
    mag_pogson = m_zero_flux - 2.5 * np.log10(f_over_f0)
    
    return mag_pogson

def convert_sdss_to_delve(g_sdss, r_sdss, i_sdss, z_sdss):
    """
    Convert SDSS magnitudes to DELVE (DES) system.
    Equations from DES DR2 documentation.
    """
    # Calculate color terms
    g_i_sdss = g_sdss - i_sdss
    r_i_sdss = r_sdss - i_sdss
    
    # Apply transformations (DES = DELVE)
    g_delve = g_sdss - 0.061 * g_i_sdss + 0.008
    r_delve = r_sdss - 0.155 * r_i_sdss - 0.007
    i_delve = i_sdss - 0.166 * r_i_sdss + 0.032
    z_delve = z_sdss - 0.056 * r_i_sdss + 0.027
    
    return g_delve, r_delve, i_delve, z_delve

def convert_ps1_to_delve(g_ps1, r_ps1, i_ps1, z_ps1):
    """
    Convert Pan-STARRS1 magnitudes to DELVE (DES) system.
    Equations from DES DR2 documentation.
    """
    # Calculate color terms
    g_i_ps1 = g_ps1 - i_ps1
    r_i_ps1 = r_ps1 - i_ps1
    
    # Apply transformations (DES = DELVE)
    g_delve = g_ps1 + 0.028 * g_i_ps1 + 0.020
    r_delve = r_ps1 - 0.142 * r_i_ps1 - 0.010
    i_delve = i_ps1 - 0.155 * r_i_ps1 + 0.015
    z_delve = z_ps1 - 0.114 * r_i_ps1 - 0.010
    
    return g_delve, r_delve, i_delve, z_delve

# Read data from file
data = ascii.read('FL_23_all_photo.dat')

# Get magnitudes from DECALS
r_mag_decals = data['mag_r_decals']
i_mag_decals = data['mag_i']
z_mag_decals = data['mag_z_decals']

# Get magnitudes from DELVE
r_mag_delve = data['mag_auto_r_delve']
i_mag_delve = data['mag_auto_i_delve']
z_mag_delve = data['mag_auto_z_delve']

# Get magnitudes from SDSS (cmodelmag) - asinh magnitudes
print("\n=== SDSS Data Diagnostics ===")
print(f"Total rows in data: {len(data)}")
print(f"SDSS columns available: {[col for col in data.colnames if 'cmodelmag' in col.lower()]}")

try:
    r_mag_sdss_asinh = data['cmodelmag_r']
    i_mag_sdss_asinh = data['cmodelmag_i']
    z_mag_sdss_asinh = data['cmodelmag_z']
    
    print(f"SDSS r-band: {np.sum(np.isfinite(r_mag_sdss_asinh))} valid values out of {len(r_mag_sdss_asinh)}")
    print(f"SDSS i-band: {np.sum(np.isfinite(i_mag_sdss_asinh))} valid values out of {len(i_mag_sdss_asinh)}")
    print(f"SDSS z-band: {np.sum(np.isfinite(z_mag_sdss_asinh))} valid values out of {len(z_mag_sdss_asinh)}")
    print(f"Sample SDSS r asinh values: {r_mag_sdss_asinh[:5]}")
    print(f"SDSS r range: {np.nanmin(r_mag_sdss_asinh):.2f} to {np.nanmax(r_mag_sdss_asinh):.2f}")
    
    # Convert SDSS asinh magnitudes to Pogson magnitudes
    print("\nConverting SDSS asinh to Pogson...")
    r_mag_sdss = sdss_asinh_to_pogson(r_mag_sdss_asinh, 'r')
    i_mag_sdss = sdss_asinh_to_pogson(i_mag_sdss_asinh, 'i')
    z_mag_sdss = sdss_asinh_to_pogson(z_mag_sdss_asinh, 'z')
    print(f"After conversion - r: {np.sum(np.isfinite(r_mag_sdss))} valid, i: {np.sum(np.isfinite(i_mag_sdss))} valid, z: {np.sum(np.isfinite(z_mag_sdss))} valid")
    print(f"Sample converted r Pogson values: {r_mag_sdss[:5]}")
    print(f"Converted r range: {np.nanmin(r_mag_sdss):.2f} to {np.nanmax(r_mag_sdss):.2f}")
    
except KeyError as e:
    print(f"ERROR: SDSS column not found: {e}")
    print(f"Available columns containing 'mag': {[col for col in data.colnames if 'mag' in col.lower()][:20]}")
    # Create dummy arrays to prevent errors
    r_mag_sdss = np.full(len(data), np.nan)
    i_mag_sdss = np.full(len(data), np.nan)
    z_mag_sdss = np.full(len(data), np.nan)

# Get magnitudes from PS1 (MeanPSFMag)
r_mag_ps1 = data['rMeanPSFMag']
i_mag_ps1 = data['iMeanPSFMag']
z_mag_ps1 = data['zMeanPSFMag']

# ============================================================================
# Convert SDSS and PS1 magnitudes to DELVE system
# ============================================================================

print("\n=== Converting SDSS to DELVE system ===")
# Get g-band SDSS if available for color calculation
if 'cmodelmag_g' in data.colnames:
    g_mag_sdss_asinh = data['cmodelmag_g']
    g_mag_sdss = sdss_asinh_to_pogson(g_mag_sdss_asinh, 'g')
    print("g-band SDSS available, converting SDSS to DELVE...")
    _, r_mag_sdss_delve, i_mag_sdss_delve, z_mag_sdss_delve = convert_sdss_to_delve(
        g_mag_sdss, r_mag_sdss, i_mag_sdss, z_mag_sdss)
    print(f"SDSS converted to DELVE system")
else:
    print("g-band SDSS not available, using r,i,z only (approximate conversion)")
    # For r,i,z we can still use r-i color
    r_i_sdss = r_mag_sdss - i_mag_sdss
    r_mag_sdss_delve = r_mag_sdss - 0.155 * r_i_sdss - 0.007
    i_mag_sdss_delve = i_mag_sdss - 0.166 * r_i_sdss + 0.032
    z_mag_sdss_delve = z_mag_sdss - 0.056 * r_i_sdss + 0.027

print("\n=== Converting PS1 to DELVE system ===")
# Get g-band PS1 if available
if 'gMeanPSFMag' in data.colnames:
    g_mag_ps1 = data['gMeanPSFMag']
    print("g-band PS1 available, converting PS1 to DELVE...")
    _, r_mag_ps1_delve, i_mag_ps1_delve, z_mag_ps1_delve = convert_ps1_to_delve(
        g_mag_ps1, r_mag_ps1, i_mag_ps1, z_mag_ps1)
    print(f"PS1 converted to DELVE system")
else:
    print("g-band PS1 not available, using r,i,z only (approximate conversion)")
    # For r,i,z we can still use r-i color
    r_i_ps1 = r_mag_ps1 - i_mag_ps1
    r_mag_ps1_delve = r_mag_ps1 - 0.142 * r_i_ps1 - 0.010
    i_mag_ps1_delve = i_mag_ps1 - 0.155 * r_i_ps1 + 0.015
    z_mag_ps1_delve = z_mag_ps1 - 0.114 * r_i_ps1 - 0.010

# ============================================================================
# Calculate magnitude differences (now all in DELVE system)
# ============================================================================

# DELVE - DECALS differences
diff_r_delve_decals = r_mag_delve - r_mag_decals
diff_i_delve_decals = i_mag_delve - i_mag_decals
diff_z_delve_decals = z_mag_delve - z_mag_decals

# DELVE - SDSS differences (SDSS converted to DELVE system)
diff_r_delve_sdss = r_mag_delve - r_mag_sdss_delve
diff_i_delve_sdss = i_mag_delve - i_mag_sdss_delve
diff_z_delve_sdss = z_mag_delve - z_mag_sdss_delve

# DELVE - PS1 differences (PS1 converted to DELVE system)
diff_r_delve_ps1 = r_mag_delve - r_mag_ps1_delve
diff_i_delve_ps1 = i_mag_delve - i_mag_ps1_delve
diff_z_delve_ps1 = z_mag_delve - z_mag_ps1_delve

# ============================================================================
# Filter out invalid values (NaN, inf, and very large values)
# ============================================================================

def filter_valid(diff_array):
    """Filter out NaN, inf, and extreme values."""
    valid = np.isfinite(diff_array) & (np.abs(diff_array) < 10.0)
    return diff_array[valid]

# Filter all differences
diff_r_delve_decals_clean = filter_valid(diff_r_delve_decals)
diff_i_delve_decals_clean = filter_valid(diff_i_delve_decals)
diff_z_delve_decals_clean = filter_valid(diff_z_delve_decals)

diff_r_delve_sdss_clean = filter_valid(diff_r_delve_sdss)
diff_i_delve_sdss_clean = filter_valid(diff_i_delve_sdss)
diff_z_delve_sdss_clean = filter_valid(diff_z_delve_sdss)

print(f"\nAfter filtering DELVE-SDSS differences:")
print(f"  r-band: {len(diff_r_delve_sdss_clean)} valid values")
print(f"  i-band: {len(diff_i_delve_sdss_clean)} valid values")
print(f"  z-band: {len(diff_z_delve_sdss_clean)} valid values")

diff_r_delve_ps1_clean = filter_valid(diff_r_delve_ps1)
diff_i_delve_ps1_clean = filter_valid(diff_i_delve_ps1)
diff_z_delve_ps1_clean = filter_valid(diff_z_delve_ps1)

# ============================================================================
# Plot 1: DELVE - DECALS differences
# ============================================================================

fig1, axes1 = plt.subplots(3, 1, figsize=(10, 12))
fig1.suptitle('DELVE - DECALS Magnitude Differences', fontsize=16, fontweight='bold')

# r-band
axes1[0].hist(diff_r_delve_decals_clean, bins=50, alpha=0.7, color='red', edgecolor='black')
axes1[0].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes1[0].set_xlabel('r-band: DELVE DR2 - DECaLS DR10', fontsize=12)
axes1[0].set_ylabel('Number of sources', fontsize=12)
axes1[0].set_xlim(-0.6, 0.6)
axes1[0].grid(True, alpha=0.3)

# i-band
axes1[1].hist(diff_i_delve_decals_clean, bins=50, alpha=0.7, color='green', edgecolor='black')
axes1[1].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes1[1].set_xlabel('i-band: DELVE DR2 - DECaLS DR10', fontsize=12)
axes1[1].set_ylabel('Number of sources', fontsize=12)
axes1[1].set_xlim(-0.6, 0.6)
axes1[1].grid(True, alpha=0.3)

# z-band
axes1[2].hist(diff_z_delve_decals_clean, bins=50, alpha=0.7, color='blue', edgecolor='black')
axes1[2].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes1[2].set_xlabel('z-band: DELVE DR2 - DECaLS DR10', fontsize=12)
axes1[2].set_ylabel('Number of sources', fontsize=12)
axes1[2].set_xlim(-0.6, 0.6)
axes1[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('delve_decals_differences.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# Plot 2: DELVE - SDSS differences (all Pogson)
# ============================================================================

fig2, axes2 = plt.subplots(3, 1, figsize=(10, 12))
fig2.suptitle('DELVE - SDSS Magnitude Differences', fontsize=16, fontweight='bold')

# r-band
axes2[0].hist(diff_r_delve_sdss_clean, bins=50, alpha=0.7, color='red', edgecolor='black')
axes2[0].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes2[0].set_xlabel('r-band: DELVE DR2 - SDSS DR17', fontsize=12)
axes2[0].set_ylabel('Number of sources', fontsize=12)
axes2[0].set_xlim(-0.6, 0.6)
axes2[0].grid(True, alpha=0.3)

# i-band
axes2[1].hist(diff_i_delve_sdss_clean, bins=50, alpha=0.7, color='green', edgecolor='black')
axes2[1].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes2[1].set_xlabel('i-band: DELVE DR2 - SDSS DR17', fontsize=12)
axes2[1].set_ylabel('Number of sources', fontsize=12)
axes2[1].set_xlim(-0.6, 0.6)
axes2[1].grid(True, alpha=0.3)

# z-band
axes2[2].hist(diff_z_delve_sdss_clean, bins=50, alpha=0.7, color='blue', edgecolor='black')
axes2[2].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes2[2].set_xlabel('z-band: DELVE DR2 - SDSS DR17', fontsize=12)
axes2[2].set_ylabel('Number of sources', fontsize=12)
axes2[2].set_xlim(-0.6, 0.6)
axes2[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('delve_sdss_differences.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# Plot 3: DELVE - PS1 differences
# ============================================================================

fig3, axes3 = plt.subplots(3, 1, figsize=(10, 12))
fig3.suptitle('DELVE - PS1 Magnitude Differences', fontsize=16, fontweight='bold')

# r-band
axes3[0].hist(diff_r_delve_ps1_clean, bins=50, alpha=0.7, color='red', edgecolor='black')
axes3[0].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes3[0].set_xlabel('r-band: DELVE DR2 - PS1', fontsize=12)
axes3[0].set_ylabel('Number of sources', fontsize=12)
axes3[0].set_xlim(-0.6, 0.6)
axes3[0].grid(True, alpha=0.3)

# i-band
axes3[1].hist(diff_i_delve_ps1_clean, bins=50, alpha=0.7, color='green', edgecolor='black')
axes3[1].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes3[1].set_xlabel('i-band: DELVE DR2 - PS1', fontsize=12)
axes3[1].set_ylabel('Number of sources', fontsize=12)
axes3[1].set_xlim(-0.6, 0.6)
axes3[1].grid(True, alpha=0.3)

# z-band
axes3[2].hist(diff_z_delve_ps1_clean, bins=50, alpha=0.7, color='blue', edgecolor='black')
axes3[2].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes3[2].set_xlabel('z-band: DELVE DR2 - PS1', fontsize=12)
axes3[2].set_ylabel('Number of sources', fontsize=12)
axes3[2].set_xlim(-0.6, 0.6)
axes3[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('delve_ps1_differences.png', dpi=300, bbox_inches='tight')
plt.close()

print("All plots saved successfully!")
print(f"  - delve_decals_differences.png")
print(f"  - delve_sdss_differences.png")
print(f"  - delve_ps1_differences.png")

# ============================================================================
# TM Photo Data Plots
# ============================================================================

print("\n=== Processing final_TM_photo.dat ===")

# Read data from TM file
data_tm = ascii.read('final_TM_photo.dat')

print(f"Total rows in TM data: {len(data_tm)}")

# Get magnitudes from DECALS
r_mag_decals_tm = data_tm['mag_r_decals']
i_mag_decals_tm = data_tm['mag_i_decals']
z_mag_decals_tm = data_tm['mag_z_decals']

# Get magnitudes from DELVE
r_mag_delve_tm = data_tm['mag_r_delve']
i_mag_delve_tm = data_tm['mag_i_delve']
z_mag_delve_tm = data_tm['mag_z_delve']

# Get magnitudes from PS1 (MeanPSFMag)
r_mag_ps1_tm = data_tm['rMeanPSFMag']
i_mag_ps1_tm = data_tm['iMeanPSFMag']
z_mag_ps1_tm = data_tm['zMeanPSFMag']

# Convert PS1 to DELVE system for TM data
print("\n=== Converting PS1 to DELVE system (TM) ===")
if 'gMeanPSFMag' in data_tm.colnames:
    g_mag_ps1_tm = data_tm['gMeanPSFMag']
    print("g-band PS1 available, converting PS1 to DELVE...")
    _, r_mag_ps1_delve_tm, i_mag_ps1_delve_tm, z_mag_ps1_delve_tm = convert_ps1_to_delve(
        g_mag_ps1_tm, r_mag_ps1_tm, i_mag_ps1_tm, z_mag_ps1_tm)
    print(f"PS1 converted to DELVE system")
else:
    print("g-band PS1 not available, using r,i,z only (approximate conversion)")
    r_i_ps1_tm = r_mag_ps1_tm - i_mag_ps1_tm
    r_mag_ps1_delve_tm = r_mag_ps1_tm - 0.142 * r_i_ps1_tm - 0.010
    i_mag_ps1_delve_tm = i_mag_ps1_tm - 0.155 * r_i_ps1_tm + 0.015
    z_mag_ps1_delve_tm = z_mag_ps1_tm - 0.114 * r_i_ps1_tm - 0.010

# Calculate magnitude differences
diff_r_delve_decals_tm = r_mag_delve_tm - r_mag_decals_tm
diff_i_delve_decals_tm = i_mag_delve_tm - i_mag_decals_tm
diff_z_delve_decals_tm = z_mag_delve_tm - z_mag_decals_tm

# Filter valid values for DECALS
diff_r_delve_decals_tm_clean = filter_valid(diff_r_delve_decals_tm)
diff_i_delve_decals_tm_clean = filter_valid(diff_i_delve_decals_tm)
diff_z_delve_decals_tm_clean = filter_valid(diff_z_delve_decals_tm)

# Calculate and filter PS1 differences (PS1 converted to DELVE system)
diff_r_delve_ps1_tm = r_mag_delve_tm - r_mag_ps1_delve_tm
diff_i_delve_ps1_tm = i_mag_delve_tm - i_mag_ps1_delve_tm
diff_z_delve_ps1_tm = z_mag_delve_tm - z_mag_ps1_delve_tm

diff_r_delve_ps1_tm_clean = filter_valid(diff_r_delve_ps1_tm)
diff_i_delve_ps1_tm_clean = filter_valid(diff_i_delve_ps1_tm)
diff_z_delve_ps1_tm_clean = filter_valid(diff_z_delve_ps1_tm)

print(f"After filtering DELVE-DECALS TM differences:")
print(f"  r-band: {len(diff_r_delve_decals_tm_clean)} valid values")
print(f"  i-band: {len(diff_i_delve_decals_tm_clean)} valid values")
print(f"  z-band: {len(diff_z_delve_decals_tm_clean)} valid values")

print(f"After filtering DELVE-PS1 TM differences:")
print(f"  r-band: {len(diff_r_delve_ps1_tm_clean)} valid values")
print(f"  i-band: {len(diff_i_delve_ps1_tm_clean)} valid values")
print(f"  z-band: {len(diff_z_delve_ps1_tm_clean)} valid values")

# ============================================================================
# Plot 1: DELVE - DECALS differences (TM)
# ============================================================================

fig1_tm, axes1_tm = plt.subplots(3, 1, figsize=(10, 12))
fig1_tm.suptitle('DELVE - DECALS Magnitude Differences (TM)', fontsize=16, fontweight='bold')

# r-band
axes1_tm[0].hist(diff_r_delve_decals_tm_clean, bins=50, alpha=0.7, color='red', edgecolor='black')
axes1_tm[0].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes1_tm[0].set_xlabel('r-band: DELVE DR2 - DECaLS DR10', fontsize=12)
axes1_tm[0].set_ylabel('Number of sources', fontsize=12)
axes1_tm[0].set_xlim(-0.6, 0.6)
axes1_tm[0].grid(True, alpha=0.3)

# i-band
axes1_tm[1].hist(diff_i_delve_decals_tm_clean, bins=50, alpha=0.7, color='green', edgecolor='black')
axes1_tm[1].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes1_tm[1].set_xlabel('i-band: DELVE DR2 - DECaLS DR10', fontsize=12)
axes1_tm[1].set_ylabel('Number of sources', fontsize=12)
axes1_tm[1].set_xlim(-0.6, 0.6)
axes1_tm[1].grid(True, alpha=0.3)

# z-band
axes1_tm[2].hist(diff_z_delve_decals_tm_clean, bins=50, alpha=0.7, color='blue', edgecolor='black')
axes1_tm[2].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes1_tm[2].set_xlabel('z-band: DELVE DR2 - DECaLS DR10', fontsize=12)
axes1_tm[2].set_ylabel('Number of sources', fontsize=12)
axes1_tm[2].set_xlim(-0.6, 0.6)
axes1_tm[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('delve_decals_TM.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# Plot 2: DELVE - PS1 differences (TM)
# ============================================================================

fig2_tm, axes2_tm = plt.subplots(3, 1, figsize=(10, 12))
fig2_tm.suptitle('DELVE - PS1 Magnitude Differences (TM)', fontsize=16, fontweight='bold')

# r-band
axes2_tm[0].hist(diff_r_delve_ps1_tm_clean, bins=50, alpha=0.7, color='red', edgecolor='black')
axes2_tm[0].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes2_tm[0].set_xlabel('r-band: DELVE DR2 - PS1', fontsize=12)
axes2_tm[0].set_ylabel('Number of sources', fontsize=12)
axes2_tm[0].set_xlim(-0.6, 0.6)
axes2_tm[0].grid(True, alpha=0.3)

# i-band
axes2_tm[1].hist(diff_i_delve_ps1_tm_clean, bins=50, alpha=0.7, color='green', edgecolor='black')
axes2_tm[1].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes2_tm[1].set_xlabel('i-band: DELVE DR2 - PS1', fontsize=12)
axes2_tm[1].set_ylabel('Number of sources', fontsize=12)
axes2_tm[1].set_xlim(-0.6, 0.6)
axes2_tm[1].grid(True, alpha=0.3)

# z-band
axes2_tm[2].hist(diff_z_delve_ps1_tm_clean, bins=50, alpha=0.7, color='blue', edgecolor='black')
axes2_tm[2].axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
axes2_tm[2].set_xlabel('z-band: DELVE DR2 - PS1', fontsize=12)
axes2_tm[2].set_ylabel('Number of sources', fontsize=12)
axes2_tm[2].set_xlim(-0.6, 0.6)
axes2_tm[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('delve_ps1_TM.png', dpi=300, bbox_inches='tight')
plt.close()

print("\nTM plots saved successfully!")
print(f"  - delve_decals_TM.png")
print(f"  - delve_ps1_TM.png")
