#python3 rcsed2_coverage.py [optional_input_file.fits]
# If no argument provided, uses default file: final_TM.fits

from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from matplotlib import patheffects
import matplotlib.pyplot as plt
import sys
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS
import scipy.ndimage as ndimage

# Default input file path for band coverage plots
default_input_file = '/Users/dathev/PhDProjects/high-z-qso-selection/output_catalogs/results_catalogs/final_TM.fits'
# File for survey coverage plot
survey_input_file = '/Users/dathev/PhDProjects/high-z-qso-selection/output_catalogs/results_catalogs/initial_output_old_930.fits'

# Use command line argument if provided, otherwise use default
if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    input_file = default_input_file

# ========== Survey Coverage Plot ==========
print(f'Opening survey file: {survey_input_file}')
table_survey = fits.open(survey_input_file)
print('Survey file opened')
print(f'Total sources in survey file: {len(table_survey[1].data)}')

# Get all source coordinates for survey plot
all_ra = table_survey[1].data['ra']
all_dec = table_survey[1].data['dec']

# Filter out NaN values
valid_mask = ~(np.isnan(all_ra) | np.isnan(all_dec))
all_ra = all_ra[valid_mask]
all_dec = all_dec[valid_mask]
print(f'Valid sources (after NaN filter): {len(all_ra)}')

# Convert to SkyCoord and wrap for Aitoff
c_all = SkyCoord(ra=all_ra, dec=all_dec, frame=ICRS, unit=(u.deg, u.deg))
all_ra_rad = c_all.ra.wrap_at(180 * u.deg).radian
all_dec_rad = c_all.dec.radian

# Create single plot with all surveys overlaid
fig = plt.figure(figsize=(12, 6))
ax = plt.subplot(111, projection="aitoff")

print('Plotting surveys...')

# Create histogram once (all surveys use same data)
counts, xbins, ybins = np.histogram2d(-all_ra_rad, all_dec_rad, bins=300)
dd = ndimage.gaussian_filter(counts, sigma=3.0, order=0)
data = np.zeros((np.shape(dd)[0]+2, np.shape(dd)[1]+2))
data[1:-1, 1:-1] = dd

# Plot in order: WISE (back) -> DELVE (middle) -> VHS/VIKING (front)
# WISE - drawn first, with "//" hatching pattern (back layer) - lower opacity
cf1 = plt.contourf(data.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
                   levels=[0.1, 50000], colors='blue', alpha=0.08, zorder=1)
# Apply hatching to filled contours - access collections from axes
for i, collection in enumerate(ax.collections):
    if i < len(cf1.allsegs):  # Check if this is from cf1
        collection.set_hatch('//')
        collection.set_edgecolor('blue')
plt.contour(data.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
           levels=[0.1], colors='blue', linewidths=1.5, zorder=1)

# DELVE - drawn second, with "||" hatching pattern (middle layer) - medium opacity
cf2 = plt.contourf(data.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
                   levels=[0.1, 50000], colors='purple', alpha=0.20, zorder=2)
# Apply hatching - get the last collection added (should be from cf2)
n_collections_before = len(ax.collections)
for collection in ax.collections[n_collections_before-len(cf2.allsegs):n_collections_before]:
    collection.set_hatch('||')
    collection.set_edgecolor('purple')
plt.contour(data.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
           levels=[0.1], colors='purple', linewidths=1.5, zorder=2)

# VHS/VIKING - drawn last, with "--" hatching pattern (front layer) - higher opacity
cf3 = plt.contourf(data.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
                   levels=[0.1, 50000], colors='orange', alpha=0.30, zorder=3)
# Apply hatching - get the last collection added (should be from cf3)
n_collections_before = len(ax.collections)
for collection in ax.collections[n_collections_before-len(cf3.allsegs):n_collections_before]:
    collection.set_hatch('--')
    collection.set_edgecolor('orange')
plt.contour(data.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
           levels=[0.1], colors='orange', linewidths=2.0, zorder=3)

# Create legend with hatching patterns
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='purple', edgecolor='purple', hatch='||', label='DELVE', alpha=0.6),
    Patch(facecolor='orange', edgecolor='orange', hatch='--', label='VHS/VIKING', alpha=0.6),
    Patch(facecolor='blue', edgecolor='blue', hatch='//', label='WISE', alpha=0.6)
]

plt.legend(handles=legend_elements, loc='upper right', frameon=True, fontsize='large')
plt.grid(visible=True, color='k', linestyle=':', linewidth=1, zorder=100)

print('Saving...')
# plt.gca().set_xticklabels(['']*10)
plt.savefig('sky_rcsed_phot_v2.png', format='png')
plt.close()

# ========== Band Coverage Plots by Survey ==========
print('\n=== Creating Band Coverage Plots by Survey ===')
print(f'Opening band coverage file: {input_file}')

# Read data for band coverage (use final_TM.fits)
table = fits.open(input_file)
print('Band coverage file opened')
print(f'Total sources in band coverage file: {len(table[1].data)}')

# Read data for band coverage
data = table[1].data
ra_all = data['ra']
dec_all = data['dec']

# Filter out NaN coordinates
valid_coords = ~(np.isnan(ra_all) | np.isnan(dec_all))
ra_all = ra_all[valid_coords]
dec_all = dec_all[valid_coords]

# Helper function to create band coverage plot with overlaid bands
def plot_band_coverage(bands_info, survey_name, filename, total_sources, make_visible=False):
    """Create a single band coverage plot for a survey with all bands overlaid"""
    fig = plt.figure(figsize=(12, 6))
    ax = plt.subplot(111, projection="aitoff")
    
    # Create legend entries
    legend_elements = []
    
    # Plot each band on the same plot
    for band_name, mask, color in bands_info:
        if np.sum(mask) > 0:
            # Get coordinates for sources with this band
            ra_band_coords = ra_all[mask]
            dec_band_coords = dec_all[mask]
            
            # Convert to SkyCoord and wrap for Aitoff
            c_band = SkyCoord(ra=ra_band_coords, dec=dec_band_coords, frame=ICRS, unit=(u.deg, u.deg))
            ra_band_rad = c_band.ra.wrap_at(180 * u.deg).radian
            dec_band_rad = c_band.dec.radian
            
            # Create histogram and smooth
            counts, xbins, ybins = np.histogram2d(-ra_band_rad, dec_band_rad, bins=200)
            dd = ndimage.gaussian_filter(counts, sigma=2.5, order=0)
            data_smooth = np.zeros((np.shape(dd)[0]+2, np.shape(dd)[1]+2))
            data_smooth[1:-1, 1:-1] = dd
            
            # Adjust alpha and linewidth for better visibility if requested
            if make_visible:
                alpha_val = 0.25  # More opaque for better visibility
                linewidth_val = 2.5  # Thicker lines
            else:
                alpha_val = 0.15
                linewidth_val = 1.5
            
            # Plot filled contours
            plt.contourf(data_smooth.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
                        levels=[0.1, 50000], colors=color, alpha=alpha_val)
            # Plot contour lines
            plt.contour(data_smooth.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
                       levels=[0.1], colors=color, linewidths=linewidth_val)
            
            # Calculate percentage
            percentage = (np.sum(mask) / total_sources) * 100
            
            # Add to legend
            from matplotlib.lines import Line2D
            legend_elements.append(Line2D([0], [0], color=color, lw=2, 
                                         label=f'{band_name} ({percentage:.1f}%)'))
    
    plt.grid(visible=True, color='k', linestyle=':', linewidth=0.5, zorder=100)
    plt.legend(handles=legend_elements, loc='upper right', frameon=True, fontsize=12)
    plt.title(f'{survey_name} Band Coverage', fontsize=18, fontweight='bold', pad=20)
    
    plt.savefig(filename, format='png', dpi=300, bbox_inches='tight')
    print(f'Saved: {filename}')
    plt.close()

# Check which sources have detections in each band
# DELVE bands (optical) - check both delve and decals columns
# Try DELVE columns first, fall back to DECaLS if not available
if 'mag_r_delve' in data.columns.names:
    has_r = ~np.isnan(data['mag_r_delve'][valid_coords]) & (data['mag_r_delve'][valid_coords] != 99.0)
    has_i = ~np.isnan(data['mag_i_delve'][valid_coords]) & (data['mag_i_delve'][valid_coords] != 99.0)
    has_z = ~np.isnan(data['mag_z_delve'][valid_coords]) & (data['mag_z_delve'][valid_coords] != 99.0)
elif 'mag_r_decals' in data.columns.names:
    has_r = ~np.isnan(data['mag_r_decals'][valid_coords]) & (data['mag_r_decals'][valid_coords] != 99.0)
    has_i = ~np.isnan(data['mag_i_decals'][valid_coords]) & (data['mag_i_decals'][valid_coords] != 99.0)
    has_z = ~np.isnan(data['mag_z_decals'][valid_coords]) & (data['mag_z_decals'][valid_coords] != 99.0)
else:
    # Fallback to old column names
    has_r = ~np.isnan(data['mag_auto_r'][valid_coords]) & (data['mag_auto_r'][valid_coords] != 99.0)
    has_i = ~np.isnan(data['mag_auto_i'][valid_coords]) & (data['mag_auto_i'][valid_coords] != 99.0)
    has_z = ~np.isnan(data['mag_auto_z'][valid_coords]) & (data['mag_auto_z'][valid_coords] != 99.0)

# VHS bands (infrared)
has_y = ~np.isnan(data['ypetromag'][valid_coords]) & (data['ypetromag'][valid_coords] != -9.9999949E8)
has_j = ~np.isnan(data['jpetromag'][valid_coords]) & (data['jpetromag'][valid_coords] != -9.9999949E8)
has_h = ~np.isnan(data['hpetromag'][valid_coords]) & (data['hpetromag'][valid_coords] != -9.9999949E8)
has_ks = ~np.isnan(data['kspetromag'][valid_coords]) & (data['kspetromag'][valid_coords] != -9.9999949E8)

# WISE bands
has_w1 = ~np.isnan(data['w1mpro'][valid_coords]) & (data['w1mpro'][valid_coords] != 99.0)
has_w2 = ~np.isnan(data['w2mpro'][valid_coords]) & (data['w2mpro'][valid_coords] != 99.0)

total_valid = len(ra_all)

print(f'DELVE bands: R={np.sum(has_r)} ({100*np.sum(has_r)/total_valid:.1f}%), I={np.sum(has_i)} ({100*np.sum(has_i)/total_valid:.1f}%), Z={np.sum(has_z)} ({100*np.sum(has_z)/total_valid:.1f}%)')
print(f'VHS bands: Y={np.sum(has_y)} ({100*np.sum(has_y)/total_valid:.1f}%), J={np.sum(has_j)} ({100*np.sum(has_j)/total_valid:.1f}%), H={np.sum(has_h)} ({100*np.sum(has_h)/total_valid:.1f}%), Ks={np.sum(has_ks)} ({100*np.sum(has_ks)/total_valid:.1f}%)')
print(f'WISE bands: W1={np.sum(has_w1)} ({100*np.sum(has_w1)/total_valid:.1f}%), W2={np.sum(has_w2)} ({100*np.sum(has_w2)/total_valid:.1f}%)')

# Create DELVE plot (make more visible)
delve_bands = [
    ('DELVE R', has_r, 'green'),
    ('DELVE I', has_i, 'orange'),
    ('DELVE Z', has_z, 'red')
]
plot_band_coverage(delve_bands, 'DELVE', 'delve_band_coverage.png', total_valid, make_visible=True)

# Create VHS plot
vhs_bands = [
    ('VHS Y', has_y, 'purple'),
    ('VHS J', has_j, 'cyan'),
    ('VHS H', has_h, 'magenta'),
    ('VHS Ks', has_ks, 'yellow')
]
plot_band_coverage(vhs_bands, 'VHS', 'vhs_band_coverage.png', total_valid, make_visible=False)

# Create WISE plot (make more visible, no CatWISE)
wise_bands = [
    ('WISE W1', has_w1, 'blue'),
    ('WISE W2', has_w2, 'red')
]
plot_band_coverage(wise_bands, 'WISE', 'wise_band_coverage.png', total_valid, make_visible=True)

print('\nDone!')

# ========== Manual Survey Coverage Plot (No Input File) ==========
def plot_manual_survey_coverage():
    """
    Plot survey coverage based on manual specifications:
    - DELVE: -84° < dec < +18° AND |b| > 18° (Galactic latitude) - southern hemisphere
    - DECaLS: -18° < dec < +84° AND |b| > 18° (Galactic latitude)
    - VHS: dec < 0°, with gaps at dec ~ -5°, -15°, -30°, and -40°
    - VIKING: three equatorial strips along 0° declination (-5° to +5°): 
      Western (RA 210°-240°), Central (RA 330°-30°), Eastern (RA 120°-150°)
      Plus rectangular region: 0° ≲ RA ≲ 50° and -35° ≲ dec ≲ -10°
    - WISE: all-sky coverage
    - Initial catalog: coverage from initial_output_old_930.fits
    """
    print('\n=== Creating Manual Survey Coverage Plot ===')
    
    # Read initial catalog coverage
    print('Reading initial catalog file...')
    initial_file = '/Users/dathev/PhDProjects/high-z-qso-selection/output_catalogs/results_catalogs/initial_output_old_930.fits'
    table_initial = fits.open(initial_file)
    initial_ra = table_initial[1].data['ra']
    initial_dec = table_initial[1].data['dec']
    
    # Filter out NaN values
    valid_initial = ~(np.isnan(initial_ra) | np.isnan(initial_dec))
    initial_ra = initial_ra[valid_initial]
    initial_dec = initial_dec[valid_initial]
    print(f'Initial catalog sources: {len(initial_ra)}')
    
    # Generate synthetic coordinate arrays based on specifications
    # Use a fine grid for smooth coverage representation
    ra_grid = np.linspace(0, 360, 720)  # 0.5 degree resolution
    dec_grid = np.linspace(-90, 90, 360)  # 0.5 degree resolution
    
    # Create meshgrid
    RA, DEC = np.meshgrid(ra_grid, dec_grid)
    
    # Initialize coverage masks
    delve_mask = np.zeros_like(RA, dtype=bool)
    decals_mask = np.zeros_like(RA, dtype=bool)
    vhs_mask = np.zeros_like(RA, dtype=bool)
    viking_mask = np.zeros_like(RA, dtype=bool)
    wise_mask = np.ones_like(RA, dtype=bool)  # All-sky
    
    # Convert to Galactic coordinates to check |b| > 18° (used for both DELVE and DECaLS)
    all_coords = SkyCoord(ra=RA.flatten()*u.deg, dec=DEC.flatten()*u.deg, frame=ICRS)
    all_gal = all_coords.galactic
    b_gal = all_gal.b.deg.reshape(RA.shape)
    
    # DELVE: -84° < dec < +18° AND |b| > 18° (Galactic latitude)
    # Same shape as northern but mirrored to southern hemisphere
    delve_mask = (
        (DEC > -84.0) & (DEC < 18.0) &  # Celestial declination bounds (southern hemisphere)
        (np.abs(b_gal) > 18.0)  # Galactic latitude |b| > 18°
    )
    
    # DECaLS: -18° < dec < +84° AND |b| > 18° (Galactic latitude)
    decals_mask = (
        (DEC > -18.0) & (DEC < 84.0) &  # Celestial declination bounds
        (np.abs(b_gal) > 18.0)  # Galactic latitude |b| > 18°
    )
    
    # VHS: dec < 0°, with gaps at dec ~ -5°, -15°, -30°, and -40°
    # Create horizontal stripe gaps (with some width)
    vhs_base = (DEC < 0.0)
    gap_width = 2.0  # degrees
    gaps = (
        (DEC > -5.0 - gap_width/2) & (DEC < -5.0 + gap_width/2) |
        (DEC > -15.0 - gap_width/2) & (DEC < -15.0 + gap_width/2) |
        (DEC > -30.0 - gap_width/2) & (DEC < -30.0 + gap_width/2) |
        (DEC > -40.0 - gap_width/2) & (DEC < -40.0 + gap_width/2)
    )
    vhs_mask = vhs_base & ~gaps
    
    # VIKING: three equatorial strips along 0° declination (-5° to +5°)
    # Western strip: RA -150° to -120° (210° to 240°)
    # Central strip: RA -30° to 30° (330° to 30°, wrapping at 360°)
    # Eastern strip: RA 120° to 150°
    # Plus: rectangular region 0° ≲ RA ≲ 50° and -35° ≲ dec ≲ -10°
    viking_equatorial = (
        # Western strip: RA 210° to 240°
        ((RA >= 210.0) & (RA <= 240.0) & (DEC >= -5.0) & (DEC <= 5.0)) |
        # Central strip: RA 330° to 360° OR 0° to 30°
        (((RA >= 330.0) & (RA <= 360.0)) | ((RA >= 0.0) & (RA <= 30.0))) & 
        (DEC >= -5.0) & (DEC <= 5.0) |
        # Eastern strip: RA 120° to 150°
        ((RA >= 120.0) & (RA <= 150.0) & (DEC >= -5.0) & (DEC <= 5.0))
    )
    # Rectangular region filling VHS gap
    viking_rectangular = (
        (RA >= 0.0) & (RA <= 50.0) &
        (DEC >= -35.0) & (DEC <= -10.0)
    )
    viking_mask = viking_equatorial | viking_rectangular
    
    
    # Convert to coordinate arrays for plotting
    # Sample points from the masks to create coordinate lists
    np.random.seed(42)  # For reproducibility
    n_points_per_deg2 = 10  # Density of points
    
    def mask_to_coords(mask, ra_grid, dec_grid):
        """Convert a 2D mask to RA/Dec coordinate arrays"""
        ra_coords = []
        dec_coords = []
        for i in range(len(dec_grid)):
            for j in range(len(ra_grid)):
                if mask[i, j]:
                    # Add multiple points per pixel for density
                    n_local = max(1, int(n_points_per_deg2 * 0.5 * 0.5))  # ~0.5 deg resolution
                    for _ in range(n_local):
                        ra_coords.append(ra_grid[j] + np.random.uniform(-0.25, 0.25))
                        dec_coords.append(dec_grid[i] + np.random.uniform(-0.25, 0.25))
        return np.array(ra_coords), np.array(dec_coords)
    
    # Generate coordinates for each survey
    print('Generating DELVE coordinates...')
    delve_ra, delve_dec = mask_to_coords(delve_mask, ra_grid, dec_grid)
    
    print('Generating DECaLS coordinates...')
    decals_ra, decals_dec = mask_to_coords(decals_mask, ra_grid, dec_grid)
    
    print('Generating VHS coordinates...')
    vhs_ra, vhs_dec = mask_to_coords(vhs_mask, ra_grid, dec_grid)
    
    print('Generating VIKING coordinates...')
    viking_ra, viking_dec = mask_to_coords(viking_mask, ra_grid, dec_grid)
    
    print('Generating WISE coordinates...')
    wise_ra, wise_dec = mask_to_coords(wise_mask, ra_grid, dec_grid)
    
    # Create the plot
    fig = plt.figure(figsize=(12, 6))
    ax = plt.subplot(111, projection="aitoff")
    
    print('Plotting surveys...')
    
    # Helper function to plot a survey
    def plot_survey_coverage(ra_coords, dec_coords, color, hatch, alpha, zorder, label):
        if len(ra_coords) == 0:
            return None
        
        # Clip coordinates to valid ranges
        ra_coords = np.clip(ra_coords, 0, 360)
        dec_coords = np.clip(dec_coords, -90, 90)
        
        # Convert to SkyCoord and wrap for Aitoff
        c = SkyCoord(ra=ra_coords*u.deg, dec=dec_coords*u.deg, frame=ICRS)
        ra_rad = c.ra.wrap_at(180 * u.deg).radian
        dec_rad = c.dec.radian
        
        # Create histogram and smooth
        counts, xbins, ybins = np.histogram2d(-ra_rad, dec_rad, bins=300)
        dd = ndimage.gaussian_filter(counts, sigma=3.0, order=0)
        data = np.zeros((np.shape(dd)[0]+2, np.shape(dd)[1]+2))
        data[1:-1, 1:-1] = dd
        
        # Plot filled contours
        cf = plt.contourf(data.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
                         levels=[0.1, 50000], colors=color, alpha=alpha, zorder=zorder)
        
        # Apply hatching
        n_collections_before = len(ax.collections)
        for collection in ax.collections[n_collections_before-len(cf.allsegs):n_collections_before]:
            collection.set_hatch(hatch)
            collection.set_edgecolor(color)
        
        # Plot contour lines
        plt.contour(data.T, extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], 
                   levels=[0.1], colors=color, linewidths=1.5, zorder=zorder)
        
        return cf
    
    # Plot in order: WISE (back) -> DECaLS -> DELVE -> VHS -> VIKING -> Initial catalog (very front)
    plot_survey_coverage(wise_ra, wise_dec, 'blue', '//', 0.08, 1, 'WISE')
    plot_survey_coverage(decals_ra, decals_dec, 'cyan', '++', 0.18, 2, 'DECaLS')
    plot_survey_coverage(delve_ra, delve_dec, 'purple', '||', 0.20, 2.5, 'DELVE')
    plot_survey_coverage(vhs_ra, vhs_dec, 'orange', '--', 0.30, 3, 'VHS')
    
    # Plot VIKING as yellow filled region (same style as other surveys)
    print('Plotting VIKING coverage...')
    plot_survey_coverage(viking_ra, viking_dec, 'yellow', 'xx', 0.25, 4, 'VIKING')
    
    # Plot initial catalog coverage - very front layer
    print('Plotting initial catalog coverage...')
    plot_survey_coverage(initial_ra, initial_dec, 'green', '..', 0.35, 5, 'Initial output')
    
    # Create legend with hatching patterns - order: DELVE, DECaLS, VHS, WISE, VIKING, Initial catalog
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='purple', edgecolor='purple', hatch='||', label='DELVE', alpha=0.6),
        Patch(facecolor='cyan', edgecolor='cyan', hatch='++', label='DECaLS', alpha=0.6),
        Patch(facecolor='orange', edgecolor='orange', hatch='--', label='VHS', alpha=0.6),
        Patch(facecolor='blue', edgecolor='blue', hatch='//', label='WISE', alpha=0.6),
        Patch(facecolor='yellow', edgecolor='yellow', hatch='xx', label='VIKING', alpha=0.6),
        Patch(facecolor='green', edgecolor='green', hatch='..', label='Initial output', alpha=0.6)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right', frameon=True, fontsize='large')
    plt.grid(visible=True, color='k', linestyle=':', linewidth=1, zorder=100)
    
    print('Saving...')
    plt.savefig('sky_manual_coverage.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print('Saved: sky_manual_coverage.png')

# Uncomment the line below to generate the manual coverage plot
plot_manual_survey_coverage()
