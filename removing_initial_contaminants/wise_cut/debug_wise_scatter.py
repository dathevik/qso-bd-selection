import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import random

# Set random seed for reproducibility
random.seed(42)

# Read data
data_bd_wise = ascii.read('BD_wise_SpT.dat', encoding='latin-1')

# Brown Dwarfs list (Bañados+2016)
W1_bd_raw = data_bd_wise['w1mpro']
W2_bd_raw = data_bd_wise['w2mpro']
Sp_T_bd = data_bd_wise['optical_type']

# Convert to AB magnitudes
W1_bd = W1_bd_raw + 2.699
W2_bd = W2_bd_raw + 3.399

# Calculate W1-W2 color
W1_W2_bd = W1_bd - W2_bd

# Check for NaN values
print("=== DATA QUALITY CHECK ===")
print(f"Total BD sources: {len(W1_bd)}")
print(f"W1 NaN values: {np.sum(np.isnan(W1_bd))}")
print(f"W2 NaN values: {np.sum(np.isnan(W2_bd))}")
print(f"W1-W2 NaN values: {np.sum(np.isnan(W1_W2_bd))}")

# Remove NaN values for analysis
mask_valid = ~(np.isnan(W1_bd) | np.isnan(W2_bd) | np.isnan(W1_W2_bd))
W1_bd_clean = W1_bd[mask_valid]
W2_bd_clean = W2_bd[mask_valid]
W1_W2_bd_clean = W1_W2_bd[mask_valid]
Sp_T_bd_clean = Sp_T_bd[mask_valid]

print(f"Valid BD sources (no NaN): {len(W1_bd_clean)}")

# Check spectral type distribution
print("\n=== SPECTRAL TYPE DISTRIBUTION ===")
M_type_mask = np.char.startswith(np.array(Sp_T_bd_clean, dtype=str), 'M')
L_type_mask = np.char.startswith(np.array(Sp_T_bd_clean, dtype=str), 'L')
T_type_mask = np.char.startswith(np.array(Sp_T_bd_clean, dtype=str), 'T')
other_type_mask = ~M_type_mask & ~L_type_mask & ~T_type_mask

print(f"M-type BDs: {np.sum(M_type_mask)}")
print(f"L-type BDs: {np.sum(L_type_mask)}")
print(f"T-type BDs: {np.sum(T_type_mask)}")
print(f"Other types: {np.sum(other_type_mask)}")

# Calculate statistics for clean data
bd_in_box_clean = np.sum((W1_W2_bd_clean >= -0.6) & (W1_W2_bd_clean <= 0.6))
bd_percent_clean = (bd_in_box_clean / len(W1_W2_bd_clean)) * 100

print(f"\n=== SELECTION BOX STATISTICS (CLEAN DATA) ===")
print(f"BD sources in box (-0.6 to 0.6): {bd_in_box_clean} out of {len(W1_W2_bd_clean)} ({bd_percent_clean:.1f}%)")

# Check W1-W2 distribution
print(f"\n=== W1-W2 DISTRIBUTION ===")
print(f"W1-W2 range: {np.min(W1_W2_bd_clean):.3f} to {np.max(W1_W2_bd_clean):.3f}")
print(f"W1-W2 mean: {np.mean(W1_W2_bd_clean):.3f}")
print(f"W1-W2 median: {np.median(W1_W2_bd_clean):.3f}")

# Check distribution in different ranges
print(f"\n=== DISTRIBUTION BY RANGES ===")
print(f"W1-W2 < -0.6: {np.sum(W1_W2_bd_clean < -0.6)} ({np.sum(W1_W2_bd_clean < -0.6)/len(W1_W2_bd_clean)*100:.1f}%)")
print(f"-0.6 <= W1-W2 <= 0.6: {bd_in_box_clean} ({bd_percent_clean:.1f}%)")
print(f"W1-W2 > 0.6: {np.sum(W1_W2_bd_clean > 0.6)} ({np.sum(W1_W2_bd_clean > 0.6)/len(W1_W2_bd_clean)*100:.1f}%)")

# Create diagnostic plot
plt.figure(figsize=(12, 8))

# Plot all BD points
plt.scatter(W1_bd_clean, W1_W2_bd_clean, 
           label=f"BD Sample (n={len(W1_bd_clean)})", s=20, alpha=0.3, marker="o", color="green", zorder=1)

# Highlight points in selection box
mask_in_box = (W1_W2_bd_clean >= -0.6) & (W1_W2_bd_clean <= 0.6)
plt.scatter(W1_bd_clean[mask_in_box], W1_W2_bd_clean[mask_in_box], 
           label=f"BD in selection box (n={np.sum(mask_in_box)})", s=30, alpha=0.7, marker="o", color="red", zorder=2)

# Add WISE color cut lines
plt.axhline(y=-0.6, color='black', linestyle='--', linewidth=2, label='WISE color cut')
plt.axhline(y=0.6, color='black', linestyle='--', linewidth=2)

# Labels and formatting
plt.xlabel("W1 (AB)", fontsize=14)
plt.ylabel("W1-W2 (AB)", fontsize=14)
plt.title(f"BD WISE Color Distribution\n{bd_percent_clean:.1f}% in selection box", fontsize=16)
plt.legend(fontsize=12)
plt.ylim(-1, 3.5)
plt.xlim(14, 21)
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("debug_bd_wise_distribution.png", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

print(f"\nDiagnostic plot saved as 'debug_bd_wise_distribution.png'") 