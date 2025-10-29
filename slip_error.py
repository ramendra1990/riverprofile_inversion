# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 14:28:35 2025

@author: ramendra
"""

import math

# -----------------------
# Input Parameters
# -----------------------

# Terrace height data (meters)
h1 = 8.66      # Paleo terrace height
h1_err = 3.95  # Uncertainty in h1

h2 = 18.00     # Present-day terrace height
h2_err = 0.50  # Uncertainty in h2

# Dip angles (degrees)
delta1_deg = 6         # Dip of MHT (assumed constant)
delta2_deg = 14        # Dip of frontal ramp
delta2_err_deg = 4     # Uncertainty in delta2

# Radiocarbon age (years BP)
age_yr = 2356
age_err_yr = 32


# -----------------------
# Step 1: Calculate Slip
# -----------------------

# Convert degrees to radians
delta1_rad = math.radians(delta1_deg)
delta2_rad = math.radians(delta2_deg)
delta2_err_rad = math.radians(delta2_err_deg)

# Sines of dip angles
sin_delta1 = math.sin(delta1_rad)
sin_delta2 = math.sin(delta2_rad)
delta_sin = sin_delta1 - sin_delta2

# Elevation difference
delta_h = h1 - h2
delta_h_err = math.sqrt(h1_err**2 + h2_err**2)

# Slip (in meters)
slip = delta_h / delta_sin

# Monte Carlo propagation
import numpy as np
N = 200000  # number of Monte Carlo samples
rng = np.random.default_rng(1)
# Draw samples (assume Gaussian, independent)
h1_samp = rng.normal(h1, h1_err, size=N)
h2_samp = rng.normal(h2, h2_err, size=N)
delta1_samp = rng.normal(delta1_rad, 0, size=N)
delta2_samp = rng.normal(delta2_rad, delta2_err_rad, size=N)

sindiff = np.sin(delta1_samp) - np.sin(delta2_samp)
eps = 1e-4
frac_small = np.mean(np.abs(sindiff) < eps)
print(f"Fraction of samples with |sin_diff| < {eps:.1e}: {frac_small:.3%}")

# mask_bad = np.abs(sindiff) < eps
mask_bad = sindiff > 0
X = np.full(N, np.nan)
X[~mask_bad] = (h1_samp[~mask_bad] - h2_samp[~mask_bad]) / sindiff[~mask_bad]

# Summary statistics (robust)
valid = ~np.isnan(X)
print("N valid:", valid.sum())
median = np.median(X[valid])
p2p5, p97p5 = np.percentile(X[valid], [2.5, 97.5])
mean = np.mean(X[valid])
std = np.std(X[valid], ddof=1)

print(f"Median = {median:.3g}")
print(f"Mean ± std = {mean:.3g} ± {std:.3g}")
print(f"95% CI (percentiles) = [{p2p5:.3g}, {p97p5:.3g}]")

slip_sample = (h1_samp - h2_samp) / (np.sin(delta1_samp) - np.sin(delta2_samp))

# Summarize:
slip_mean = np.mean(slip_sample)
slip_median = np.median(slip_sample)
slip_std = np.std(slip_sample, ddof=1)
slip_lo, slip_hi = np.percentile(slip_sample, [2.5, 97.5])  # 95% CI

print(f"Mean = {slip_mean:.5f}, Median = {slip_median:.5f}, 1σ = {slip_std:.5f}")
print(f"95% CI = [{slip_lo:.5f}, {slip_hi:.5f}]")

# Error in sin(delta2)
cos_delta2 = math.cos(delta2_rad)
sin_delta2_err = cos_delta2 * delta2_err_rad

# Error in delta_sin (only delta2 has uncertainty)
delta_sin_err = sin_delta2_err

# Slip uncertainty
slip_err = abs(slip) * math.sqrt(
    (delta_h_err / delta_h)**2 + (delta_sin_err / delta_sin)**2
)

# -----------------------
# Step 2: Calculate Slip Rate (mm/yr)
# -----------------------

# Slip rate in mm/yr
slip_rate = (slip * 1000) / age_yr

# Slip rate uncertainty
slip_rate_err = slip_rate * math.sqrt(
    (slip_err / slip)**2 + (age_err_yr / age_yr)**2
)

# -----------------------
# Output
# -----------------------

print(f"Slip = {slip:.2f} ± {slip_err:.2f} m")
print(f"Slip Rate = {slip_rate:.2f} ± {slip_rate_err:.2f} mm/yr")
