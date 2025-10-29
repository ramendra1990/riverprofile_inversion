# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 09:46:26 2025

@author: ramendra
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
df = pd.read_excel(r'D:\ramendra_EPM103\Rajeeb\s357_PourPointResults_v3.xlsx', sheet_name='Sheet2')

# -----------------------------
# Example data (replace with yours)
ksn = df['Mean_Ksn'] # mean ksn
erosion = df['erosion_rate_mmky']   # mean erosion rate

sigmaE = df['error_E_mmky'] / 4 # uncertainty in erosion
# sigmaE = np.zeros((len(df), )) # uncertainty in erosion
sigmaksn = df['STD_Ksn'] # uncertainty in ksn
# sigmaksn = np.zeros((len(df), )) # uncertainty in erosion

# idx = df['Area_km2'] > 10
# ksn = ksn[idx]
# erosion = erosion[idx]
# sigmaE = sigmaE[idx]

# max_error = np.argmax(sigmaE/erosion) # error is more than half (50%) of the erosion rate.
# max_ksn = np.argmax(ksn)
# idx_del = [max_error]
# ksn = np.delete(ksn, idx_del)
# erosion = np.delete(erosion, idx_del)
# sigmaE = np.delete(sigmaE, idx_del)
# sigmaksn = np.delete(sigmaksn, idx_del)

# -----------------------------
# Bootstrap regression
n_boot = 100000
slopes, intercepts, r2_values = [], [], []

for _ in range(n_boot):
    # Resample erosion rate within uncertainty
    erosion_sample = erosion + np.random.normal(0, sigmaE)
    ksn_sample = ksn + np.random.normal(0, sigmaksn)
    
    # Linear fit (y = m*x + b)
    m, b = np.polyfit(ksn, erosion_sample, 1)
    slopes.append(m)
    intercepts.append(b)
    
    # Predicted values
    y_pred = m * ksn + b
    
    # Compute R²
    ss_res = np.sum((erosion_sample - y_pred) ** 2)
    ss_tot = np.sum((erosion_sample - np.mean(erosion_sample)) ** 2)
    r2 = 1 - ss_res/ss_tot
    r2_values.append(r2)

# Convert to arrays
slopes = np.array(slopes)
intercepts = np.array(intercepts)
r2_values = np.array(r2_values)

import seaborn as sns
# plt.figure()
# plt.subplot(2,2,1)
# sns.kdeplot(slopes)
# plt.subplot(2,2,2)
# sns.kdeplot(intercepts)
# plt.subplot(2,2,3)
# sns.kdeplot(r2_values)

# Best estimates
slope_mean, slope_std = np.mean(slopes), np.std(slopes)
intercept_mean, intercept_std = np.mean(intercepts), np.std(intercepts)
r2_mean, r2_std = np.mean(r2_values), np.std(r2_values)

print(f"Slope: {slope_mean:.3f} ± {slope_std:.3f}")
print(f"Intercept: {intercept_mean:.3f} ± {intercept_std:.3f}")
print(f"R²: {r2_mean:.3f} ± {r2_std:.3f}")

opt_m = 0.45 # From straghtaining of chi-plots
K = round(slope_mean / 1e2, 2) # (x10^(-4)) m^(1-2*opt_m)/yr
std_K = round(slope_std / 1e2, 2)
r2_score = round(r2_mean, 1)
# -----------------------------
# Plot regression results. Main plot to be exported for manuscript.
mpl.rcParams['svg.fonttype'] = 'none' # For saving figures as .svg. Most compatble with inkscape
mpl.rcParams['font.family'] = 'Arial'  # or another common system font

# plt.figure(figsize=(5,4))
plt.figure(figsize=(3,2.5))
# Plot observed data with error bars
plt.errorbar(ksn, erosion, xerr = sigmaksn, yerr=sigmaE, fmt='none', capsize=4,
             lw = 1, color="black", alpha = 0.2)
plt.scatter(ksn, erosion, marker='o', s=35, facecolors='none', edgecolors = 'k')
# Define x-range for fit line
x_fit = np.linspace(min(ksn)-1, max(ksn)+1, 200)

# Best-fit line
y_fit = slope_mean * x_fit + intercept_mean
# plt.plot(x_fit, y_fit, 'r-', label=f"$K={K} \pm {std_K} \times 10^{{-4}}$")
# plt.plot(x_fit, y_fit, 'r-', label=r"$K=\left(%s \pm %s \right) \times 10^{-4} \ m^{0.7}/yr$" "\n" r"$R^{2}=%s$" %(K, std_K, r2_score))
plt.plot(x_fit, y_fit, 'k-')

# 95% confidence interval from bootstrap
y_boot = np.array([m * x_fit + b for m, b in zip(slopes, intercepts)])
y_lower = np.percentile(y_boot, 2.5, axis=0)
y_upper = np.percentile(y_boot, 97.5, axis=0)
plt.fill_between(x_fit, y_lower, y_upper, color='k', alpha=0.15)

# insert a text with slope and R2 values
text = r"$K=\left(%s \pm %s \right) \times 10^{-4} \ m^{%s}/yr$" "\n" r"$R^{2}=%s$" %(K, std_K, round(1-2*opt_m, 1), r2_score)
plt.text(0.05, 0.80, text, transform=plt.gca().transAxes, color="k", fontsize = 8)

# Labels
plt.xlabel(r"$K_{sn} (m^{%s})$" %(round(2*opt_m, 1)), fontsize = 12)
plt.ylabel(r"$E (m \ Myr^{-1})$", fontsize = 12)
plt.xticks(fontsize=12)   # x-axis tick labels
plt.yticks(fontsize=12)
# plt.title("Erosion rate vs. ksn (bootstrap regression)")
# plt.legend()
# plt.grid(True)
plt.tight_layout()

# Export the plot
plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\EvsKsn.svg')
# plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\EvsKsn.pdf')
# plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\EvsKsn.png', dpi = 600)

# -----------------------------
# # Plot histogram of R² values
# plt.figure(figsize=(6,4))
# plt.hist(r2_values, bins=30, color='skyblue', edgecolor='black')
# plt.axvline(r2_mean, color='red', linestyle='--', label=f"Mean R² = {r2_mean:.2f}")
# plt.xlabel("R²")
# plt.ylabel("Frequency")
# plt.title("Bootstrap R² distribution")
# plt.legend()
# plt.tight_layout()
# plt.show()
