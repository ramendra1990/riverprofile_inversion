#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 17:56:27 2023

@author: cds
"""

def fdir_to(idx, idy, fdir, method = 'arcgis'):
    methods = ['qgis', 'arcgis']
    if method not in methods:
        raise ValueError("Invalid sim type. Expected one of: %s" % methods)
    
    if method == 'qgis':
        if fdir == 0:
            idx += 0
            idy += -1
        elif fdir == 1:
            idx += 1
            idy += -1
        elif fdir == 2:
            idx += 1
            idy += 0
        elif fdir == 3:
            idx += 1
            idy += 1
        elif fdir == 4:
            idx += 0
            idy += 1
        elif fdir == 5:
            idx += -1
            idy += 1
        elif fdir == 6:
            idx += -1
            idy += 0
        elif fdir == 7:
            idx += -1
            idy += -1
    else:
        if fdir == 64:
            idx += 0
            idy += -1
        elif fdir == 128:
            idx += 1
            idy += -1
        elif fdir == 1:
            idx += 1
            idy += 0
        elif fdir == 2:
            idx += 1
            idy += 1
        elif fdir == 4:
            idx += 0
            idy += 1
        elif fdir == 8:
            idx += -1
            idy += 1
        elif fdir == 16:
            idx += -1
            idy += 0
        elif fdir == 32:
            idx += -1
            idy += -1
    
    return idy, idx

def CreateChi(str_subset, facc_subset, fdir_subset, 
              dem_subset, X_coord_2d, Y_coord_2d, m, 
              flow_method = 'qgis'):
    methods = ['qgis', 'arcgis']
    if flow_method not in methods:
        raise ValueError("Invalid sim type. Expected one of: %s" % methods)
        
    import numpy as np
    
    ind_1 = np.where(str_subset == 1)
    cost_array_1d = np.zeros(ind_1[0].shape)
    # str_1d = str_subset.ravel()
    facc_1d = facc_subset[ind_1]    
    dem_1d = dem_subset[ind_1]
    
    X_coord_1d = X_coord_2d[ind_1]
    Y_coord_1d = Y_coord_2d[ind_1]
    
    nrows = str_subset.shape[0]
    ncols = str_subset.shape[1]
    
    node_mat = np.zeros((len(ind_1[0]), 9))
    node_mat[:] = np.nan
    
    # Calculation of the 1d indices for the stream pixels
    for i in range(len(ind_1[0])):
        idx = ind_1[1][i]
        idy = ind_1[0][i]
        
        node_mat[i, 0] = idy * ncols + idx
            
        if np.isnan(fdir_subset[idy, idx]):
            node_mat[i, 1] = np.nan
            
        else:
            fdir_id = fdir_to(idx, idy, fdir_subset[idy, idx], method = flow_method)
            
            if fdir_id[0] < 0 or fdir_id[0] >= nrows or fdir_id[1] < 0 or fdir_id[1] >= ncols:
                node_mat[i, 1] = np.nan
            elif facc_subset[fdir_id] < facc_subset[idy, idx] or np.isnan(facc_subset[fdir_id]):
                node_mat[i, 1] = np.nan
            elif fdir_id[0] >= 0 and fdir_id[0] < nrows and fdir_id[1] >= 0 and fdir_id[1] < ncols:
                node_mat[i, 1] = fdir_id[0] * ncols + fdir_id[1]
            
        print(i)
    
    node_mat[:, 2] = facc_1d
    node_mat[:, 3] = cost_array_1d
    node_mat[:, 4] = X_coord_1d
    node_mat[:, 5] = Y_coord_1d
    node_mat[:, 8] = dem_1d
    
    while np.sum(np.isnan(node_mat[:, 7])) > 0:
        indx_node = np.argwhere(node_mat[:, 2] == np.max(node_mat[np.isnan(node_mat[:, 7]), 2]))
        
        if len(indx_node) > 1:
            for i in indx_node:
                if np.isnan(node_mat[i, 7]):
                    indx_new = i
            indx_node = indx_new
            
            
        if np.isnan(node_mat[indx_node[0], 1]):
            node_mat[indx_node[0], 3] = 0
            node_mat[indx_node[0], 6] = 0
            node_mat[indx_node[0], 7] = 0
        else:
            target_idx = np.argwhere(node_mat[:, 0] == node_mat[indx_node[0], 1])
            x_j = node_mat[target_idx[0], 4]
            y_j = node_mat[target_idx[0], 5]
            dist_j = node_mat[target_idx[0], 3]
            chi_j = node_mat[target_idx[0], 7]
            x_i = node_mat[indx_node[0], 4]
            y_i = node_mat[indx_node[0], 5]
            dist_ij = np.sqrt(((x_i - x_j) ** 2) + ((y_i - y_j) ** 2))
            Area_i = node_mat[indx_node[0], 2]
            chi_i = chi_j + (dist_ij / (Area_i ** m))
            dist_i = dist_j + dist_ij
            node_mat[indx_node[0], 6] = dist_ij
            node_mat[indx_node[0], 7] = chi_i
            node_mat[indx_node[0], 3] = dist_i
            
        x = node_mat[indx_node[0], 4]
        y = node_mat[indx_node[0], 5]
        chi = node_mat[indx_node[0], 7]
        dist = node_mat[indx_node[0], 3]
        
        while len(indx_node) > 0:
            indx_node = np.argwhere(node_mat[:, 1] == node_mat[indx_node[0], 0])
            if len(indx_node) > 0:
                x_i = node_mat[indx_node[0], 4]
                y_i = node_mat[indx_node[0], 5]
                dist_ij = np.sqrt(((x_i - x) ** 2) + ((y_i - y) ** 2))
                Area_i = node_mat[indx_node[0], 2]
                chi_i = chi + (dist_ij / (Area_i ** m))
                dist_i = dist + dist_ij
                node_mat[indx_node[0], 3] = dist_i
                node_mat[indx_node[0], 6] = dist_ij
                node_mat[indx_node[0], 7] = chi_i
                x = node_mat[indx_node[0], 4]
                y = node_mat[indx_node[0], 5]
                chi = node_mat[indx_node[0], 7]
                dist = node_mat[indx_node[0], 3]
                # print(indx_node[0], len(indx_node))
            else:
                pass
            
        print('/////-----------------\\\\\\')    
        print(np.sum(np.isnan(node_mat[:, 7])))
        print('/////-----------------\\\\\\')

    return node_mat

def trunk(str_subset, fdir_subset, facc_subset):
    
    import numpy as np
    from scipy.ndimage import convolve
    
    endpoint_kernel = np.uint8([[1, 1, 1],[1,10,1],[1,1,1]])
    str_filtered = convolve(str_subset.astype(np.uint8), endpoint_kernel, mode='nearest')
    endpoints_i, endpoints_j = np.where(str_filtered == 11)
    
    cost_list = []
    path_list = []
    
    indx = np.where(facc_subset == np.nanmax(facc_subset))
    ind_outlet = (indx[0][0], indx[1][0])
    
    for i in range(len(endpoints_i)):
        start_point = (endpoints_i[i], endpoints_j[i])
        path = [start_point]
        cost = 0
        while not (start_point == ind_outlet):
            idx = start_point[1]
            idy = start_point[0]
            start_point = fdir_to(idx, idy, fdir_subset[idy, idx], method = 'arcgis')
            if np.log(fdir_subset[idy, idx]) / np.log(2) % 2 == 0:
                cost += 1
            else:
                cost += 1.41
            path.append(start_point)
        
        path_list.append(path)
        cost_list.append(cost)
     
    indx = np.argmax(np.array(cost_list))
    trunk_path = path_list[indx]
    str_trunk = str_subset * 0
    for p in trunk_path:
        str_trunk[p[0], p[1]] = 1
        
    return str_trunk


def smoothIT(ar, s, w):
    t = (((w - 1) / 2) - 0.5) / s
    import numpy as np
    from scipy import ndimage
    ar1 = ar.copy()
    ar1[np.isnan(ar1)] = 0
    ar_smooth = ndimage.gaussian_filter(ar1, sigma = s, truncate = t, mode = 'constant')
    ar_smooth[np.isnan(ar)] = np.nan
    return ar_smooth

# %%
from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['svg.fonttype'] = 'none'   # keep text as text
mpl.rcParams['font.family'] = 'Arial'  # or another common system font
from scipy import ndimage
# import skimage.graph

# %% Data read

# Foor mohand anticline
str_in = gdal.Open(r'D:\ramendra_EPM103\Rajeeb\Data\S_clipped.tif')
facc_in = gdal.Open(r'D:\ramendra_EPM103\Rajeeb\Data\A_clipped.tif')
fdir_in = gdal.Open(r'D:\ramendra_EPM103\Rajeeb\Data\FD_clipped.tif')
dem_in = gdal.Open(r'D:\ramendra_EPM103\Rajeeb\Data\DEM_clipped_filled.tif')

# Get georeferencing info
gt = str_in.GetGeoTransform()
xmin = gt[0]
xres = gt[1]
ymax = gt[3]
yres = gt[5]
xmax = xmin + str_in.RasterXSize * xres
ymin = ymax + str_in.RasterYSize * yres
extent = [xmin, xmax, ymin, ymax]  # extent in map coordinates

"""
# For havelock
str_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/Str11Sm_1e6.tif')
facc_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/DEM11Sm_facc_m2.tif')
fdir_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/DEM11Sm_fdir.tif')
dem_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/DEM11Sm_fill.tif')
# slope_in = gdal.Open(r'E:\paper_underPreparation\Active\Saikat\Slope11_mm.tif')

# # For henry
# str_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/Str_henry_clipSM_1e6.tif')
# facc_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/henry_clipSM_facc.tif')
# fdir_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/henry_clipSM_fdir.tif')
# dem_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/henry_clipSM_fill.tif')

# For south andaman
str_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/Str_South_IslandSM_1e6.tif')
facc_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/South_andamanSM_facc.tif')
fdir_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/South_andamanSM_fdir.tif')
dem_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files3/South_andamanSM_fill.tif')

# # ttbox files
# str_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/ttbox_files/S11_1e6.tif')
# facc_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/ttbox_files/DA_grid11.tif')
# fdir_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/ttbox_files/FLOWobj11.tif')
# dem_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/ttbox_files/filled_DEM11.tif')
# slope_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/processed_files2/Saikat/Slope11_mm.tif')

# str_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/WG_processed2/WG_clipped1_str5e6.tif')
# facc_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/WG_processed2/WG_clipped1_facc.tif')
# fdir_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/WG_processed2/WG_clipped1_fdir.tif')
# dem_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/WG_processed2/WG_clipped1_fill.tif')
# slope_in = gdal.Open('/home/cds/Desktop/ramendra/saikat/WG_processed2/WG_clipped1_slp_mm.tif')
"""

# %% Data process
facc_in_re = gdal.Warp('', facc_in, xRes = str_in.GetGeoTransform()[1], 
                        yRes = -str_in.GetGeoTransform()[5], format = 'VRT')
fdir_in_re = gdal.Warp('', fdir_in, xRes = str_in.GetGeoTransform()[1], 
                        yRes = -str_in.GetGeoTransform()[5], format = 'VRT')
dem_in_re = gdal.Warp('', dem_in, xRes = str_in.GetGeoTransform()[1], 
                        yRes = -str_in.GetGeoTransform()[5], format = 'VRT')
# slope_in_re = gdal.Warp('', slope_in, xRes = str_in.GetGeoTransform()[1], 
#                         yRes = -str_in.GetGeoTransform()[5], format = 'VRT')

str_array = str_in.ReadAsArray()
facc_array = facc_in_re.ReadAsArray()
fdir_array = fdir_in_re.ReadAsArray()
dem_array = dem_in_re.ReadAsArray()
# slope_array = slope_in_re.ReadAsArray()

str_array = str_array.astype('float')
str_array[str_array < 0] = np.nan
facc_array[facc_array < 0] = np.nan
fdir_array = fdir_array.astype('float')
fdir_array[(fdir_array <= 0) | ((fdir_array > 128))] = np.nan
dem_array = dem_array.astype('float')
dem_array[(dem_array <= -9999) | ((dem_array > 9999))] = np.nan
# slope_array[slope_array == slope_array.min()] = np.nan



# plt.figure()
# plt.subplot(221)
# plt.imshow(str_array)
# plt.subplot(222)
# plt.imshow(fdir_array)
# plt.subplot(223)
# plt.imshow(facc_array ** 0.2)
# plt.subplot(224)
# plt.imshow(dem_array)


structure = np.ones((3, 3)).astype('int')
ar1 = str_array.copy()
ar1[np.isnan(ar1)] = 0 # only zeros and ones in the array
labeled, ncomponents = ndimage.label(ar1, structure) # Number of connected components in the stream network

sizes = ndimage.sum(ar1, labeled, range(ncomponents + 1))
k = len(sizes) - 1

# Make dictionary with sorted sequence of the connected compnents as the keys
CC_dic = dict.fromkeys(range(1, k + 1))
for i in range(1, k + 1):
    indx = np.argsort(sizes, axis=0)[-i]
    CC = labeled == indx
    CC_dic[i] = CC
    
largestCC = labeled == np.argmax(sizes)

ulx, xres, xskew, uly, yskew, yres  = str_in.GetGeoTransform()
lrx = ulx + (str_in.RasterXSize * xres)
lry = uly + (str_in.RasterYSize * yres)
X = np.linspace(ulx, lrx, num = str_in.RasterXSize)
Y = np.linspace(uly, lry, num = str_in.RasterYSize)

X_coord_2d = np.zeros(str_array.shape)
for i in range(str_array.shape[0]):
    X_coord_2d[i, :] = X
    
Y_coord_2d = np.zeros(str_array.shape)
for i in range(str_array.shape[1]):
    Y_coord_2d[:, i] = Y

#%% CC_dic has all the connected components. User has to decide how many components need to be considered for the inversion. Accordingly, all the input arrays have to be modified as follows. Here, for example, we are putting a constraint that the size of the component has to be larger than 10, which can be changed according to the scenario at hand.
# Following are the line of codes for modification of the arrays based on the sizes of the connected components
# Note: The constraints can be put in terms of number of components also
str_subset = 0
fdir_subset = 0
facc_subset = 0
dem_subset = 0
slope_subset = 0
X_coord_2d_subset = 0
Y_coord_2d_subset = 0

id_list = [2, 3, 5, 6, 7, 8, 9, 19, 39, 68, 77]
# --------- For a specific DN (for instance, in Mohand anticline, khajnawar basin is the 9th largest across the MA)
trunk_stream = 'n'
i = 6
cc = CC_dic[id_list[i]]
print(int(np.nansum(cc)))
# To modify DEM for each connected component. Minimum elevation has to be zero
F = facc_array * cc # We will use maximum flow accumulation to find out the base level
S = str_array * cc
if trunk_stream == 'y':
    S = trunk(S, fdir_array, F)
D = dem_array * cc
r, c = np.where(F == np.nanmax(F))
D[S == 1] = D[S == 1] - D[r, c]
str_subset += S
fdir_subset += fdir_array * cc
facc_subset += F
dem_subset += D
X_coord_2d_subset += X_coord_2d * cc
Y_coord_2d_subset += Y_coord_2d * cc

plt.figure()
# plt.imshow(dem_array, cmap='grey')
plt.imshow(str_subset, extent=extent, origin="upper")
plt.figure()
plt.subplot(221)
plt.imshow(str_subset)
plt.subplot(222)
plt.imshow(fdir_subset)
plt.subplot(223)
plt.imshow(facc_subset ** 0.2)
plt.subplot(224)
plt.imshow(dem_subset)

# import geopandas as gpd
# gdf = gpd.read_file(r'D:\ramendra_EPM103\Rajeeb\s357_CRNresults_Be10_pourpoints.shp')
# gdf.plot(marker = '+', ax = plt.gca())

 
#%% 
m_try = np.arange(0.0, 0.8, 0.05)  

nbins = 25
scat_metric = np.zeros((len(m_try),))
plt.figure()
for i, m in enumerate(m_try):
    node_mat = CreateChi(str_subset, facc_subset, fdir_subset, dem_subset, 
                          X_coord_2d_subset, Y_coord_2d_subset, m, 
                          flow_method='arcgis')
    plt.subplot(4, 4, i + 1)
    plt.scatter(node_mat[:, 7], node_mat[:, 8], s = 1, 
                marker = '.', alpha = 0.8)
    plt.tick_params(axis = 'both', which = 'both', bottom = False, 
                    top = False, left = False, right = False, 
                    labelbottom = False, 
                    labelleft = False)
    textstr = 'm = %.2f' % (m_try[i])
    ax = plt.gca()
    plt.text(0.1, 0.8, textstr, fontsize = 8, 
              transform = ax.transAxes)
    
    chi_sorted = node_mat[node_mat[:, 7].argsort()]
    val_in_bins = int(np.floor(len(chi_sorted) / nbins))
    chi_scat = np.zeros((nbins, 1))
    for k in range(1, nbins):
        chi_scat[k, 0] = np.std(chi_sorted[val_in_bins * (k - 1) : val_in_bins * k, 8])
        
    chi_scat[-1] = np.std(chi_sorted[val_in_bins * (nbins-1):,8])
    scat_metric[i] = np.mean(chi_scat)

plt.subplots_adjust(wspace = 0.1, hspace = 0.1)
plt.tight_layout()

# plt.figure(figsize=(5, 4))    
plt.figure(figsize=(3, 2.5))   
plt.plot(m_try, scat_metric, 'ko-')
plt.xlabel('m', fontsize = 12)
# plt.ylabel('mean of z scatter in bins', fontsize = 12)
plt.ylabel(r'$\langle |\Delta z| \rangle$', fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.tight_layout()

plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\opt_m.svg')

# print(np.mean(m_try[np.argmin(scat_metric)]))

# -----------------------------------------------------------------------------
from scipy import stats
opt_m = 0.45
node_mat = CreateChi(str_subset, facc_subset, fdir_subset, dem_subset, 
                      X_coord_2d_subset, Y_coord_2d_subset, opt_m, 
                      flow_method='arcgis')
X, Y = node_mat[:, 7], node_mat[:, 8]
slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y)
X_fit = np.linspace(min(X), max(X), 100)
Y_fit = slope * X_fit + intercept

plt.figure(figsize=(5, 4))
plt.scatter(X, Y, label="Data", edgecolors = 'k', facecolors = 'none', alpha = 0.1)
plt.plot(X_fit, Y_fit, color="k", 
         label=f"K={slope:0.3f} \n $\\Delta$k={std_err:0.3f} \n R^2={r_value**2:0.2f}")
plt.xlim(0,)
plt.ylim(0,)
plt.xlabel("x")
plt.ylabel("y")
plt.legend()

# %% For bootstrapping error. there is whole new routine we have to write
# We need a couple more function
def fit_once(tau, z, q = 10):
    # Returns the optimal theta (uplft rate, Uq) and the fitted elevation values
    tau_mat = np.vstack((tau, z)).T
    tau_mat_sorted = tau_mat[tau_mat[:, 0].argsort()] # Sort wrt to tau values
    tau = tau_mat_sorted[:, 0]
    z = tau_mat_sorted[:, 1] - np.min(tau_mat_sorted[:, 1])

    ind_nz = z != 0
    z = z[ind_nz]
    tau = tau[ind_nz]

    N = len(z)

    # Define q, i.e. the number of delta_t
    rem = N % q
    Nq = int((N - rem) / q) # Number of data points in each time division

    t_vec = np.zeros((q+1, ))
    t_vec[0] = 0
    indx = Nq - 1
    for i in range(1, q):
        t_vec[i] = tau[int(indx)]
        indx += Nq
        
    t_vec[q] = tau[-1 - rem]

    delta_t = np.diff(t_vec)

    # Define R
    R = np.zeros((N - rem,))
    for i in range(0, N - rem):
        j = int(np.floor(i / Nq))
        if j == 0:
            R[i] = tau[i]
        else:
            delta_t_sum = 0
            for a in range(j):
                delta_t_sum += delta_t[a]
            R[i] = tau[i] - delta_t_sum
            
    # Define A matrix (Forward operator matrix)
    if rem == 0:
        N1 = N
    else:
        N1 = N + 1 - rem
        
    A = np.zeros((N1, q))
    for i in range(N1):
        j = int(np.floor(i / Nq))
        if j == 0:
            A[i, j] = R[i]
        else:
            A[i, 0 : j] = delta_t[0 : j]
            if j < q:
                A[i, j] = R[i]
                
    AT = A.T
    ATA = np.dot(AT, A)

    # Linear regression for calculation of U vector
                
    # Define apriori for Uplift rate
    U_pri = np.zeros((q,))
    # sum_r = 0

    sum_U_pri = 0
    for i in range(N1):
        sum_tau = np.sum(A[i, :])
        if sum_tau > 0:
            sum_U_pri = sum_U_pri +  (z[i] / sum_tau)
            
    U_pri = np.ones((q, )) * sum_U_pri / N1


    epsilon = 1e3
    ATA_damped_inv = np.linalg.inv(ATA + ((epsilon ** 2) * np.identity(q)))
    ATA_damped_inv_AT = np.dot(ATA_damped_inv, AT)
    U = U_pri + np.dot(ATA_damped_inv_AT, z[0 : N1] - np.dot(A, U_pri)) # Model parameter
    # Non-negative uplift not allowed. Positivity constraint
    # U[U < 0] = 0 
    z_inf = np.dot(A, U) # Fitted results

    # Time capsule and plotting of U vs T
    T = np.zeros((q + 1,))
    for i in range(1, q + 1):
        T[i] = T[i - 1] + delta_t[i - 1]
        
    U1 = np.zeros((q + 1,))
    U1[0 : -1] = U
    U1[-1] = U[-1]

    ur_data = np.vstack((T, U1)).T
    
    z_transformed, tau_transformed = z[0: N1], tau[0: N1]
    return ur_data, z_inf, z_transformed, tau_transformed

# Mammen wild multipliers (good choice)
def mammen(n, rng):
    u = rng.random(n)
    a = (np.sqrt(5) + 1)/2
    b = (np.sqrt(5) - 1)/2
    v = np.where(u < (np.sqrt(5) - 1)/(2*np.sqrt(5)),
                 -a, b)  # mean 0, var 1
    return v

# %%
from matplotlib.lines import Line2D # For extra legend handling
import matplotlib.patches as mpatches
import pandas as pd

# """
opt_m = 0.45
node_mat = CreateChi(str_subset, facc_subset, fdir_subset, dem_subset,
                     X_coord_2d_subset, Y_coord_2d_subset, opt_m,
                     flow_method='arcgis')
df1 = pd.DataFrame({'X' : node_mat[:, 4], 
                   'Y' : node_mat[:, 5], 
                   'Elevation' : node_mat[:, 8] - node_mat[:, 8].min(), 
                   'Distance' : node_mat[:, 3], 
                   'chi' : node_mat[:, 7]})
z = df1['Elevation']

# For different ke
plt.figure(figsize=(9, 3))

opt_ke = 1.26e-4 # For khajnawar basin (mohand anticline)
tau = df1['chi'] / opt_ke
theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
plt.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
         color = 'gray', ls = ':', alpha = 0.3)

opt_ke = 1.43e-4 # For khajnawar basin (mohand anticline)
tau = df1['chi'] / opt_ke
theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
plt.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
         color = 'b', ls = '-', alpha = 1)

opt_ke = 1.60e-4 # For khajnawar basin (mohand anticline)
tau = df1['chi'] / opt_ke
theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
plt.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
         color = 'k', ls = '--', alpha = 0.3)

# Plot monsoon data from Prell and Kutzbach, 1987
P_df = pd.read_excel(r'D:\ramendra_EPM103\Rajeeb\terrace_dates.xlsx', sheet_name='1987_Prell')
ax = plt.gca()
ax1 = ax.twinx()
idx = P_df['time'] <= 90
ax1.plot(P_df.loc[idx, 'time'], P_df.loc[idx, 'del_P(0-30)'], 'orange', alpha = 0.3)


ax.set_xlim(0,)
ax.set_xlabel("Time (ka)", fontsize = 15)
ax.set_ylabel("Uplift rate (mm/yr)", fontsize = 15)
ax1.set_ylabel(r'$\Delta P \ (\%)$', fontsize = 15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax1.tick_params(axis='both', which='major', labelsize=15)
# ax.set_yticks(fontsize=15)
# ax1.set_xticks(fontsize=15)
# ax1.set_yticks(fontsize=15)

# legend_elements = [
#     Line2D([0], [0], color='gray', ls = ':', lw=1, label=r'$K=1.26\times 10^{-4} \ m^{0.1}/yr$'),
#     Line2D([0], [0], color='b', lw=1, ls = '-', label=r'$K=1.43\times 10^{-4} \ m^{0.1}/yr$'),
#     Line2D([0], [0], color='k', lw=1, ls = '--', label=r'$K=1.60\times 10^{-4} \ m^{0.1}/yr$')
# ]
legend_elements = [
    Line2D([0], [0], color='gray', ls = ':', lw=1, label=r'$U_{K-\Delta K}$'),
    Line2D([0], [0], color='b', lw=1, ls = '-', label=r'$U_{K}$'),
    Line2D([0], [0], color='k', lw=1, ls = '--', label=r'$U_{K+\Delta K}$'), 
    Line2D([0], [0], color='orange', lw=1, ls = '-', label=r'$\Delta P$')
]
plt.legend(handles=legend_elements, fontsize=10)


plt.tight_layout()

# Save the sensitivity to K figure
plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\K-sensitivity.svg')
plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\K-sensitivity&climate.svg')

# --------- a new figure ------------------
fig, axes = plt.subplots(
    2, 1,               # 2 rows, 1 column
    figsize=(9,3),
    gridspec_kw={'height_ratios': [1, 2]}  # top subplot twice the height of bottom
)

ax = axes[1]

# mean_k = 1.43e-4
# sd_k = 0.17e-4

# for i in range(2000):
#     opt_ke = np.random.normal(loc = 1.43e-4, scale=0.17e-4)
#     if (opt_ke >= mean_k - sd_k) & (opt_ke <= mean_k + sd_k):
#         tau = df1['chi'] / opt_ke
#         theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
#         ax.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
#                  color = 'grey', ls = '-', alpha = 0.01)
    
# opt_ke = 1.43e-4 # For khajnawar basin (mohand anticline)
# tau = df1['chi'] / opt_ke
# theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
# ax.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
#          color = 'b', ls = '-', alpha = 1)
    
opt_ke = 1.26e-4 # For khajnawar basin (mohand anticline)
tau = df1['chi'] / opt_ke
theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
ax.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
         color = 'gray', ls = ':', alpha = 0.3)

opt_ke = 1.43e-4 # For khajnawar basin (mohand anticline)
tau = df1['chi'] / opt_ke
theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
ax.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
         color = 'b', ls = '-', alpha = 1)

opt_ke = 1.60e-4 # For khajnawar basin (mohand anticline)
tau = df1['chi'] / opt_ke
theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
ax.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
         color = 'k', ls = '--', alpha = 0.3)

ax1 = ax.twinx()
idx = P_df['time'] <= 90
ax1.plot(P_df.loc[idx, 'time'], P_df.loc[idx, 'del_P(0-30)'], 'r', alpha = 0.5)

ax.minorticks_on()
ax.set_xlim(0,)
ax.set_xlabel("Time (ka)", fontsize = 14)
ax.set_ylabel("Uplift rate (mm/yr)", fontsize = 14)
ax1.set_ylabel(r'$\Delta P \ (\%)$', fontsize = 14)
ax.tick_params(axis='both', which='major', labelsize=14)
ax1.tick_params(axis='both', which='major', labelsize=14)

legend_elements = [
    Line2D([0], [0], color='gray', ls = ':', lw=1, label=r'$U_{K-\Delta K}$'),
    Line2D([0], [0], color='b', lw=1, ls = '-', label=r'$U_{K}$'),
    Line2D([0], [0], color='k', lw=1, ls = '--', label=r'$U_{K+\Delta K}$'), 
    Line2D([0], [0], color='orange', lw=1, ls = '-', label=r'$\Delta P$')
]
plt.legend(handles=legend_elements, fontsize=10)

ax2 = axes[0]
ax2.set_xlim(ax.get_xlim()[0], ax.get_xlim()[1])
ax2.set_ylim(0, 9)

marker_style = dict(
    # fmt='o',
    markersize=5,
    # markerfacecolor='red',
    markeredgecolor='black',
    markeredgewidth=1,
    capsize=5, 
    ecolor='k')

ax2.axvspan(17, 29, color='grey', alpha=0.3) # Philip et al. (2012)
ax2.axvspan(2, 5.8, color='grey', alpha=0.3) # Philip et al. (2012)
# ax2.errorbar(7.55, 1, xerr=0.3, fmt='*') # Kaushal etal. (2023)
ax2.scatter(7.55, 1, marker = '*', facecolor = 'brown') # Kaushal etal. (2023)
ax2.errorbar(33, 2, xerr=2.3, fmt = 'o', mfc = 'b', **marker_style) # Suresh and Kumar (2020)
ax2.errorbar(26.1, 3, xerr=2.6, fmt = 's', mfc = 'r', **marker_style) # Suresh and Kumar (2020)
ax2.errorbar(26.5, 4, xerr=3.5, fmt = '^', mfc = 'g', **marker_style) # Kordt etal. (2025)
ax2.errorbar(24, 5, xerr=1, fmt = 'D', mfc = 'orange', **marker_style) # Vassallo etal. (2015)
ax2.errorbar(42.94, 6, xerr=3.57, fmt = 'v', mfc = 'purple', **marker_style) # Thakur etal. (2014)


ax2.set_xticks([])   # remove x ticks
ax2.set_yticks([])   # remove y ticks
ax2.set_xlabel("")   # remove x label
ax2.set_ylabel("")   # remove y label

plt.tight_layout(pad=0.3)

plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\K-sensitivity&climate&uplift_v2.svg')

# """


# %% Wild bootstrap. Multiplying residual with random mammen multiplier
opt_m = 0.45
opt_ke = 1.43e-4 # For khajnawar basin (mohand anticline)
node_mat = CreateChi(str_subset, facc_subset, fdir_subset, dem_subset,
                     X_coord_2d_subset, Y_coord_2d_subset, opt_m,
                     flow_method='arcgis')
import pandas as pd
df1 = pd.DataFrame({'X' : node_mat[:, 4], 
                   'Y' : node_mat[:, 5], 
                   'Elevation' : node_mat[:, 8] - node_mat[:, 8].min(), 
                   'Distance' : node_mat[:, 3], 
                   'tau' : node_mat[:, 7] / opt_ke})

tau = df1['tau']
z = df1['Elevation']

theta_hat, z_hat, z_t, tau_t = fit_once(tau, z)
T, U = theta_hat[:, 0], theta_hat[:, 1]
# residual
r = z_t - z_hat
B = 2000
U_hat = []
Z_hat = []
TAU_hat = []

rng = np.random.default_rng(0)

for b in range(B):
    v = mammen(len(r), rng)
    z_star = z_hat + v * r  # wild bootstrap
    try:
        theta_b, z_hat_b, _, tau_star_t = fit_once(tau_t, z_star)  # warm start at theta_hat
        U_hat.append(theta_b[:, 1])
        Z_hat.append(z_hat_b)
        TAU_hat.append(tau_star_t)
    except Exception:
        pass  # handle failures if any
        
    if b%100==0:
        print(b)

U_hat = np.array(U_hat)
U_med = np.median(U_hat, axis=0)
ci_lo = np.percentile(U_hat, 2.5, axis=0)
ci_hi = np.percentile(U_hat, 97.5, axis=0)

# ========= For plotting exercise of uplift ==============================================

plt.figure(figsize=(6, 2))
plt.plot(theta_hat[:, 0] / 1e3, theta_hat[:, 1]* 1e3, drawstyle = 'steps-post', 
         color = 'b', ls = '-', alpha = 1, 
         label = 'model uplift rate')

# --------- For rug/point cloud overlay plot -----
t_half = np.zeros((len(T) - 1, )) # Point cloud will be pltted at the mid-point of each time interval
for i in range(len(T) - 1):
    t_half[i] = (T[i] + T[i + 1]) / 2

plt.figure(figsize=(6, 2.5))
for i in range(len(U_hat)):
    plt.plot(t_half/1e3, U_hat[i, :-1]*1e3, 'k.', alpha=0.03, markersize=2)    

# plt.plot(T/1e3, U_med*1e3, 'r-', lw=2, label="Median uplift")
plt.step(T/1e3, U_med*1e3, 'r-', lw=2, where="post")

# Proxy artists for legend
legend_elements = [
    Line2D([0], [0], color='k', marker='.', linestyle='None', alpha=0.5, label="Bootstrap samples"),
    Line2D([0], [0], color='r', lw=2, label="Median uplift")
]

plt.xlabel("Time (ka)")
plt.ylabel("Uplift rate (mm/yr)")
plt.xlim(0,)
plt.legend(handles=legend_elements)
plt.tight_layout()

# ---------- Fan Chart (Quantile Bands)
q5, q25, q50, q75, q95 = np.percentile(U_hat, [5,25,50, 75,95], axis=0)
plt.figure(figsize=(6, 2.5))
plt.fill_between(T/1e3, q5*1e3, q95*1e3, step="post", color="skyblue", alpha=0.3, label="5–95% range")
plt.fill_between(T/1e3, q25*1e3, q75*1e3, step="post", color="dodgerblue", alpha=0.5, label="25–75% range")
plt.step(T/1e3, U_med*1e3, 'r-', lw=2, where="post", label="Median uplift")
plt.xlabel("Time (ka)")
plt.ylabel("Uplift rate (mm/yr)")
plt.xlim(0,)
plt.legend()
plt.tight_layout()

# ----------- Spaghetti Plot
plt.figure(figsize=(6, 2.5))
idx = np.random.choice(len(U_hat), size=100, replace=False)
for i in idx:  # plot a subset for clarity
    plt.step(T/1e3, U_hat[i]*1e3, color="gray", alpha=0.02, where="post")

plt.step(T/1e3, U_med*1e3, 'b-', lw=2, where="post")

# Proxy artists
legend_elements = [
    Line2D([0], [0], color="gray", lw=1, alpha=0.5, label="Bootstrap realizations"),
    Line2D([0], [0], color="b", lw=2, label="Median uplift")
]

plt.xlabel("Time (ka)")
plt.ylabel("Uplift rate (mm/yr)")
plt.xlim(0,)
plt.legend(handles=legend_elements)
plt.tight_layout()

# ------------ Spaghetti + Fan chart (quantile bands)
# plt.figure(figsize=(6, 2.5))
plt.figure(figsize=(9, 3))
idx = np.random.choice(len(U_hat), size=100, replace=False)
for i in idx:  # plot a subset for clarity
    plt.step(T/1e3, U_hat[i]*1e3, color="gray", alpha=0.04, where="post")

q5, q25, q75, q95 = np.percentile(U_hat, [5,25,75,95], axis=0)
plt.fill_between(T/1e3, q5*1e3, q95*1e3, step="post", color="skyblue", alpha=0.3, label="5–95% range")
plt.fill_between(T/1e3, q25*1e3, q75*1e3, step="post", color="dodgerblue", alpha=0.5, label="25–75% range")
plt.step(T/1e3, U_med*1e3, 'b-', lw=1, where="post", label="Median uplift")
plt.xlabel("Time (ka)", fontsize=15)
plt.ylabel("Uplift rate (mm/yr)", fontsize=15)
plt.xticks(fontsize=15)   # x-axis tick labels
plt.yticks(fontsize=15)
plt.xlim(0,)

# Proxy artists
legend_elements = [
    mpatches.Patch(color="skyblue", alpha=0.3, label="5–95% range"),
    mpatches.Patch(color="dodgerblue", alpha=0.5, label="25–75% range"),
    Line2D([0], [0], color="gray", lw=1, alpha=0.2, label="Bootstrap realizations"),
    Line2D([0], [0], color="b", lw=1, label="Median uplift")
]
plt.legend(handles=legend_elements, fontsize = 12, ncols=1)
plt.tight_layout()

# Uplift rate
plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\estimated_u_rate.svg')
plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\estimated_u_rate.png', dpi = 600)

# ========= For plotting exercise of fitted_Z vs observed Z ==============================================

# Fan Chart (quantile shading)
A0 = 1 # unit of m2
chi_t = tau_t * opt_ke * (A0 ** opt_m) # unit of meter
chi_star_t = tau_star_t * opt_ke * (A0 ** opt_m)
# plt.figure(figsize=(5,3.5))

plt.figure(figsize=(3,2.5))
idx = np.random.choice(len(Z_hat), size=100, replace=False)
for i in idx:  # plot a subset for clarity
    plt.plot(chi_star_t, Z_hat[i], color="dodgerblue", alpha=0.05)

z_q5, z_q25, z_q50, z_q75, z_q95 = np.percentile(np.array(Z_hat), [5,25,50,75,95], axis=0)
plt.scatter(chi_t, z_t, marker = '.', s = 10, c = 'k', 
            alpha = 0.1, label = 'Observation')
# plt.fill_between(tau_star_t/1e3, z_q5, z_q95, color="skyblue", alpha=0.3, label="5–95% range")
# plt.fill_between(tau_star_t/1e3, z_q25, z_q75, color="dodgerblue", alpha=0.5, label="25–75% range")
plt.plot(chi_star_t, z_q50, 'k-', lw=1, label="Median fit")

legend_elements = [
    # Line2D([0], [0], color='k', marker='.', linestyle='None', alpha=0.5, label="Observations"),
    Line2D([0], [0], color="b", lw=1, alpha=0.5, label="Bootstrap \n realizations"),
    Line2D([0], [0], color="k", lw=1, label="Median fit")
]

plt.ylabel('Elevation (m)', fontsize = 12)
plt.xlabel(r'$\chi$ (m)', fontsize = 12)
plt.legend(fontsize = 10, handles=legend_elements)
plt.xlim(0,)
plt.ylim(0,)
plt.xticks(fontsize=12)   # x-axis tick labels
plt.yticks(fontsize=12)
import matplotlib.ticker as ticker
plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
plt.tight_layout()

# ------------- Save the figures ------------------------
# obs Z and fitted Z
plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\Zfit.svg')
plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\Zfit.png', dpi = 600)

# ---- If cumulative uplift ----
# For cumulative uplift (baselevel fall)
def cum_U(T, U):
    UT = np.vstack((T, U)).T
    UT_sorted = UT[UT[:, 0].argsort()[::-1]]
    U_cum = np.zeros((U.shape))
    for i in range(len(U) - 1):
        aoc = (UT_sorted[i, 0] - UT_sorted[i + 1, 0]) * (UT_sorted[i, 1])
        U_cum[i + 1] = U_cum[i] + aoc
    return U_cum

U_cum_hat = []
for i in range(U_hat.shape[0]):
    U_cum = cum_U(T, U_hat[i, :])
    U_cum_hat.append(U_cum)

U_cum_hat = np.array(U_cum_hat)
mean_cum = np.mean(U_cum_hat, axis=0)
lower_cum = np.percentile(U_cum_hat, 2.5, axis=0)
upper_cum = np.percentile(U_cum_hat, 97.5, axis=0)

U_cum_med = cum_U(T, q50)
U_cum_li = cum_U(T, np.percentile(U_hat, 2.5, axis=0))
U_cum_hi = cum_U(T, np.percentile(U_hat, 97.5, axis=0))

plt.figure(figsize = (9, 3))
plt.fill_between(np.flipud(T / 1e3), U_cum_li, U_cum_hi, step="pre", color="skyblue", 
                 alpha=0.3, label=r"$2\sigma \ bound$")
# plt.fill_between(T/1e3, q25*1e3, q75*1e3, step="post", color="dodgerblue", alpha=0.5, label="25–75% range")
plt.step(np.flipud(T / 1e3), U_cum_med, 'b-', lw=1, where="pre", label="Cumulative uplift")
plt.xlabel("Time (ka)", fontsize=14)
plt.ylabel("Total relief (m)", fontsize=14)
plt.xticks(fontsize=14)   # x-axis tick labels
plt.yticks(fontsize=14)
plt.xlim(0,)
plt.legend(fontsize = 10, ncols=1)
plt.tight_layout()

plt.savefig(r'D:\ramendra_EPM103\Rajeeb\figures\cum_uplift.svg')

# # Proxy artists
# legend_elements = [
#     mpatches.Patch(color="skyblue", alpha=0.3, label="5–95% range"),
#     mpatches.Patch(color="dodgerblue", alpha=0.5, label="25–75% range"),
#     Line2D([0], [0], color="gray", lw=1, alpha=0.2, label="Bootstrap realizations"),
#     Line2D([0], [0], color="b", lw=1, label="Median uplift")
# ]
# plt.legend(handles=legend_elements, fontsize = 12, ncols=1)
# plt.tight_layout()

# UT = np.vstack((T, U_med, q5, q95)).T
# UT_sorted = UT[UT[:, 0].argsort()[::-1]]
# U_cum = np.zeros((U_med.shape))
# U5_cum = np.zeros((U_med.shape))
# U50_cum = np.zeros((U_med.shape))
# U50_cum = np.zeros((U_med.shape))
# for i in range(len(U_med) - 1):
#     # aoc = 0.5 * (UT_sorted[i, 0] - UT_sorted[i + 1, 0]) * (UT_sorted[i, 1] + UT_sorted[i + 1, 1])
#     aoc = (UT_sorted[i, 0] - UT_sorted[i + 1, 0]) * (UT_sorted[i, 1])
#     U_cum[i + 1] = U_cum[i] + aoc

# # from scipy.integrate import cumtrapz
# # U_cum = cumtrapz(U, x=T[0: q], initial=0)
# plt.figure()
# plt.plot(np.flipud(T / 1e3), U_cum, drawstyle = 'steps-pre')
    

# %% Another fun thing to try out..bayesian inversion with basin NUT MCMC sampling using pymc
import pymc as pm
import arviz as az
n_bins = 10 # we have 10 time intervals
z_obs = z[0:A.shape[0]]

# --------------------------------------------------
# PyMC Bayesian inversion model. Here theta is our U vector.
# --------------------------------------------------
# Taking a prior in the form of Gaussian random walk (correlated thetas with predefined sigmas)
np.random.seed(42)
with pm.Model() as model:
    # Smoothness hyperparameter for random walk 
    sigma_theta = pm.HalfNormal("sigma_theta", sigma=8.0) #d_theta = 2e-3 m/yr (2 mm/yr)

    # Unconstrained Gaussian Random Walk
    theta_raw = pm.GaussianRandomWalk("theta_raw", sigma=sigma_theta, shape=n_bins)

    # Transform uplift rates to [0,10] mm/yr
    theta = pm.Deterministic("theta", 10e-3 * pm.math.sigmoid(theta_raw))

    prior_samples = pm.sample_prior_predictive(samples=2000)

    # Forward model: predicted elevations
    mu = pm.math.dot(A, theta)

    # Observation error
    sigma_obs = pm.HalfNormal("sigma_obs", sigma=5.0)

    # Likelihood: observed elevations
    y_obs = pm.Normal("y_obs", mu=mu, sigma=sigma_obs, observed=z_obs)

    # --------------------------------------------------
    # Run MCMC
    # --------------------------------------------------
    trace = pm.sample(draws=2000, tune=2000, target_accept=0.9, chains=4, random_seed=42)
    
    # --------------------------------------------------
    # Posterior Predictive Checks
    # --------------------------------------------------
    ppc = pm.sample_posterior_predictive(trace, var_names=["y_obs"], random_seed=42)

# Enforce piecewise-constant thetas without any prior correlation
np.random.seed(42)
with pm.Model() as piecewise_model:
    # independent (non-smooth) prior per bin:
    theta = pm.TruncatedNormal("theta", mu=5e-3, sigma=3e-3, lower=0.0, upper=1e-2, shape=n_bins)
    # Forward model: predicted elevations
    mu = pm.math.dot(A, theta)

    sigma_obs = pm.HalfNormal("sigma_obs", sigma=5.0)
    y_obs = pm.Normal("y_obs", mu=mu, sigma=sigma_obs, observed=z_obs)

    trace_pw = pm.sample(2000, tune=2000, target_accept=0.9, chains=4, random_seed=42)
    
    # Posterior Predictive Checks
    # --------------------------------------------------
    ppc = pm.sample_posterior_predictive(trace_pw, var_names=["y_obs"], random_seed=42)
# --------------------------------------------------
# Summaries
# --------------------------------------------------
trace = trace_pw
print(az.summary(trace, var_names=["theta", "sigma_theta", "sigma_obs"]))

# Trace and posterior plots
az.plot_trace(trace, var_names=["theta"])

az.plot_posterior(trace, var_names=["theta"])
plt.show()

# Extract posterior uplift rates (theta) ----------------------
theta_samples = trace.posterior["theta"].values
theta_samples = theta_samples.reshape(-1, theta_samples.shape[-1]) # (n_samples, n_bins)
last_col = theta_samples[:, [9]]
theta_extended = np.hstack([theta_samples, last_col])
# Compute summaries
median = np.median(theta_extended, axis=0)
# hdi = az.hdi(theta_extended, hdi_prob=0.95)

q5, q25, q50, q75, q95 = np.percentile(theta_extended, [5,25,50,75,95], axis=0)
plt.figure(figsize=(6, 2.5))
plt.fill_between(T/1e3, q5*1e3, q95*1e3, step="post", color="skyblue", alpha=0.3, label="5–95% range")
plt.fill_between(T/1e3, q25*1e3, q75*1e3, step="post", color="dodgerblue", alpha=0.5, label="25–75% range")
plt.step(T/1e3, q50*1e3, 'r-', lw=2, where="post", label="Median uplift")
plt.xlabel("Time (ka)")
plt.ylabel("Uplift rate (mm/yr)")
plt.xlim(0,)
plt.legend()
plt.tight_layout()


# --------------------------------------------------
# Plot Posterior Predictive vs Observed
# --------------------------------------------------
ppc_y = ppc.posterior_predictive["y_obs"].values   # ndarray
ppc_y = ppc_y.reshape(-1, ppc_y.shape[-1])        # (n_samples, n_obs)

plt.figure(figsize=(5,4))

# Plot observed river profile
plt.plot(z_obs, "ko", label="Observed elevations", markersize=2)

# Plot posterior predictive draws (spaghetti plot, thin alpha)
for i in range(100):
    plt.plot(ppc_y[i], color="skyblue", alpha=0.04)

# Mean posterior predictive
plt.plot(np.median(ppc_y, axis=0), color="blue", lw=2, label="Posterior mean fit")

plt.xlabel("Observation index (along profile)")
plt.ylabel("Elevation (m)")
plt.title("Posterior Predictive Check: River Profile Fit")
plt.legend()
plt.xlim(0,)
plt.tight_layout()

# ----- plot prior vs posterior for the uplift rate ---------------
theta_prior = prior_samples.prior["theta"]  # shape (n_prior_draws, 10)
theta_post  = trace.posterior["theta"].stack(sample=("chain","draw")).values.T  # shape (10, n_post_draws)

fig, axes = plt.subplots(2, 5, figsize=(15,6), sharey=True)
axes = axes.ravel()

import seaborn as sns
for i in range(10):
    sns.kdeplot(theta_prior[0, :, i], ax=axes[i], label="Prior", color="blue")
    sns.kdeplot(theta_post[:, i], ax=axes[i], label="Posterior", color="red")
    axes[i].set_title(f"Bin {i+1}")
    if i == 0:
        axes[i].legend()

plt.tight_layout()


# %%
plt.savefig('/home/cds/Desktop/ramendra/saikat/results/zVSzinf_southAnd.pdf')
plt.savefig('/home/cds/Desktop/ramendra/saikat/results/uplift_rate_henry.pdf')
plt.savefig('/home/cds/Desktop/ramendra/saikat/results/Cum_U_southAnd.pdf')
plt.savefig('/home/cds/Desktop/ramendra/saikat/results/Cum_U_v2.pdf')



