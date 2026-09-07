# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 13:30:06 2026

@author: adhh334
"""

import os
from pathlib import Path
workdir = Path(r"D:\GitHub Repos\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from skimage.measure import marching_cubes
from mpl_toolkits.mplot3d import Axes3D
import dask.array as da
from tqdm import tqdm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import colors
import scipy.io
from matplotlib.markers import MarkerStyle
from matplotlib.ticker import FormatStrFormatter

def export_legend(legend, filename="legend.pdf", expand=[-5,-5,5,5]):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)



ZARR_DIR = Path(r"D:\GitHub Repos\HeLa_Cell_Data\CIL50051\Zarr")

# Read Nucleus and Mitochondria data
data_dir = Path(r'D:\GitHub Repos\HeLa_Cell_Data\CIL50051\GeneratedData')
mitochondria_dir = Path(f'{data_dir}/mitochondria_props_isotropic')

ROIs = np.unique(['_'.join(ROI.split('.')[0].split('_')[:-1])
        for ROI in os.listdir(ZARR_DIR)
        if ROI.endswith('.zarr')])

nucleus_data = pd.read_csv(f'{data_dir}/nucleus_props_isotropic/nucleus_properties.csv')
mito_data_list = []
for ROI_name in ROIs:
    mito_data = pd.read_csv(f'{mitochondria_dir}/{ROI_name}_properties.csv')
    mito_data['ROI'] = ROI_name
    mito_data_list.append(mito_data)

all_mito_data = pd.concat(mito_data_list)

# Append nucleus properties to all_mito_data
all_mito_data = all_mito_data.merge(nucleus_data, suffixes=("", "_nucleus") ,on='ROI')


all_mito_data = all_mito_data.drop(['distance_to_current_mito'], axis=1, errors='ignore')


###############################################################################
# Orientation of mitochondria with respect to the CENTROID by THRESHOLDS

n_bins = 9
normalised_heights = True
bar_width = 1


distance_types = ['centroid_distance', 'surface_distance']
angle_types = ['centroid_angle', 'surface_angle']

for d_type in distance_types:
    for angle in angle_types:
        distances = all_mito_data['centroid_distance']
        figure_dir = Path(f'D:\GitHub Repos\HeLa_Cell_Data\CIL50051\Figures/ROSE_isotropic/{angle}/by_{d_type}')
        os.makedirs(figure_dir, exist_ok=True)
    
        P = [25, 50, 100]
        percentiles = np.percentile(distances, P)
        k=0
        for p in percentiles:
            fig = plt.figure(figsize=(10,10))
            mito_data = all_mito_data[distances <= p]
            hist = np.histogram(mito_data[angle], bins=n_bins, range=(0, np.pi))
    
            radii = hist[0]
            if normalised_heights:
                radii = radii/np.sum(radii)
    
            theta = [(hist[1][k]+hist[1][k+1])/2 for k in range(len(hist[1])-1)]
            width = bar_width*(np.pi/n_bins)
    
            # Coloring
            norm = colors.Normalize(vmin=0, vmax=0.4)
            cmap = plt.cm.jet
            C = cmap(norm(radii))
    
            # Plotting
            ax = plt.subplot(projection='polar')
            chart = ax.bar(theta, radii, width=width, bottom=0.0,
               color=C, alpha=1, edgecolor='k', zorder=10, lw=3)
    
    
            ax.set_title(f'{angle} of mitos with {d_type} <= {round(p,2)}')
            # ax.set_yticks([0, 0.2, 0.4])
            # ax.set_ylim(0, 0.4)
            ax.set_thetalim(0, np.pi)
            ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
            ax.tick_params(axis='y', labelsize=36, pad=10)
            ax.tick_params(axis='x', labelsize=42, pad=30)
    
            ax.grid('both', lw=3, zorder=-10, alpha=0.8)
            for spine in ax.spines.values():
                spine.set_linewidth(3)
            plt.savefig(f'{figure_dir}/{d_type}_{P[k]}.pdf')
            plt.savefig(f'{figure_dir}/{d_type}_{P[k]}.png', dpi=300, bbox_inches='tight', transparent=True)
    
            # --- Separate colorbar export ---
            fig_cb, ax_cb = plt.subplots(figsize=(1, 6))  # (width, height)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
    
            cbar = plt.colorbar(sm, cax=ax_cb)
            cbar.set_label('Normalised Frequency', rotation=270, labelpad=15)
    
            # Save only the colorbar
            # fig_cb.savefig(f'{figure_dir}/colorbar.pdf', bbox_inches='tight')
            # fig_cb.savefig(f'{figure_dir}/colorbar.png', dpi=1200, bbox_inches='tight', transparent=True)
            plt.close(fig_cb)
    
            k += 1


###############################################################################
# Orientation of mitochondria with respect to the CENTROID by CELL

n_bins = 9
normalised_heights = True
bar_width = 1

angle_types = ['centroid_angle', 'surface_angle']


for angle in angle_types:
    figure_dir = Path(f'D:/GitHub Repos/HeLa_Cell_Data/CIL50051/Figures/ROSE_isotropic/{angle}/by_cell')
    os.makedirs(figure_dir, exist_ok=True)
    for roi_name in ROIs:
        fig = plt.figure(figsize=(10,10))
        mito_data = all_mito_data.loc[all_mito_data['ROI'] == roi_name]
        hist = np.histogram(mito_data[angle], bins=n_bins, range=(0, np.pi))
    
        radii = hist[0]
        if normalised_heights:
            radii = radii/np.sum(radii)
    
        theta = [(hist[1][k]+hist[1][k+1])/2 for k in range(len(hist[1])-1)]
        width = bar_width*(np.pi/n_bins)
    
        # Coloring
        norm = colors.Normalize(vmin=0, vmax=0.6)
        cmap = plt.cm.jet
        C = cmap(norm(radii))
    
        # Plotting
        ax = plt.subplot(projection='polar')
        chart = ax.bar(theta, radii, width=width, bottom=0.0,
               color=C, alpha=1, edgecolor='k', lw=3, zorder=10)
    
        ax.set_title(f'{roi_name} {angle}')
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
        ax.set_ylim(0, 0.4)
        ax.set_thetalim(0, np.pi)
        ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
        ax.tick_params(axis='y', labelsize=36, pad=10)
        ax.tick_params(axis='x', labelsize=42, pad=30)
        ax.grid('both', lw=2, zorder=-10, alpha=0.8)
        for spine in ax.spines.values():
            spine.set_linewidth(3)
    
        # --- Separate colorbar export ---
        fig_cb, ax_cb = plt.subplots(figsize=(1, 6))  # (width, height)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
    
        cbar = plt.colorbar(sm, cax=ax_cb)
        cbar.set_label('Normalised Frequency', rotation=270, labelpad=15)
    
        # Save only the colorbar
        fig_cb.savefig(f'{figure_dir}/colorbar.pdf', dpi=300, bbox_inches='tight', transparent=True)
        plt.close(fig_cb)
    
        plt.savefig(f'{figure_dir}/{roi_name[:11]}_{angle}.pdf')


