# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 14:48:02 2026

@author: adhh334
"""

import os
from pathlib import Path
workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
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
import hela_utils



# Morphology comparisons

ZARR_DIR = Path(r"E:\HeLa\Data\CIL50051\Zarr")
# New Data
data_dir = Path(r'E:\HeLa\Data\CIL50051\GeneratedData')
mitochondria_dir = Path(f'{data_dir}/mitochondria_props_isotropic')

ROIs = np.unique(['_'.join(ROI.split('.')[0].split('_')[:-1])
        for ROI in os.listdir(mitochondria_dir)
        if ROI.endswith('.csv')])

# nucleus_data = pd.read_csv(f'{data_dir}/nucleus_props_isotropic/nucleus_properties.csv')
mito_data_list = []
for ROI_name in ROIs:
    mito_data = pd.read_csv(f'{mitochondria_dir}/{ROI_name}_properties.csv')
    mito_data['ROI'] = ROI_name
    mito_data_list.append(mito_data)

# all_mito_data = pd.concat(mito_data_list)

# # Append nucleus properties to all_mito_data
# all_mito_data = all_mito_data.merge(nucleus_data, suffixes=("", "_nucleus") ,on='ROI')


# all_mito_data = all_mito_data.drop(['distance_to_current_mito'], axis=1, errors='ignore')