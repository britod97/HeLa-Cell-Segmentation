# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 14:28:06 2026

@author: brito
"""

import os
from pathlib import Path
workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import numpy as np
import dask.array as da
import hela_utils
import pandas as pd



SPACING = np.array([3.6, 3.6, 60])
chunk_size = (64,64,255)

ZARR_DIR = Path(r"E:\HeLa\Data\CIL50051\Zarr")
MITO_DIR = Path(r"E:\HeLa\Data\CIL50051\GeneratedData\mitochondria_props_isotropic")
NUCLEUS_DIR = Path(r"E:\HeLa\Data\CIL50051\GeneratedData\nucleus_props_isotropic")
os.makedirs(MITO_DIR, exist_ok=True)
os.makedirs(NUCLEUS_DIR, exist_ok=True)

ROIs = np.unique(['_'.join(ROI.split('.')[0].split('_')[:-1])
        for ROI in os.listdir(ZARR_DIR)
        if ROI.endswith('.zarr')])


###############################################################################
nucleus_results = []
k=0
for ROI_name in ROIs:
    print(f'{k+1}/{len(ROIs)}')
    k+=1
    print(f'ROI: {ROI_name}')
    nucleus_lazy = da.from_zarr(f"{ZARR_DIR}/{ROI_name}_Nuclei.zarr")
    bbox = hela_utils.dask_bbox(nucleus_lazy, pad=2)
    slices = tuple(slice(lo, hi) for lo, hi in bbox)
    offset = np.array([s.start for s in slices])

    nucleus = nucleus_lazy[slices].compute().astype('uint8')

    mito_labels = hela_utils.compute_labels(ZARR_DIR, ROI_name,
                                             chunk_size=chunk_size,
                                             SPACING=SPACING)
    nucleus_props = hela_utils.compute_nucleus_properties(nucleus, ROI_name,
                                                            SPACING=SPACING,
                                                            offset=offset)

    mito_df = hela_utils.compute_mitochondria_properties(mito_labels,
                                                         SPACING=SPACING
                                                         )

    mito_df.to_csv(f"{MITO_DIR}/{ROI_name}_properties.csv", index=False)

    nucleus_results.append(nucleus_props)


nucleus_df = pd.DataFrame(nucleus_results)
nucleus_df.to_csv(f"{NUCLEUS_DIR}/nucleus_properties.csv", index=False)

###############################################################################
nucleus_df = pd.read_csv(f"{NUCLEUS_DIR}/nucleus_properties.csv")

for ROI_name in ROIs:
    print(ROI_name)
    df = pd.read_csv(f"{MITO_DIR}/{ROI_name}_properties.csv")

    # Append all features
    df = hela_utils.append_angle_relative_to_nucleus_centroid(df, nucleus_df, ROI_name)
    # print("Centroid angle done.")
    df = hela_utils.append_angle_relative_to_nucleus_surface(df, nucleus_df, ROI_name,
                                                             ZARR_DIR,
                                                             chunk_size=chunk_size,
                                                             SPACING=SPACING)
    # print("Surface angle done.")
    df = hela_utils.append_spherical_coordinates(df, nucleus_df, ROI_name)
    print("Spherical coordinates.")
    
    df.to_csv(f"{MITO_DIR}/{ROI_name}_properties.csv", index=False)
    
###############################################################################
nucleus_df = pd.read_csv(f"{NUCLEUS_DIR}/nucleus_properties.csv")

for idx, row in nucleus_df.iterrows():
    df = pd.read_csv(f"{MITO_DIR}/{row['ROI']}_properties.csv")
    dic = hela_utils.compute_mito_cloud_properties(df, row)
    for key, value in dic.items():
        nucleus_df.at[idx, key] = value

nucleus_df.to_csv(f"{NUCLEUS_DIR}/nucleus_properties.csv", index=False)

# ###############################################################################
for ROI_name in ROIs:
    df = pd.read_csv(f"{MITO_DIR}/{ROI_name}_properties.csv")
    
    df = hela_utils.append_PCA_distribution(df, nucleus_df, ROI_name)
    print("Global direction")
    df = hela_utils.append_angle_relative_to_global_directions(df, nucleus_df, ROI_name)
    
    df.to_csv(f"{MITO_DIR}/{ROI_name}_properties.csv", index=False)