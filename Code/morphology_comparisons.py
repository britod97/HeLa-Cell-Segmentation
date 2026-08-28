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
from scipy import stats
from sklearn.covariance import EllipticEnvelope
from sklearn.preprocessing import StandardScaler


def compute_cell_summary(mito_df, nucleus_df):
    """One row per ROI: nucleus features + aggregated mitochondria population features."""
    mito_agg = mito_df.groupby('ROI').agg(
        mito_count=('volume', 'count'),
        mito_volume_mean=('volume', 'mean'),
        mito_volume_total=('volume', 'sum'),
        mito_surface_area_mean=('surface_area', 'mean'),
        mito_centroid_distance_mean=('centroid_distance', 'mean'),
        mito_surface_distance_mean=('surface_distance', 'mean'),
        mito_centroid_angle_mean=('centroid_angle', 'mean'),
        mito_surface_angle_mean=('surface_angle', 'mean'),
    ).reset_index()

    nucleus_cols = ['ROI', 'volume', 'surface_area', 'anisotropy',
                     'mitocloud_L1', 'mitocloud_L2', 'mitocloud_L3',
                     'distance_of_centroids']
    nucleus_cols = [c for c in nucleus_cols if c in nucleus_df.columns]
    nuc = nucleus_df[nucleus_cols].rename(
        columns={c: f'nucleus_{c}' for c in nucleus_cols if c != 'ROI'}
    )

    cell_summary = nuc.merge(mito_agg, on='ROI', how='left')
    cell_summary['mito_volume_fraction'] = (
        cell_summary['mito_volume_total'] / cell_summary['nucleus_volume']
    )
    return cell_summary

def flag_feature_outliers(new_summary, old_summary, feature_cols=None, z_thresh=3, robust=True):
    if feature_cols is None:
        feature_cols = [c for c in new_summary.columns
                         if c != 'ROI' and pd.api.types.is_numeric_dtype(new_summary[c])
                         and c in old_summary.columns]

    results = new_summary[['ROI']].copy()
    flags = pd.DataFrame(index=new_summary.index)

    for col in feature_cols:
        old_vals = old_summary[col].dropna()
        if robust:
            center = old_vals.median()
            spread = stats.median_abs_deviation(old_vals, scale='normal')
        else:
            center = old_vals.mean()
            spread = old_vals.std()
        if spread == 0 or np.isnan(spread):
            continue
        z = (new_summary[col] - center) / spread
        results[f'{col}_z'] = z
        flags[col] = z.abs() > z_thresh

    results['n_flagged_features'] = flags.sum(axis=1)
    results['flagged_features'] = flags.apply(lambda r: [c for c in flags.columns if r[c]], axis=1)
    results['is_outlier_cell'] = results['n_flagged_features'] > 0
    return results, flags

def flag_multivariate_outliers(new_summary, old_summary, feature_cols, contamination=0.05):
    old_clean = old_summary[feature_cols].dropna()
    scaler = StandardScaler().fit(old_clean)
    old_scaled = scaler.transform(old_clean)

    detector = EllipticEnvelope(contamination=contamination, random_state=0).fit(old_scaled)

    new_clean = new_summary[feature_cols].dropna()
    new_scaled = scaler.transform(new_clean)
    scores = detector.score_samples(new_scaled)  # higher = more "normal"
    preds = detector.predict(new_scaled)          # -1 = outlier

    out = new_summary.loc[new_clean.index, ['ROI']].copy()
    out['mahalanobis_score'] = -scores
    out['is_multivariate_outlier'] = preds == -1
    return out


if __name__ == "__main__":
    # New Data
    new_data_dir = Path(r'E:\HeLa\Data\CIL50051\GeneratedData')
    new_mitochondria_dir = Path(f'{new_data_dir}/mitochondria_props_isotropic')
    new_nucleus_data_path = Path(f'{new_data_dir}/nucleus_props_isotropic/nucleus_properties.csv')
    
    new_mito_data, new_nucleus_data = hela_utils.load_mito_data(new_mitochondria_dir, new_nucleus_data_path)
    
    
    # Old data
    old_data_dir = Path(r'E:\HeLa\Data\EMPIAR-10094\GeneratedData')
    old_mitochondria_dir = Path(f'{old_data_dir}/mitochondria_props_isotropic')
    old_nucleus_data_path = Path(f'{old_data_dir}/nucleus_props_isotropic/nucleus_properties.csv')
    
    old_mito_data, old_nucleus_data = hela_utils.load_mito_data(old_mitochondria_dir, old_nucleus_data_path)

    new_cell_summary = compute_cell_summary(new_mito_data, new_nucleus_data)
    old_cell_summary = compute_cell_summary(old_mito_data, old_nucleus_data)
    
    
    
    feature_cols = [c for c in new_cell_summary.columns
                     if c != 'ROI' and pd.api.types.is_numeric_dtype(new_cell_summary[c])]
    
    univariate, flags = flag_feature_outliers(new_cell_summary, old_cell_summary, feature_cols)
    multivariate = flag_multivariate_outliers(new_cell_summary, old_cell_summary,
                                                feature_cols=['mito_count', 'mito_volume_mean',
                                                               'mito_volume_fraction', 'nucleus_volume',
                                                               'nucleus_anisotropy'])
    
    summary = univariate.merge(multivariate, on='ROI', how='left')
    flagged = summary[summary['is_outlier_cell'] | summary['is_multivariate_outlier']]
    print(flagged[['ROI', 'n_flagged_features', 'flagged_features', 'mahalanobis_score']])