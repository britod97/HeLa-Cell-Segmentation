# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 14:48:02 2026

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
import hela_utils
from scipy import stats
from sklearn.covariance import EllipticEnvelope
from sklearn.preprocessing import StandardScaler
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)  # don't wrap based on terminal width


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
    """
    Flag cells in a new dataset whose per-feature values deviate strongly from
    the distribution of the same features in an old (reference) dataset.

    For each feature, a z-like score is computed:

        z = (new_value - center) / spread

    where `center` and `spread` are estimated from `old_summary` only:
      - robust=True  (default): center = median, spread = MAD scaled to be
        comparable to a standard deviation (`scipy.stats.median_abs_deviation`,
        `scale='normal'`). Less sensitive to outliers/skew in the old
        reference population than mean/std.
      - robust=False: center = mean, spread = std (classic z-score).

    A cell is flagged on a given feature if |z| > z_thresh. A cell is flagged
    overall (`is_outlier_cell`) if it exceeds the threshold on at least one
    feature. Features are evaluated independently (univariate) - this does
    NOT account for correlations between features. Use
    `flag_multivariate_outliers` for that.

    Features where the old-data spread is zero or NaN (e.g. too few ROIs, or
    a constant column) are silently skipped and excluded from flagging, since
    z-scores would be undefined or infinite.

    Parameters
    ----------
    new_summary : pd.DataFrame
    old_summary : pd.DataFrame
    feature_cols : list of str, optional
    z_thresh : float, default 3
        Absolute z-score threshold above which a feature is flagged.
    robust : bool, default True
        If True, use median/MAD (scaled) for center/spread. If False, use
        mean/std.

    Returns
    -------
    results : pd.DataFrame
        One row per ROI in `new_summary`, with columns:
          - 'ROI'
          - '{col}_z' for each tested feature: the signed z-score
          - 'n_flagged_features': count of features with |z| > z_thresh
          - 'flagged_features': list of feature names that were flagged
          - 'is_outlier_cell': True if n_flagged_features > 0
    flags : pd.DataFrame
        Boolean matrix (same index as `new_summary`), one column per tested
        feature, True where that feature was flagged for that cell. Columns
        with zero/NaN spread in `old_summary` are absent (not just False).
    """
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
    z_thresh = 2
    
    # New Data
    new_data_dir = Path(r'D:\GitHub Repos\HeLa_Cell_Data\CIL50051\GeneratedData')
    new_mitochondria_dir = Path(f'{new_data_dir}/mitochondria_props_isotropic')
    new_nucleus_data_path = Path(f'{new_data_dir}/nucleus_props_isotropic/nucleus_properties.csv')
    
    new_mito_data, new_nucleus_data = hela_utils.load_mito_data(new_mitochondria_dir, new_nucleus_data_path)
    
    
    # Old data
    old_data_dir = Path(r'D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\GeneratedData')
    old_mitochondria_dir = Path(f'{old_data_dir}/mitochondria_props_isotropic')
    old_nucleus_data_path = Path(f'{old_data_dir}/nucleus_props_isotropic/nucleus_properties.csv')
    
    old_mito_data, old_nucleus_data = hela_utils.load_mito_data(old_mitochondria_dir, old_nucleus_data_path)

    new_cell_summary = compute_cell_summary(new_mito_data, new_nucleus_data)
    old_cell_summary = compute_cell_summary(old_mito_data, old_nucleus_data)
    
    
    
    feature_cols = [c for c in new_cell_summary.columns
                     if c != 'ROI' and pd.api.types.is_numeric_dtype(new_cell_summary[c])]
    
    univariate, flags = flag_feature_outliers(new_cell_summary,
                                              old_cell_summary,
                                              feature_cols,
                                              z_thresh=z_thresh)
    multivariate = flag_multivariate_outliers(new_cell_summary,
                                              old_cell_summary,
                                                feature_cols=['mito_count', 'mito_volume_mean',
                                                               'mito_volume_fraction', 'nucleus_volume',
                                                               'nucleus_anisotropy'])
    
    summary = univariate.merge(multivariate, on='ROI', how='left')
    flagged = summary[summary['is_outlier_cell'] | summary['is_multivariate_outlier']]
    print(flagged[['ROI', 'n_flagged_features', 'flagged_features', 'mahalanobis_score']])