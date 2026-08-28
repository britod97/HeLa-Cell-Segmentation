# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 14:21:22 2026

@author: brito
"""

import pandas as pd
import dask.array as da
import numpy as np
from skimage.measure import regionprops
from skimage.measure import marching_cubes
import os
import sys
from sklearn.decomposition import PCA
import pickle
from scipy import ndimage
from tqdm import tqdm
from scipy.ndimage import distance_transform_edt

def dask_bbox(dask_arr, pad=2):
    """Find bounding box of nonzero voxels without materializing the full array."""
    nz = dask_arr > 0
    axes = list(range(dask_arr.ndim))
    bbox = []
    for ax in axes:
        other_axes = tuple(a for a in axes if a != ax)
        proj = nz.any(axis=other_axes).compute()  # small 1D result
        idx = np.nonzero(proj)[0]
        lo, hi = int(idx.min()), int(idx.max()) + 1
        bbox.append((max(0, lo - pad), hi + pad))
    return bbox


def load_mito_data(mitochondria_dir, nucleus_data_path):
    nucleus_data = pd.read_csv(nucleus_data_path)
    
    ROIs = np.unique(['_'.join(ROI.split('.')[0].split('_')[:-1])
            for ROI in os.listdir(mitochondria_dir)
            if ROI.endswith('.csv')])
    
    mito_data_list = []
    for ROI_name in ROIs:
        mito_data = pd.read_csv(f'{mitochondria_dir}/{ROI_name}_properties.csv')
        mito_data['ROI'] = ROI_name
        mito_data_list.append(mito_data)

    all_mito_data = pd.concat(mito_data_list)

    # # Append nucleus properties to all_mito_data
    all_mito_data = all_mito_data.merge(nucleus_data, suffixes=("", "_nucleus") ,on='ROI')
    all_mito_data = all_mito_data.drop(['distance_to_current_mito'], axis=1, errors='ignore')
    
    return all_mito_data, nucleus_data


################################################################
'''''''''''''''''
Basic Geometric Functions
'''''''''''''''''
def centroid_from_binary(mask, SPACING = np.array([10,10,50])):
    # voxel indices
    coords = np.nonzero(mask)   # tuple of (z_idx, y_idx, x_idx)
    # compute mean per axis
    centroid_vox = [np.mean(c) for c in coords]
    # convert to physical units
    centroid_phys = np.array(centroid_vox)*np.array(SPACING)
    return centroid_phys

def compute_angle_between_vectors(v1, v2):
    # Ensure v1 and v2 are 2D arrays for batch processing
    # If they are 1D (single vectors), convert them to 2D for consistent logic
    if v1.ndim == 1:
        v1 = v1[np.newaxis, :]
    if v2.ndim == 1:
        v2 = v2[np.newaxis, :]

    # Calculate row-wise dot product
    dot_product = np.sum(v1 * v2, axis=1)

    norm_v1 = np.linalg.norm(v1, axis=1)
    norm_v2 = np.linalg.norm(v2, axis=1)

    # Compute cosine of the angle, handling division by zero and clipping for arccos domain
    cosine_angle = dot_product / (norm_v1 * norm_v2 + 1e-12)
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0) # Ensure values are within [-1, 1] for arccos

    angle = np.arccos(cosine_angle)
    return angle.squeeze()

################################################################
'''''''''''''''''
Nuclear Properties
'''''''''''''''''
def build_nucleus_distance_map(ROI_name, ZARR_DIR, SPACING = np.array([10,10,50])):

    nucleus_path = f"{ZARR_DIR}/{ROI_name}_Nuclei.zarr"

    nucleus = da.from_zarr(nucleus_path).compute()

    dist = distance_transform_edt(
        nucleus == 0,
        sampling=SPACING
    )

    return nucleus, dist

def compute_nucleus_properties(nucleus_array, ROI_name,
                               SPACING = np.array([10,10,50]),
                               offset=np.array([0, 0, 0])):
    
    voxel_volume = np.prod(SPACING)
    prop = regionprops(nucleus_array)[0]
    centroid = centroid_from_binary(nucleus_array, SPACING=SPACING)
    centroid = centroid + offset * SPACING  # shift back to global coords

    volume = prop.area * voxel_volume

    verts, faces, _, _ = marching_cubes(
        nucleus_array,
        level=0.5,
        spacing=SPACING
    )
    verts = verts + offset * SPACING  # shift mesh vertices too

    tri_verts = verts[faces]
    vec1 = tri_verts[:,1] - tri_verts[:,0]
    vec2 = tri_verts[:,2] - tri_verts[:,0]
    cross = np.cross(vec1, vec2)
    surface_area = np.linalg.norm(cross, axis=1).sum() / 2

    return {
        "ROI": ROI_name,
        "centroid_x": centroid[0],
        "centroid_y": centroid[1],
        "centroid_z": centroid[2],
        "volume": volume,
        "surface_area": surface_area
    }

################################################################
'''''''''''''''''
Mitochondria Properties
'''''''''''''''''

# Collection 2
def compute_labels(zarr_dir, ROI_name, chunk_size=(64,64,300), SPACING = np.array([10,10,50])):
    mitochondria_path = (f'{zarr_dir}/{ROI_name}_Mitochondria.zarr')
    cell_path = (f'{zarr_dir}/{ROI_name}_Cell.zarr')

    mito_da = da.from_zarr(mitochondria_path, chunks=chunk_size).astype(np.int16)
    cell_da = da.from_zarr(cell_path, chunks=chunk_size).astype(np.int8)

    # Lazy multiplication
    relevant_mito = mito_da * cell_da
    labels = relevant_mito.compute()

    return(labels)

def major_axis_from_binary(mask, SPACING = np.array([10,10,50])):
    labels, n = ndimage.label(mask)
    if n > 1:
        sizes = np.bincount(labels.ravel())
        sizes[0] = 0
        mask = labels == sizes.argmax()

    coords = np.argwhere(mask).astype(float)
    coords_phys = coords * SPACING

    centroid = coords_phys.mean(axis=0)

    pca = PCA(n_components=3)
    pca.fit(coords_phys - centroid)

    v1, v2, v3 = pca.components_
    eigvals = pca.explained_variance_

    if v1[2] < 0:
        v1 *= -1
        v2 *= -1
    return v1, v2, v3, eigvals


def compute_mitochondria_properties(labels, SPACING = np.array([10,10,50])):

    voxel_volume = np.prod(SPACING)
    props = regionprops(labels)
    results = []

    for prop in props:
        if prop.area < 10:
            continue
        minx,miny,minz,maxx,maxy,maxz = prop.bbox

        cropped = labels[minx:maxx,miny:maxy,minz:maxz] == prop.label

        major0, _, minor0, eigvals0 = major_axis_from_binary(cropped, SPACING=SPACING)

        volume = prop.area * voxel_volume

        centroid = np.array(prop.centroid) * SPACING

        try:
            verts, faces, _, _ = marching_cubes(
                cropped,
                level=0.5,
                spacing=SPACING
            )

            tri_verts = verts[faces]

            vec1 = tri_verts[:,1] - tri_verts[:,0]
            vec2 = tri_verts[:,2] - tri_verts[:,0]

            cross = np.cross(vec1, vec2)

            surface_area = np.linalg.norm(cross, axis=1).sum() / 2

        except:
            surface_area = np.nan

        results.append({
            "label": prop.label,

            "centroid_x": centroid[0],
            "centroid_y": centroid[1],
            "centroid_z": centroid[2],

            "volume": volume,
            "surface_area": surface_area,

            "orientation_x": major0[0],
            "orientation_y": major0[1],
            "orientation_z": major0[2],

            "minor_axis_x": minor0[0],
            "minor_axis_y": minor0[1],
            "minor_axis_z": minor0[2],

            "eigvals": eigvals0
        })

    return pd.DataFrame(results)


def append_angle_relative_to_nucleus_centroid(df, nucleus_df, ROI_name):

    row = nucleus_df[nucleus_df["ROI"] == ROI_name].iloc[0]

    # Nucleus centroid
    C = np.array([
        row["centroid_x"],
        row["centroid_y"],
        row["centroid_z"]
    ])

    # Mitochondria centroid
    coords = df[["centroid_x","centroid_y","centroid_z"]].values
    w = coords - C

    orientation = df[[
        "orientation_x",
        "orientation_y",
        "orientation_z"
    ]].values

    df["centroid_angle"] = compute_angle_between_vectors(w, orientation)
    df["centroid_distance"] = np.linalg.norm(w, axis=1)

    return df

def append_angle_relative_to_global_directions(df, nucleus_df, ROI_name):
    # row = nucleus_df[nucleus_df["ROI"] == ROI_name].iloc[0]

    orientation = df[[
        "orientation_x",
        "orientation_y",
        "orientation_z"
    ]].values

    df["global_angle_x"] = compute_angle_between_vectors(orientation, np.array([1,0,0]))
    df["global_angle_y"] = compute_angle_between_vectors(orientation, np.array([0,1,0]))
    df["global_angle_z"] = compute_angle_between_vectors(orientation, np.array([0,0,1]))

    return df

# def append_alignment_score(df, nucleus_df, ROI_name):
#     row = nucleus_df[nucleus_df["ROI"] == ROI_name].iloc[0]

#     orientation = df[[
#         "orientation_x",
#         "orientation_y",
#         "orientation_z"
#     ]].values

#     A = np.sum(np.cos(orientation))
#     B = np.sum(np.sin(orientation))
#     R = 1/len(mito_data)*np.sqrt(A**2 + B**2)


def append_spherical_coordinates(df, nucleus_df, ROI_name):
    # Nucleus centroid
    row = nucleus_df[nucleus_df["ROI"] == ROI_name].iloc[0]
    C = np.array([
        row["centroid_x"],
        row["centroid_y"],
        row["centroid_z"]
    ])

    # Mitochondria centroid
    coords = df[["centroid_x","centroid_y","centroid_z"]].values
    rel = coords - C

    r = np.linalg.norm(rel, axis=1)
    rho = np.linalg.norm(rel[:,:2], axis=1)

    theta_xy = np.arctan2(rel[:,1], rel[:,0])
    theta_xz = np.arctan2(rel[:,2], rel[:,0])
    theta_yz = np.arctan2(rel[:,2], rel[:,1])

    phi = np.arctan2(rho, rel[:,2])

    df["centroid_spherical_r"] = r
    df["centroid_cylindrical_rho"] = rho
    df["centroid_theta_xy"] = theta_xy
    df["centroid_theta_xz"] = theta_xz
    df["centroid_theta_yz"] = theta_yz
    df["centroid_phi"] = phi

    return df

######

def compute_closest_point_on_nucleus(mito_centroid, nucleus_centroid, nucleus_da, SPACING = np.array([10,10,50])):

    mito_vox = np.array([mito_centroid[0], mito_centroid[1], mito_centroid[2]]) / SPACING
    nuc_vox  = np.array([nucleus_centroid[0], nucleus_centroid[1], nucleus_centroid[2]]) / SPACING

    min_x = int(np.floor(min(nuc_vox[0], mito_vox[0])))
    max_x = int(np.ceil( max(nuc_vox[0], mito_vox[0])))
    min_y = int(np.floor(min(nuc_vox[1], mito_vox[1])))
    max_y = int(np.ceil( max(nuc_vox[1], mito_vox[1])))
    min_z = int(np.floor(min(nuc_vox[2], mito_vox[2])))
    max_z = int(np.ceil( max(nuc_vox[2], mito_vox[2])))


    cropped = nucleus_da[min_x:max_x, min_y:max_y, min_z:max_z].compute()

    coords = np.argwhere(cropped)

    if coords.shape[0] == 0:
        print('fallback')
        nucleus = nucleus_da.compute()
        coords = np.argwhere(nucleus)
        coords_global = coords
    else:
        coords_global = coords + np.array([min_x, min_y, min_z])

    coords_phys = coords_global * SPACING

    distances = np.linalg.norm(coords_phys - mito_centroid, axis=1)

    idx = np.argmin(distances)

    return coords_phys[idx], distances[idx]


def append_angle_relative_to_nucleus_surface(df, nucleus_df, ROI_name, ZARR_DIR,
                                             chunk_size=(64,64,300), SPACING=[10,10,50]):
    row = nucleus_df[nucleus_df["ROI"] == ROI_name].iloc[0]

    nucleus_centroid = np.array([
        row["centroid_x"],
        row["centroid_y"],
        row["centroid_z"]
    ])
    nucleus_da = da.from_zarr(f'{ZARR_DIR}/{ROI_name}_Nuclei.zarr', chunks=(64,64,300)).astype(np.int8)

    L = len(df)
    for idx, row in tqdm(df.iterrows(), total=L):
        mitochondria_centroid = np.array([row['centroid_x'], row['centroid_y'], row['centroid_z']])
        P, d = compute_closest_point_on_nucleus(mitochondria_centroid,
                                                nucleus_centroid, nucleus_da,
                                                SPACING=SPACING)
        v1 = np.array([mitochondria_centroid[0] - P[0], mitochondria_centroid[1] - P[1], mitochondria_centroid[2] - P[2]])
        v2 = np.array([row['orientation_x'], row['orientation_y'], row['orientation_z']])
        angle = compute_angle_between_vectors(v1, v2)
        df.at[idx, 'surface_angle'] = angle
        df.at[idx, 'surface_distance'] = d

    return df

def append_PCA_distribution(df, nucleus_df, ROI_name):
    row = nucleus_df[nucleus_df["ROI"] == ROI_name].iloc[0]
    v1 = row[["mitocloud_v1_x", "mitocloud_v1_y", "mitocloud_v1_z"]].values
    v2 = row[["mitocloud_v2_x", "mitocloud_v2_y", "mitocloud_v2_z"]].values
    v3 = row[["mitocloud_v3_x", "mitocloud_v3_y", "mitocloud_v3_z"]].values

    vector_dict = {
        'v1':v1,
        'v2':v2,
        'v3':v3
    }

    coords = df[["centroid_x","centroid_y","centroid_z"]].values

    cx = row['centroid_x'].item()
    cy = row['centroid_y'].item()
    cz = row['centroid_z'].item()

    coords = coords - np.array([cx, cy, cz])

    # Project onto PCA plane
    # coords_v1 = coords @ v1
    # coords_v2 = coords @ v2
    # coords_v3 = coords @ v3

    for u_key,w_key in [('v1', 'v2'), ('v1','v3'), ('v2','v3')]:
        theta = []
        rho = []

        u = vector_dict[u_key]
        w = vector_dict[w_key]

        coords_u = coords @ u
        coords_w = coords @ w

        # Angle in intrinsic PCA plane
        for x, y in zip(coords_u, coords_w):
            theta.append(np.arctan2(x, y))
            # Optional radial distance in that plane
            rho.append(np.sqrt(x**2 + y**2))
        df[f'centroid_theta_{u_key}{w_key}'] = theta
        df[f'centroid_rho_{u_key}{w_key}'] = rho

    return df

###################

def compute_mito_cloud_properties(df, nucleus_row):
    # Mitochondria coordinates
    coords = df[["centroid_x","centroid_y","centroid_z"]].values
    centroid = coords.mean(axis=0)

    # Nucleus centroid coordinates
    nucleus_centroid = np.array([
        nucleus_row["centroid_x"],
        nucleus_row["centroid_y"],
        nucleus_row["centroid_z"]
    ])

    pca = PCA(n_components=3)
    pca.fit(coords - centroid)

    eigvals = pca.explained_variance_
    v1, v2, v3 = pca.components_
    # print(v1)
    anisotropy = eigvals[0]/(eigvals[1] + 1e-12)

    return {
        "mitocloud_centroid_x": centroid[0],
        "mitocloud_centroid_y": centroid[1],
        "mitocloud_centroid_z": centroid[2],
        "distance_of_centroids": np.linalg.norm(centroid - nucleus_centroid),
        "mitocloud_L1": eigvals[0],
        "mitocloud_L2": eigvals[1],
        "mitocloud_L3": eigvals[2],
        "anisotropy": anisotropy,
        "mitocloud_v1_x": v1[0],
        "mitocloud_v1_y": v1[1],
        "mitocloud_v1_z": v1[2],
        "mitocloud_v2_x": v2[0],
        "mitocloud_v2_y": v2[1],
        "mitocloud_v2_z": v2[2],
        "mitocloud_v3_x": v3[0],
        "mitocloud_v3_y": v3[1],
        "mitocloud_v3_z": v3[2]
    }