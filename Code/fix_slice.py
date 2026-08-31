# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 10:51:23 2026

@author: brito
"""

from scipy.ndimage import distance_transform_edt
import numpy as np
from pathlib import Path
import dask.array as da
from scipy.ndimage import median_filter
from scipy.io import savemat

def smooth_binary_median(mask, size=(3, 3, 3)):
    return median_filter(mask.astype(np.uint8), size=size) > 0


def interpolate_binary_slice(slice_below, slice_above):
    """
    Shape-based interpolation between two binary masks.
    Returns a binary mask representing the 'halfway' shape between them.
    """
    # Signed distance transform: positive inside, negative outside
    def signed_dist(mask):
        return distance_transform_edt(mask) - distance_transform_edt(~mask)

    d_below = signed_dist(slice_below.astype(bool))
    d_above = signed_dist(slice_above.astype(bool))

    d_mid = 0.5 * (d_below + d_above)
    return (d_mid > 0)


ROI_name = 'ROI_01'
ZARR_DIR = Path(r"D:\GitHub Repos\HeLa_Cell_Data\CIL50051\Zarr")
cell_path = (f'{ZARR_DIR}/{ROI_name}_Cell.zarr')

chunk_size = (64,64,255)
cell_da = da.from_zarr(cell_path, chunks=chunk_size).astype(np.int8)


k = 62  # bad slice index — replace with the actual one

cell_np = cell_da.compute()  # if it's not already materialized

new_slice = interpolate_binary_slice(cell_np[:, :, k-1], cell_np[:, :, k+1])
cell_np[:, :, k] = new_slice.astype(cell_np.dtype)

cell_np = smooth_binary_median(cell_np, size=(3, 3, 3))

# Push back into a dask array if needed downstream
cell_da_fixed = da.from_array(cell_np, chunks=cell_da.chunksize)

savemat(
    f'{ZARR_DIR}/{ROI_name}_Cell_fixed.mat',
    {'mask': cell_np.astype(np.uint8)}  # key name is what you'll reference in MATLAB
)

out_path = f'{ZARR_DIR}/{ROI_name}_Cell_fixed.zarr'
cell_da_fixed.to_zarr(out_path, overwrite=True)