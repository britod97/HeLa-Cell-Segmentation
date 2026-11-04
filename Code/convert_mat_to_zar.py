# -*- coding: utf-8 -*-
import os
from pathlib import Path
os.chdir(r"E:\HeLa\HeLa-Cell-Segmentation\Code")

import scipy.io
import numpy as np
from tqdm import tqdm
import h5py
import zarr
from numcodecs import Blosc

# Location of .mat files
matlab_folder = Path(r"E:\HeLa\Data\EMPIAR-10094\Noisy\Matlab\Hela_ROI_08_30_6001_6460_221")

# Location to place .zarr files
zarr_folder = Path(r"E:\HeLa\Data\EMPIAR-10094\Noisy\Zarr")
zarr_folder.mkdir(parents=True, exist_ok=True)

compressor = Blosc(cname='zstd', clevel=7, shuffle=Blosc.BITSHUFFLE)

ROI = "ROI_08"
VALID_NOISE = {"000", "025", "050", "100", "125", "150"}   # noise level suffixes


def load_volume(file_path):
    """Return (volume, key). Works for v7.3 (HDF5) and older .mat files.
    The file contains a single volume, so the variable name is auto-detected."""
    try:
        with h5py.File(file_path, 'r') as f:
            # keep only real datasets (skip '#refs#' etc.), pick the largest
            cands = [k for k, v in f.items()
                     if isinstance(v, h5py.Dataset) and not k.startswith('#') and v.ndim >= 3]
            if not cands:
                raise ValueError(f"No 3D dataset found (keys: {list(f.keys())})")
            key = max(cands, key=lambda k: f[k].size)
            vol = f[key][()]
        # h5py reads MATLAB arrays with reversed axes
        vol = np.transpose(vol, tuple(range(vol.ndim - 1, -1, -1)))
        return vol, key
    except (OSError, NotImplementedError):
        pass  # not HDF5 -> fall back to scipy

    mat = scipy.io.loadmat(file_path)
    cands = [k for k, v in mat.items()
             if not k.startswith('__') and isinstance(v, np.ndarray) and v.ndim >= 3]
    if not cands:
        raise ValueError(f"No 3D array found (keys: {list(mat.keys())})")
    key = max(cands, key=lambda k: mat[k].size)
    return mat[key], key


for mat_file in tqdm(sorted(os.listdir(matlab_folder))):
    if not mat_file.endswith('.mat') or ROI not in mat_file:
        continue

    stem = Path(mat_file).stem
    noise = stem.split('_')[-1]
    if noise not in VALID_NOISE:
        print(f"Skipping {mat_file}: unexpected noise suffix '{noise}'")
        continue

    file_path = matlab_folder / mat_file
    try:
        dset, key = load_volume(file_path)
    except Exception as e:
        print(f"Could not read {file_path}: {e}")
        continue

    if dset.dtype == bool:
        dset = dset.astype(np.uint8)
    print(f"Loaded {mat_file} (var='{key}', shape={dset.shape}, dtype={dset.dtype}, noise={noise})")

    z1 = zarr.open(
        str(zarr_folder / f"{stem}.zarr"),
        mode='w',
        shape=dset.shape,
        chunks=(10, 256, 256),
        dtype=dset.dtype,
        compressor=compressor,
        zarr_format=2
    )

    # Copy in blocks of 10 slices (matches the chunk size along axis 0)
    for i in range(0, dset.shape[0], 10):
        z1[i:i + 10] = dset[i:i + 10]