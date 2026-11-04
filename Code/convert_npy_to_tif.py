# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:23:38 2026

@author: adhh334
"""

import os
from pathlib import Path
workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import numpy as np
import tifffile

# Directory containing .npy files
input_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\Noisy\Predictions\Hela_ROI_10_30_1965_1378_181")


for npy_file in input_dir.glob("*.npy"):
    arr = np.load(npy_file)

    tif_file = npy_file.with_suffix(".tif")
    tifffile.imwrite(tif_file, arr)

    print(f"Converted {npy_file.name} -> {tif_file.name}")

print("Done.")