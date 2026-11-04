# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 17:11:28 2026

@author: adhh334
"""

import os
from pathlib import Path
# workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
workdir = Path(r"D:\GitHub Repos\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import hela_utils
from tqdm import tqdm
import numpy as np
from PIL import Image
from tifffile import imwrite
import matplotlib.pyplot as plt


def add_noise(img, level, z, sd):
    """sd: std of the clean ROI (raw units), computed once per ROI."""
    x = img.astype(np.float32) + level * sd * z
    if np.issubdtype(img.dtype, np.integer):
        info = np.iinfo(img.dtype)
        x = np.clip(np.rint(x), info.min, info.max)
    return x.astype(img.dtype)

root_dir = Path(r"D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\PNGs")
    

''' NOISE '''
noisy_dir = Path(r"D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Tiffs")


# levels = [0, 0.25, 0.5]
levels = [1.0, 1.25, 1.5, 1.75, 2.0]


ROIs = os.listdir(root_dir)
ROIs_to_process = ["07", "08", "09", "10", "13", "26"]
ROIs = [ROI for ROI in ROIs if ROI[9:11] in ROIs_to_process]


for ROI_name in tqdm(ROIs):
    files = sorted(os.listdir(root_dir / ROI_name))

    # one global std per ROI (every 10th slice is plenty)
    sd = np.std(np.stack([np.array(Image.open(root_dir / ROI_name / f))
                          for f in files[::10]]).astype(np.float32))

    for lvl in levels:
        os.makedirs(noisy_dir / ROI_name / str(f"{round(lvl*100):03d}"), exist_ok=True)

    for i, file in enumerate(files):
        # print(file)
        img = np.array(Image.open(root_dir / ROI_name / file))
        rng = np.random.default_rng([0, i])           # different field per slice
        z = rng.standard_normal(img.shape).astype(np.float32)

        for lvl in levels:                            # same z across levels -> nested noise
            out = add_noise(img, lvl, z, sd)
            imwrite(noisy_dir / ROI_name / str(f"{round(lvl*100):03d}") / f"{file[:-4]}.tif", out)   # outputs as .tif

# Progressive series: reuse the same z so each level is the previous one plus more noise