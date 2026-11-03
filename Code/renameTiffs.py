# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 13:31:07 2026

@author: adhh334
"""

import os
from pathlib import Path
from tqdm import tqdm

workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

rootdir = Path(r"E:\HeLa\Data\EMPIAR-10094\Noisy\Tiffs")



for subdir in tqdm(os.listdir(rootdir)):
    for noisedir in tqdm(os.listdir(rootdir/subdir)):
        for idx, file in enumerate(sorted(os.listdir(rootdir/subdir/noisedir))):
            sl = f"{int(idx) + 1:04d}"
            
            new_file = f"Hela_ROI_{subdir[15:]}_z{sl}.tif"
            
            os.rename(rootdir/subdir/noisedir/file, rootdir/subdir/noisedir/new_file)