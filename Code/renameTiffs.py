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

rootdir = Path(r"E:\HeLa\Data\EMPIAR-10094\Cropped")

for subdir in tqdm(os.listdir(rootdir)):
    for file in os.listdir(rootdir/subdir):
        sl = f"{int(file[-8:-4]) + 1:04d}"
        
        new_file = f"Hela_ROI_{subdir[15:]}_z{sl}.tif"
        
        os.rename(rootdir/subdir/file, rootdir/subdir/new_file)