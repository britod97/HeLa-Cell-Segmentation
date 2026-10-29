# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:10:48 2026

@author: adhh334
"""

import os
from pathlib import Path
workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import hela_utils
from tqdm import tqdm

root_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\Cropped")

for ROI_name in tqdm(os.listdir(root_dir)):
    input_dir = root_dir/ROI_name
    output_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\CroppedFiltered") / ROI_name
    
    # print(input_dir)
    # print(output_dir)
    # print("------------------")
    hela_utils.apply_gaussian_filter_dir(input_dir, output_dir)