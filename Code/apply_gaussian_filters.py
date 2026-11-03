# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:10:48 2026

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

root_dir = Path(r"D:\GitHub Repos\HeLa_Cell_Data\PNGs")

# ''' SMOOTH '''

# for ROI_name in tqdm(os.listdir(root_dir)):
#     input_dir = root_dir/ROI_name
#     output_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\CroppedFiltered") / ROI_name
    
#     # print(input_dir)
#     # print(output_dir)
#     # print("------------------")
#     hela_utils.apply_gaussian_filter_dir(input_dir, output_dir)
  


# fig, ax = plt.subplots(3,2, figsize=(10,10))
# ROI = "Hela_ROI_01_30_876_1665_81"

# for i, l in enumerate(levels):
#     img = Image.open(noisy_dir / ROI / str(l) / "ROI_876_1665_81_z0001.png")
    
#     ax[i//2, i%2].imshow(img, cmap='gray')
#     ax[i//2, i%2].set_title(f"level {l}")
    