# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:19:15 2026

@author: brito
"""

import os
from pathlib import Path
workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

# import scipy.io 
import hela_utils
import numpy as np
import matplotlib.pyplot as plt
import time
import perturb_utils
from tqdm import tqdm
from scipy.io import savemat

def ROI_name(mat_volume_filename):
    return '_'.join(mat_volume_filename.split('_')[:-1])

savedir = Path(r"E:\HeLa\Data\EMPIAR-10094\Perturbations/mitochondria_scaling")
# os.makedirs(savedir, exist_ok=True)

# spacing = (10,10,50)
spacing = (1,1,5)


basedir = Path(r"E:\HeLa\Data\EMPIAR-10094\Matlab")
files = os.listdir(basedir)

ROIs = np.unique([ROI_name(filename) for filename in files if filename.endswith(".mat")])

for factor in [1.2, 0.8]:
    outdir = savedir / f"{factor}"
    os.makedirs(outdir, exist_ok=True)
    
    for ROI in tqdm(ROIs[1:]):
        
        mitochondria = hela_utils.load_MAT_file(basedir / (ROI+"_Mitochondria.mat"),
                                            "mitochondria")
        print("Loaded Mito")
        cell = hela_utils.load_MAT_file(basedir/ (ROI+"_Cell.mat"),
                                                "Hela_cell")
        print("Loaded Cell")
        nucleus = hela_utils.load_MAT_file(basedir/ (ROI+"_Nuclei.mat"),
                                                "Hela_nuclei")
        print("Loaded Nuc")
        scaled = perturb_utils.scale_mitochondria(mitochondria, cell, nucleus, factor)
    
        savefile = f"{ROI[:11]}_Mitochondria_Scaled_{factor}.mat"
        savemat(outdir/savefile, {"mitochondria":scaled})
    
    
    