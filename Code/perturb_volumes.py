# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:19:15 2026

@author: brito
"""

import os
from pathlib import Path
workdir = Path(r"D:\GitHub Repos\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import scipy.io
import hela_utils


ROI = "Hela_ROI_01_30_876_1665_81"
mat_dir = Path(r"D:\GitHub Repos\HeLa_Cell_Data\Matlab")

mitochondria = hela_utils.load_MAT_file(mat_dir/ (ROI+"_Mitochondria.mat"),
                                        "mitochondria")
cell = hela_utils.load_MAT_file(mat_dir/ (ROI+"_Cell.mat"),
                                        "Hela_cell")
nucleus = hela_utils.load_MAT_file(mat_dir/ (ROI+"_Nuclei.mat"),
                                        "Hela_nuclei")