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
from scipy import ndimage
import numpy as np

def scale_instance(labels, lid, s, spacing):
    """Shape-preserving scaling of one instance about its centroid."""
    m = labels == lid
    com = np.array(ndimage.center_of_mass(m))
    # affine_transform maps output->input coords; scale in voxel units
    # (isotropic in physical space: the scale factor is the same per axis,
    # so anisotropic spacing doesn't matter for pure scaling)
    inv = np.eye(3) / s
    offset = com - inv @ com
    out = ndimage.affine_transform(m.astype(np.float32), inv, offset=offset, order=1) > 0.5
    new = labels.copy()
    new[m] = 0
    new[out & ((labels == 0) | (labels == lid))] = lid   # don't overwrite neighbours
    return new

ROI = "Hela_ROI_07_30_7395_5114_161"
mat_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\Matlab")

mitochondria = hela_utils.load_MAT_file(mat_dir/ (ROI+"_Mitochondria.mat"),
                                        "mitochondria")
cell = hela_utils.load_MAT_file(mat_dir/ (ROI+"_Cell.mat"),
                                        "Hela_cell")
nucleus = hela_utils.load_MAT_file(mat_dir/ (ROI+"_Nuclei.mat"),
                                        "Hela_nuclei")


mitochondria_in_cell = mitochondria*cell

labels = np.unique(mitochondria_in_cell)

lid = 2277

scale_instance(labels,lid,)