# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 13:28:25 2026

@author: adhh334
"""
import os
from pathlib import Path
workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import numpy as np
from skimage import io
from empanada.config_loaders import read_yaml
from empanada_napari.inference import Engine2d
from empanada_napari.utils import get_configs  # lists downloaded models -> yaml paths

model_configs = get_configs()                  # dict like the one in your traceback
config = read_yaml(model_configs['MitoNet_v1'])

# Settings matching the widget values from your traceback
params = dict(
    inference_scale=1,           # downsampling
    confidence_thr=0.2,
    nms_threshold=0.1,           # center_confidence_thr
    nms_kernel=3,                # min_distance_object_centers
    # plus label_divisor, fine_boundaries, etc. as in run_model
)

engine = Engine2d(config, **params)

image = io.imread('ROI_6005_4739_81_z0116.png')
if image.ndim == 3:
    image = image[..., 0]                      # grayscale

seg = engine.infer(image)                      # check what run_model does with the output