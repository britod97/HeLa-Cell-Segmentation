# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 17:12:50 2026

@author: brito
"""

# -*- coding: utf-8 -*-
import os, sys, subprocess
from pathlib import Path
workdir = Path(r"D:\GitHub Repos\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import numpy as np
import torch
import hdf5storage
from skimage import io
from empanada.config_loaders import read_yaml
from empanada_napari.inference import Engine3d, stack_postprocessing
from empanada_napari.utils import get_configs

MATLAB_CODE_DIR = workdir        # folder containing segment_noise_level.m, segmentNucleiHelaEM_3D.m, ...
MATLAB_EXE = r"C:\Program Files\MATLAB\R2025b\bin\matlab.exe"   # full path to matlab.exe

noisy_dir   = Path(r"D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Tiffs")
matlab_dir  = Path(r"D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\Noisy\Matlab")

LABEL_DIVISOR = 1000000
MIN_SIZE, MIN_EXTENT = 200, 5
CHUNK = (256, 256, 256)
# NOISE_LEVELS = [0, 0.25, 0.5, 1.0, 1.25, 1.5]
NOISE_LEVELS = [0.5]

# 1. Choose the ROIs (edit here, or pass them on the command line)
ROIS = [
    # 'Hela_ROI_07_30_7395_5114_161',
    'Hela_ROI_08_30_6001_6460_221',
    # 'Hela_ROI_09_30_1859_6893_201',
    # 'Hela_ROI_10_30_1965_1378_181',
    # 'Hela_ROI_13_30_2685_3607_221',
    # 'Hela_ROI_26_30_3322_2960_401'
]

if len(sys.argv) > 1:
    ROIS = sys.argv[1:]


def load_stack(roi_dir, extension="tif"):
    files = sorted(Path(roi_dir).glob(f"*.{extension}"))
    vol = np.stack([io.imread(f) for f in files])
    if vol.ndim == 4:
        vol = vol[..., 0]
    return vol


def init_mitonet():
    config = read_yaml(get_configs()["MitoNet_v1"])
    engine = Engine3d(
        config, inference_scale=1, median_kernel_size=3, nms_kernel=3,
        nms_threshold=0.1, confidence_thr=0.5, min_size=MIN_SIZE,
        min_extent=MIN_EXTENT, label_divisor=LABEL_DIVISOR,
        fine_boundaries=False, semantic_only=False, use_gpu=True,
        use_quantized=False, save_panoptic=False, store_url=None,
        chunk_size=CHUNK,
    )
    return config, engine


def segment_mito(roi, lvl, out_path):
    config, engine = init_mitonet()
    lvl_str = f"{round(lvl*100):03d}"
    vol = load_stack(noisy_dir / roi / lvl_str)
    print(f"{roi}/{lvl_str}: volume {vol.shape}, {vol.dtype}", flush=True)

    stack, trackers = engine.infer_on_axis(vol, "xy")
    gen = stack_postprocessing.__wrapped__(
        {"xy": trackers}, None, config,
        label_divisor=LABEL_DIVISOR, min_size=MIN_SIZE, min_extent=MIN_EXTENT,
        dtype=getattr(engine, "dtype", np.uint32), chunk_size=CHUNK,
    )
    for masks, class_name, instances in gen:
        masks = np.transpose(np.asarray(masks), (1, 2, 0))      # (z,y,x) -> (y,x,z)
        hdf5storage.savemat(str(out_path), {"Hela_mito": masks.astype(np.uint32)},
                            format="7.3", matlab_compatible=True, compress=True)

    del engine, stack, trackers
    torch.cuda.empty_cache()


def run_matlab(roi, lvl_str):
    cmd = (f"addpath(genpath('{MATLAB_CODE_DIR}')); "
           f"segment_noise_level('{roi}', '{lvl_str}', '{noisy_dir}', '{matlab_dir}');")
    r = subprocess.run([MATLAB_EXE, "-batch", cmd], cwd=MATLAB_CODE_DIR)
    if r.returncode != 0:
        print(f"MATLAB step failed for {roi}/{lvl_str} (exit {r.returncode})", flush=True)


def main():
    for roi in ROIS:
        out_dir = matlab_dir / roi
        out_dir.mkdir(parents=True, exist_ok=True)
        roi_id = roi[9:11]

        for lvl in NOISE_LEVELS:
            lvl_str = f"{round(lvl*100):03d}"
            print(f"=== {roi} / noise {lvl_str} ===", flush=True)

            # Mitochondria -> .mat
            mito_out = out_dir / f"ROI_{roi_id}_{lvl_str}_Mito.mat"
            if mito_out.exists():
                print(f"Skipping mito {mito_out.name} (exists)", flush=True)
            else:
                try:
                    segment_mito(roi, lvl, mito_out)
                except Exception as e:
                    print(f"FAILED mito {roi}/{lvl_str}: {e!r}", flush=True)
                    torch.cuda.empty_cache()

            # Nuclei / background / cell -> .mat (same folder)
            cell_out   = out_dir / f"ROI_{roi_id}_{lvl_str}_Cell.mat"
            nuclei_out = out_dir / f"ROI_{roi_id}_{lvl_str}_Nuclei.mat"
            if cell_out.exists() and nuclei_out.exists():
                print(f"Skipping MATLAB {lvl_str} (exists)", flush=True)
            else:
                run_matlab(roi, lvl_str)


if __name__ == "__main__":
    main()