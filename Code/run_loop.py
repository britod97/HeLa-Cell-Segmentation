# -*- coding: utf-8 -*-
import os
from pathlib import Path
workdir = Path(r"E:\HeLa\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import numpy as np
from skimage import io
from empanada.config_loaders import read_yaml
from empanada_napari.inference import Engine3d, stack_postprocessing
from empanada_napari.utils import get_configs

# images_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\Tiffs")
# predictions_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\Predictions")

noisy_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\Noisy\Tiffs")
noisy_predictions_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\Noisy\Predictions")

LABEL_DIVISOR = 1000000          # "Max objects per class in 3D"
MIN_SIZE, MIN_EXTENT = 200, 5
CHUNK = (256, 256, 256)


def load_stack(roi_dir, extension="tif"):
    files = sorted(Path(roi_dir).glob(f"*.{extension}"))   # ..._z0001, _z0002... sort correctly
    vol = np.stack([io.imread(f) for f in files])
    if vol.ndim == 4:                             # RGB-saved images
        vol = vol[..., 0]
    return vol                                    # (z, y, x), uint8

def init_mitonet():
    config = read_yaml(get_configs()["MitoNet_v1"])
    engine = Engine3d(
        config,
        inference_scale=1,
        median_kernel_size=3,
        nms_kernel=3,
        nms_threshold=0.1,
        confidence_thr=0.5,
        min_size=MIN_SIZE,
        min_extent=MIN_EXTENT,
        label_divisor=LABEL_DIVISOR,
        fine_boundaries=False,
        semantic_only=False,
        use_gpu=True,
        use_quantized=False,
        save_panoptic=False,
        store_url=None,
        chunk_size=CHUNK,
    )
    
    return config, engine


def main():
    # predictions_dir.mkdir(parents=True, exist_ok=True)
    noisy_predictions_dir.mkdir(parents=True, exist_ok=True)
    
    # ROIs = os.listdir(images_dir)
    
    noise_levels = [0, 0.25, 0.5, 1.0, 1.25, 1.5]
    # noise_levels = [0, 1.5]
    # noisy_ROIs = sorted(p.name for p in noisy_dir.iterdir() if p.is_dir())
    noisy_ROIs = [
        # 'Hela_ROI_07_30_7395_5114_161',
        'Hela_ROI_08_30_6001_6460_221',
        # 'Hela_ROI_09_30_1859_6893_201',
        'Hela_ROI_10_30_1965_1378_181',
        # 'Hela_ROI_13_30_2685_3607_221',
        # 'Hela_ROI_26_30_3322_2960_401'
    ]
    
    extension = "tif"
    for ROI in noisy_ROIs:
        for lvl in noise_levels:
            
            # Initialise MitoNet
            config, engine = init_mitonet()
            
            # Load Volume
            lvl_str = str(f"{round(lvl*100):03d}")
            subdir = noisy_dir / ROI / lvl_str
            # print(subdir)
            vol = load_stack(subdir, extension=extension)
            print(f"{ROI}: volume {vol.shape}, {vol.dtype}")
            
            # Make out_dir
            out_dir = noisy_predictions_dir / ROI #/ lvl_str
            out_dir.mkdir(parents=True, exist_ok=True)
            out = out_dir / f"ROI_{ROI[9:11]}_{lvl_str}.npy"

            
            # step 1: stack inference on xy (what the widget's stack_inference does)
            stack, trackers = engine.infer_on_axis(vol, "xy")
            trackers_dict = {"xy": trackers}
        
            # step 2: final label volume per class
            gen = stack_postprocessing.__wrapped__(
                trackers_dict, None, config,
                label_divisor=LABEL_DIVISOR,
                min_size=MIN_SIZE,
                min_extent=MIN_EXTENT,
                dtype=getattr(engine, "dtype", np.uint32),
                chunk_size=CHUNK,
            )
            
            # Save output
            for masks, class_name, instances in gen:
                np.save(out, np.asarray(masks).astype(np.uint32))
                # print(f"saved {out.name}: {len(instances)} instances")





if __name__ == "__main__":
    main()