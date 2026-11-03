# -*- coding: utf-8 -*-
import os
from pathlib import Path
workdir = Path(r"D:\GitHub Repos\HeLa-Cell-Segmentation\Code")
os.chdir(workdir)

import numpy as np
from skimage import io
from empanada.config_loaders import read_yaml
from empanada_napari.inference import Engine3d, stack_postprocessing
from empanada_napari.utils import get_configs

png_dir = Path(r"D:\GitHub Repos\HeLa_Cell_Data\PNGs")
noisy_png_dir = Path(r"D:\GitHub Repos\HeLa_Cell_Data\Noisy")

save_dir = Path(r"D:\GitHub Repos\HeLa_Cell_Data\EMPIAR-10094\NewMitochondriaPredictions\npy")

LABEL_DIVISOR = 1000000          # "Max objects per class in 3D"
MIN_SIZE, MIN_EXTENT = 200, 5
CHUNK = (256, 256, 256)


def load_stack(roi_dir):
    files = sorted(Path(roi_dir).glob("*.png"))   # ..._z0001, _z0002... sort correctly
    vol = np.stack([io.imread(f) for f in files])
    if vol.ndim == 4:                             # RGB-saved PNGs
        vol = vol[..., 0]
    return vol                                    # (z, y, x), uint8


def main():
    save_dir.mkdir(parents=True, exist_ok=True)

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

    ROI = "Hela_ROI_08_30_6001_6460_221"
    vol = load_stack(png_dir / ROI)
    print(f"{ROI}: volume {vol.shape}, {vol.dtype}")
    
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
    for masks, class_name, instances in gen:
        out = save_dir / f"{ROI}_{class_name}.npy"
        np.save(out, np.asarray(masks).astype(np.uint32))
        print(f"saved {out.name}: {len(instances)} instances")


if __name__ == "__main__":
    main()