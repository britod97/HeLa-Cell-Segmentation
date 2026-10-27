# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 11:38:52 2026

@author: adhh334
"""

from pathlib import Path

import numpy as np
import tifffile
from ncempy.io import dm

input_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\dm4")
output_dir = Path(r"E:\HeLa\Data\EMPIAR-10094\tiffs_8bit")
output_dir.mkdir(exist_ok=True)

for dm4_file in list(input_dir.glob("*.dm4")):
    img = dm.dmReader(str(dm4_file))["data"].astype(np.float64)

    img = (img - img.mean()) / img.std()
    img = 165 + 35 * img
    img = np.clip(img, 0, 255).astype(np.uint8)

    tifffile.imwrite(output_dir / f"{dm4_file.stem}.tif", img)

print("Done!")