# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 13:49:50 2026

@author: adhh334
"""

import mrcfile
from pathlib import Path

print('---Downsampled Images---')
# Open the MRC file safely using a with-statement
with mrcfile.open(Path(r'E:\HeLa\Data\CIL50051\downsampled_50051.mrc'),
                  header_only=True,
                  permissive=True) as mrc:
    header = mrc.header
    print("nx, ny, nz:", header.nx, header.ny, header.nz)
    print("mx, my, mz:", header.mx, header.my, header.mz)
    print("cella:", header.cella)
    voxel_size = mrc.voxel_size  # convenience property, in Å
    print("voxel_size (Å):", voxel_size)
    print("voxel_size (nm):", voxel_size.x/10, voxel_size.y/10, voxel_size.z/10)
    
    
print('---Fullres Images---')
# Open the MRC file safely using a with-statement
with mrcfile.open(Path(r'E:\HeLa\Data\CIL50051\fullres_50051.mrc'),
                  header_only=True,
                  permissive=True) as mrc:
    header = mrc.header
    print("nx, ny, nz:", header.nx, header.ny, header.nz)
    print("mx, my, mz:", header.mx, header.my, header.mz)
    print("cella:", header.cella)
    voxel_size = mrc.voxel_size  # convenience property, in Å
    print("voxel_size (Å):", voxel_size)
    print("voxel_size (nm):", voxel_size.x/10, voxel_size.y/10, voxel_size.z/10)