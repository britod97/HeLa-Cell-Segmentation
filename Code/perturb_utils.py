# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:26:00 2026

@author: adhh334
"""

from scipy import ndimage
from tqdm import tqdm
import numpy as np
import hela_utils


def _scale_instances(labels, factors, extra=2):
    """factors: {lid: s}. Scales each instance about its centroid, bbox-cropped."""
    out = labels.copy()                      # one copy total, not one per lid
    boxes = ndimage.find_objects(labels)         # one pass; boxes[lid-1] is that lid's bbox

    for lid, s in tqdm(factors.items()):
        box = boxes[lid - 1]
        if box is None:
            continue
        # pad the bbox enough to hold the grown instance
        sl = tuple(
            slice(max(b.start - p, 0), min(b.stop + p, n))
            for b, n in zip(box, labels.shape)
            for p in [int(np.ceil(max(s - 1, 0) * (b.stop - b.start) / 2)) + extra]
        )
        crop = out[sl]                       # view: edits write through to `out`
        m = crop == lid
        com = np.array(ndimage.center_of_mass(m))
        inv = np.eye(3) / s
        new_m = ndimage.affine_transform(
            m.astype(np.float32), inv, offset=com - inv @ com, order=1
        ) > 0.5
        crop[m] = 0
        crop[new_m & (crop == 0)] = lid      # don't overwrite neighbours
    return out


def scale_mitochondria(mito_volume, cell_volume, nuc_volume, factor):
    mitochondria_in_cell = mito_volume*cell_volume
    # border_ids = hela_utils.touches_border(mitochondria_in_cell)
    lid_list = list(np.unique(mitochondria_in_cell)[1:])   # First will always be 0
    factor_list = [factor for _ in lid_list]
    factors = dict(zip(lid_list, factor_list))
    
    scaled = _scale_instances(mitochondria_in_cell, factors)
    return scaled*~nuc_volume


def grow_to_factor(labels, lid, factor, spacing, tol=0.005):
    """Grow/shrink one instance to factor * original volume via EDT thresholding."""
    m = labels == lid
    v0 = m.sum()
    free = (labels == 0) | m
    d_in = ndimage.distance_transform_edt(m, sampling=spacing)
    d_out = ndimage.distance_transform_edt(~m, sampling=spacing)
    signed = d_in - d_out                      # >0 inside, <0 outside
    lo, hi = -signed.max(), signed.max()
    for _ in range(40):
        t = (lo + hi) / 2
        cand = (signed > t) & free
        r = cand.sum() / v0
        if abs(r - factor) < tol: break
        lo, hi = (t, hi) if r > factor else (lo, t)
    new = labels.copy()
    new[m] = 0
    new[cand] = lid
    return new