# -*- coding: utf-8 -*-
"""
step04_reclassify.py — Landslide Susceptibility Mapping
========================================================
Reclassifies all ten aligned input layers to a 1–5 susceptibility score
using the method appropriate to each variable's distribution. Processes
one layer at a time with gc.collect() between layers to manage memory on
large 1 m LiDAR rasters. Combines reclassified STI and TWI as an equal-
weight average to form a single hydro-erosion index (sti_twi_combined).

Reclassification method summary
--------------------------------
    percentile_reclass  : Slope, DWI/Aspect, Rainfall
                          Equal-area classes; no assumed distribution
    jenks_reclass       : TWI, TRI  (natural clusters)
    jenks_reclass (log) : STI       (power-law distribution)
    curv_reclass        : Profile Curvature  (zero-split)
    dist_reclass        : Distance to Stream (proximity inversion)
    apply_lookup        : Land cover, Soil drainage, Rock  (expert table)

Run order: step01 → step02 → step03 → step04 → step05
Requires:  step03 outputs in temp_dir
"""

import os
import gc
import ctypes
import numpy as np
from osgeo import gdal

from config import CONFIG
from utils import (read_raster, write_raster,
                   percentile_reclass, jenks_reclass,
                   dist_reclass, curv_reclass, apply_lookup)

# =============================================================================
# PATHS
# =============================================================================
base       = CONFIG["base_dir"]
output_dir = os.path.join(base, "output")
temp_dir   = os.path.join(base, "temp")

# =============================================================================
# 4.1  CLEAR MEMORY from previous steps  (Windows only)
# =============================================================================
print("Clearing memory from previous steps...")
for var in ['slope_arr', 'twi_arr', 'sti_arr', 'tri_arr', 'pcurv_arr',
            'dist_arr', 'rain_arr', 'sdrain_arr', 'lcov_arr', 'dwi_arr']:
    if var in dir():
        del globals()[var]
gc.collect()
ctypes.windll.kernel32.SetProcessWorkingSetSize(-1, -1, -1)
print("  Memory cleared.")

# =============================================================================
# 4.2  REFERENCE gt/proj from slope_aligned
# =============================================================================
slope_path = os.path.join(temp_dir, "slope_aligned.tif")
_, _, gt, proj = read_raster(slope_path)

# =============================================================================
# 4.3  RECLASSIFY EACH LAYER  (one at a time — memory safe)
# =============================================================================
print("\nStep 4.3: Reclassifying all layers to 1–5 susceptibility scores...")

# --- SLOPE  (percentile: equal-area, ascending) ---
arr, _, _, _ = read_raster(slope_path)
r = percentile_reclass(arr, ascending=True, label="Slope")
write_raster(os.path.join(output_dir, "slope_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- TWI  (Jenks: natural clusters, ascending) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "twi.tif"))
r = jenks_reclass(arr, ascending=True, label="TWI")
write_raster(os.path.join(output_dir, "twi_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- STI  (Jenks + log10: power-law distribution, ascending) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "sti.tif"))
r = jenks_reclass(arr, ascending=True, label="STI", log_transform=True)
write_raster(os.path.join(output_dir, "sti_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- TRI  (Jenks: natural clusters, ascending) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "tri.tif"))
r = jenks_reclass(arr, ascending=True, label="TRI")
write_raster(os.path.join(output_dir, "tri_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- PROFILE CURVATURE  (zero-split: concave=5, convex=1) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "pcurv_aligned.tif"))
r = curv_reclass(arr, label="Curvature")
write_raster(os.path.join(output_dir, "curv_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- DWI / ASPECT  (percentile: descending — south-facing = lowest = class 5) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "aspect_aligned.tif"))
r = percentile_reclass(arr, ascending=False, label="DWI")
write_raster(os.path.join(output_dir, "dwi_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- DISTANCE TO STREAM  (proximity inversion: near=5, far=1) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "dist_stream.tif"))
r = dist_reclass(arr, label="Dist to Stream")
write_raster(os.path.join(output_dir, "dist_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- RAINFALL  (percentile: ascending — higher rainfall = higher susceptibility) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "rainfall_aligned.tif"))
r = percentile_reclass(arr, ascending=True, label="Rainfall")
write_raster(os.path.join(output_dir, "rain_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- LAND COVER  (expert lookup table from CONFIG) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "landcov_aligned.tif"))
r = apply_lookup(arr, CONFIG["landcover_map"])
write_raster(os.path.join(output_dir, "landcover_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- SOIL DRAINAGE  (expert lookup table from CONFIG) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "soil_drain_aligned.tif"))
r = apply_lookup(arr, CONFIG["soil_drain_map"])
write_raster(os.path.join(output_dir, "soil_drain_r.tif"), r, gt, proj)
del arr, r; gc.collect()

# --- ROCK / LITHOLOGY  (expert lookup table from CONFIG) ---
arr, _, _, _ = read_raster(os.path.join(temp_dir, "rock_aligned.tif"))
r = apply_lookup(arr, CONFIG["rock_map"])
write_raster(os.path.join(output_dir, "rock_r.tif"), r, gt, proj)
del arr, r; gc.collect()

print("\n  All 11 reclassifications saved.")

# =============================================================================
# 4.4  COMBINE STI + TWI  (equal-weight average)
# =============================================================================
# TWI captures soil saturation potential; STI captures erosive transport power.
# Combining them equally avoids double-counting within the weighted sum,
# while retaining information from both hydrological mechanisms.
print("\nStep 4.4: Combining reclassified STI and TWI...")

twi_r, _, twi_gt, twi_proj = read_raster(os.path.join(output_dir, "twi_r.tif"))
sti_r, _, sti_gt, sti_proj = read_raster(os.path.join(output_dir, "sti_r.tif"))

print(f"  TWI range: {np.nanmin(twi_r):.1f} – {np.nanmax(twi_r):.1f}")
print(f"  STI range: {np.nanmin(sti_r):.1f} – {np.nanmax(sti_r):.1f}")

sti_twi_combined = (0.5 * twi_r) + (0.5 * sti_r)
nodata_mask      = np.isnan(twi_r) | np.isnan(sti_r)
sti_twi_combined = np.where(nodata_mask, np.nan, sti_twi_combined)

sti_twi_path = os.path.join(output_dir, "sti_twi_combined.tif")
write_raster(sti_twi_path, sti_twi_combined, twi_gt, twi_proj)
print(f"  Combined range: {np.nanmin(sti_twi_combined):.2f} – {np.nanmax(sti_twi_combined):.2f}")
print("  STI + TWI combined and saved.")

print("\n=== Step 04 complete ===")
print("    Next: run step05_weighted_sum.py")
