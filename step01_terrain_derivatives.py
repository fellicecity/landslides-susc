# -*- coding: utf-8 -*-
"""
step01_terrain_derivatives.py — Landslide Susceptibility Mapping
=================================================================
Computes primary terrain derivatives from the filled LiDAR DEM:

    - Slope (degrees)
    - Aspect → Directional Wetness Index  DWI = cos(aspect_radians)
    - Profile Curvature
    - Flow Accumulation and Drainage Direction (GRASS r.watershed)

Also reads and locks the DEM extent, EPSG, and cell size used by all
subsequent steps.

Run order: step01 → step02 → step03 → step04 → step05
"""

import os
import numpy as np
from osgeo import gdal, osr
from qgis.core import QgsRasterLayer, QgsProject
import processing

from config import CONFIG
from utils import (read_raster, write_raster, get_cell_size,
                   get_projection, load_raster)

# =============================================================================
# FOLDER STRUCTURE
# =============================================================================
base       = CONFIG["base_dir"]
input_dir  = os.path.join(base, "input")
output_dir = os.path.join(base, "output")
temp_dir   = os.path.join(base, "temp")

for d in [input_dir, output_dir, temp_dir]:
    os.makedirs(d, exist_ok=True)

dem_path = os.path.join(input_dir, CONFIG["dem"])

# =============================================================================
# READ DEM METADATA — locks the master grid for all steps
# =============================================================================
print("Reading DEM metadata...")

ref_ds   = gdal.Open(dem_path)
gt       = ref_ds.GetGeoTransform()
ncols    = ref_ds.RasterXSize
nrows    = ref_ds.RasterYSize
proj_wkt = ref_ds.GetProjection()
ref_ds   = None

LEFT   = gt[0]
TOP    = gt[3]
RIGHT  = LEFT + ncols * gt[1]
BOTTOM = TOP  + nrows * gt[5]   # gt[5] is negative

srs        = osr.SpatialReference()
srs.ImportFromWkt(proj_wkt)
EPSG       = srs.GetAttrValue("AUTHORITY", 1)
CELL_SIZE  = gt[1]
EXTENT_STR = f"{LEFT},{RIGHT},{BOTTOM},{TOP} [EPSG:{EPSG}]"

print(f"  EPSG       : {EPSG}")
print(f"  Cell size  : {CELL_SIZE} m")
print(f"  Extent     : {EXTENT_STR}")

# =============================================================================
# 1.1  DEM SINK FILLING  (uncomment if DEM is not pre-filled)
# =============================================================================
# If your DEM already has sinks filled, skip this block.
# CONFIG["dem"] should point to the filled DEM directly.
#
# dem_filled_path = os.path.join(temp_dir, "dem_filled.tif")
# processing.run("grass7:r.fill.dir", {
#     "input":    dem_path,
#     "format":   0,
#     "output":   dem_filled_path,
#     "direction":os.path.join(temp_dir, "dir_temp.tif"),
#     "areas":    os.path.join(temp_dir, "areas_temp.tif"),
#     "GRASS_REGION_PARAMETER": dem_path,
# })
# ─────────────────────────────────────────────────────────────────────────────
# If already filled, just point to the DEM directly:
dem_filled_path = dem_path

# =============================================================================
# 1.2  SLOPE, ASPECT, AND PROFILE CURVATURE (GRASS r.slope.aspect)
# =============================================================================
# Curvature sign convention (GRASS):
#   Negative = concave (water converges) = higher susceptibility → class 5
#   Positive = convex  (water disperses) = lower  susceptibility → class 1
print("\nStep 1.2: Calculating slope, aspect, and profile curvature...")

slope_path = os.path.join(temp_dir, "slope.tif")
aspect_path = os.path.join(temp_dir, "aspect.tif")
pcurv_path  = os.path.join(temp_dir, "pcurvature.tif")

processing.run("grass7:r.slope.aspect", {
    "elevation":  dem_filled_path,
    "slope":      slope_path,
    "aspect":     aspect_path,
    "pcurvature": pcurv_path,
    "format":     0,
    "precision":  0,
    "zscale":     1,
    "min_slope":  0,
    "GRASS_REGION_PARAMETER": dem_filled_path,
})
print("  Slope, aspect, curvature done.")

# =============================================================================
# 1.3  DIRECTIONAL WETNESS INDEX (DWI) from Aspect
# =============================================================================
# DWI = cos(aspect_radians)
# In the Southern Hemisphere (NZ):
#   North (0°)   → cos =  1  → drier  → lower susceptibility
#   South (180°) → cos = -1  → wetter → higher susceptibility (class 5)
# ascending=False in step04 maps lowest DWI → class 5 ✓
print("\nStep 1.3: Transforming aspect to Directional Wetness Index (DWI)...")

aspect_arr, _, asp_gt, asp_proj = read_raster(aspect_path)
aspect_rad = np.deg2rad(aspect_arr)
dwi_arr    = np.cos(aspect_rad)

dwi_path = os.path.join(temp_dir, "dwi.tif")
write_raster(dwi_path, dwi_arr, asp_gt, asp_proj)
print("  DWI saved.")

# =============================================================================
# 1.4  FLOW ACCUMULATION AND DRAINAGE DIRECTION (GRASS r.watershed)
# =============================================================================
# r.watershed returns signed accumulation values — negative at flow divides.
# These are cleaned in step02 before computing TWI/STI.
print("\nStep 1.4: Calculating flow accumulation and drainage direction...")

flow_acc_path = os.path.join(temp_dir, "flow_acc.tif")
flow_dir_path = os.path.join(temp_dir, "flow_dir.tif")

processing.run("grass7:r.watershed", {
    "elevation":    dem_filled_path,
    "threshold":    CONFIG["stream_threshold"],
    "accumulation": flow_acc_path,
    "drainage":     flow_dir_path,
    "-s":           True,
    "GRASS_REGION_PARAMETER": dem_filled_path,
})
print("  Flow accumulation and direction done.")

# =============================================================================
# 1.5  LOAD OUTPUTS INTO QGIS
# =============================================================================
print("\nLoading step01 outputs into QGIS...")
load_raster(slope_path,    "slope")
load_raster(aspect_path,   "aspect")
load_raster(pcurv_path,    "profile_curvature")
load_raster(dwi_path,      "dwi")
load_raster(flow_acc_path, "flow_acc")
load_raster(flow_dir_path, "flow_dir")

print("\n=== Step 01 complete ===")
print("    Next: run step02_advanced_indices.py")
