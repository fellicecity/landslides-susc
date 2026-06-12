# -*- coding: utf-8 -*-
"""
step02_advanced_indices.py — Landslide Susceptibility Mapping
==============================================================
Computes advanced terrain and erosion indices from the filled DEM and
cleaned flow accumulation, then aligns all continuous terrain layers
to the locked DEM extent.

    2.1  Flow accumulation cleaning (negative GRASS values → NaN)
    2.2  Terrain Ruggedness Index (TRI) — GDAL tool
    2.3  Topographic Wetness Index (TWI) — NumPy
    2.4  Sediment Transport Index (STI) — NumPy (log-transform recommended)
    2.5  Cost-path distance to streams — GRASS r.cost
    2.6  Alignment of all continuous terrain layers to locked extent

Run order: step01 → step02 → step03 → step04 → step05
Requires:  step01 outputs in temp_dir
"""

import os
import numpy as np
import warnings
from osgeo import gdal, osr
from qgis.core import QgsRasterLayer, QgsProject
import processing

from config import CONFIG
from utils import (read_raster, write_raster, align_and_clip,
                   verify_raster, load_raster)

# =============================================================================
# PATHS AND LOCKED GRID PARAMETERS
# =============================================================================
base       = CONFIG["base_dir"]
input_dir  = os.path.join(base, "input")
output_dir = os.path.join(base, "output")
temp_dir   = os.path.join(base, "temp")

dem_path        = os.path.join(input_dir, CONFIG["dem"])
dem_filled_path = dem_path   # update if using a separately filled DEM

# Paths from step01
slope_path    = os.path.join(temp_dir, "slope.tif")
flow_acc_path = os.path.join(temp_dir, "flow_acc.tif")

# Re-derive locked grid from DEM (consistent with step01)
ref_ds   = gdal.Open(dem_path)
gt       = ref_ds.GetGeoTransform()
ncols    = ref_ds.RasterXSize
nrows    = ref_ds.RasterYSize
proj_wkt = ref_ds.GetProjection()
ref_ds   = None

LEFT   = gt[0];  TOP   = gt[3]
RIGHT  = LEFT + ncols * gt[1]
BOTTOM = TOP  + nrows * gt[5]

srs        = osr.SpatialReference()
srs.ImportFromWkt(proj_wkt)
EPSG       = srs.GetAttrValue("AUTHORITY", 1)
CELL_SIZE  = gt[1]
EXTENT_STR = f"{LEFT},{RIGHT},{BOTTOM},{TOP} [EPSG:{EPSG}]"

# =============================================================================
# 2.1  FLOW ACCUMULATION CLEANING
# =============================================================================
# GRASS r.watershed returns negative values at drainage divides.
# Set <= 0 to NaN and floor remaining values at 1 to avoid log(0) in TWI.
print("Step 2.1: Cleaning flow accumulation...")

acc_arr, _, acc_gt, acc_proj = read_raster(flow_acc_path)
acc_clean = np.where(acc_arr <= 0, np.nan, acc_arr)
acc_clean = np.where(np.isnan(acc_clean), np.nan, np.maximum(acc_clean, 1))

acc_clean_path = os.path.join(output_dir, "flow_acc_clean.tif")
write_raster(acc_clean_path, acc_clean, acc_gt, acc_proj)
print("  Flow accumulation cleaned.")

# =============================================================================
# 2.2  TERRAIN RUGGEDNESS INDEX (TRI)
# =============================================================================
# Computed with the built-in GDAL TRI tool, then clipped to locked extent.
print("\nStep 2.2: Calculating TRI...")

tri_raw_path = os.path.join(temp_dir, "tri_raw.tif")
processing.run("gdal:triterrainruggednessindex", {
    "INPUT"   : dem_filled_path,
    "Z_Factor": 1,
    "OUTPUT"  : tri_raw_path,
})
tri_path = os.path.join(temp_dir, "tri.tif")
align_and_clip(tri_raw_path, tri_path,
               extent_str=EXTENT_STR, epsg=EPSG,
               cell_size=CELL_SIZE, categorical=False)
verify_raster(tri_path, "TRI")

# =============================================================================
# 2.3  TOPOGRAPHIC WETNESS INDEX (TWI)
# =============================================================================
# Formula: ln( (FlowAcc × cell_size) / tan(Slope) )
# tan(slope) floored at 0.001 to avoid division by zero on flat terrain.
print("\nStep 2.3: Calculating TWI...")

slope_arr, _, slp_gt, slp_proj = read_raster(slope_path)
acc_arr,   _, acc_gt, acc_proj = read_raster(acc_clean_path)

slope_rad = np.deg2rad(slope_arr)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    tan_slope_safe = np.where(np.tan(slope_rad) > 0, np.tan(slope_rad), 0.001)
    acc_safe_twi   = np.where(acc_arr > 0, acc_arr, 0.001) * CELL_SIZE
    twi_arr        = np.log(acc_safe_twi / tan_slope_safe)

twi_raw_path = os.path.join(temp_dir, "twi_raw.tif")
twi_path     = os.path.join(temp_dir, "twi.tif")
write_raster(twi_raw_path, twi_arr, acc_gt, acc_proj)
align_and_clip(twi_raw_path, twi_path,
               extent_str=EXTENT_STR, epsg=EPSG,
               cell_size=CELL_SIZE, categorical=False)
verify_raster(twi_path, "TWI")

# =============================================================================
# 2.4  SEDIMENT TRANSPORT INDEX (STI)
# =============================================================================
# Formula: ((FlowAcc × cell_size) / 22.13)^0.6 × (sin(Slope) / 0.0896)^1.3
# Constants follow Moore & Burch (1986).
# sin(slope) floored at 0.001; log_transform used in reclassification (step04).
print("\nStep 2.4: Calculating STI...")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sin_slope_safe = np.where(np.sin(slope_rad) > 0, np.sin(slope_rad), 0.001)
    acc_safe_sti   = acc_arr + 1
    term1          = np.power((acc_safe_sti * CELL_SIZE) / 22.13, 0.6)
    term2          = np.power(sin_slope_safe / 0.0896, 1.3)
    sti_arr        = term1 * term2

sti_raw_path = os.path.join(temp_dir, "sti_raw.tif")
sti_path     = os.path.join(temp_dir, "sti.tif")
write_raster(sti_raw_path, sti_arr, acc_gt, acc_proj)
align_and_clip(sti_raw_path, sti_path,
               extent_str=EXTENT_STR, epsg=EPSG,
               cell_size=CELL_SIZE, categorical=False)
verify_raster(sti_path, "STI")

# =============================================================================
# 2.5  COST-PATH DISTANCE TO STREAMS (GRASS r.cost)
# =============================================================================
# Streams = cells where cleaned flow accumulation > stream_threshold.
# Cost surface = slope in radians: steep terrain is topographically "further"
# from a stream than the same Euclidean distance over flat ground.
print("\nStep 2.5: Computing cost-path distance to streams...")

streams_arr  = np.where(acc_clean > CONFIG["stream_threshold"], 1, np.nan)
streams_path = os.path.join(temp_dir, "streams.tif")
write_raster(streams_path, streams_arr, acc_gt, acc_proj)
print(f"  Stream network defined at threshold = {CONFIG['stream_threshold']:,} cells")

cost_surface_arr  = np.where(slope_rad > 0, slope_rad, 0.001)
cost_surface_path = os.path.join(temp_dir, "cost_surface.tif")
write_raster(cost_surface_path, cost_surface_arr, slp_gt, slp_proj)

dist_raw_path = os.path.join(temp_dir, "dist_raw.tif")
processing.run("grass7:r.cost", {
    "input"                          : cost_surface_path,
    "start_raster"                   : streams_path,
    "output"                         : dist_raw_path,
    "memory"                         : 300,
    "GRASS_REGION_PARAMETER"         : dem_filled_path,
    "GRASS_REGION_CELLSIZE_PARAMETER": CELL_SIZE,
})

dist_path = os.path.join(temp_dir, "dist_stream.tif")
align_and_clip(dist_raw_path, dist_path,
               extent_str=EXTENT_STR, epsg=EPSG,
               cell_size=CELL_SIZE, categorical=False)
verify_raster(dist_path, "Distance to Stream")

# =============================================================================
# 2.6  ALIGN CONTINUOUS TERRAIN LAYERS TO LOCKED EXTENT
# =============================================================================
print("\nStep 2.6: Aligning continuous terrain layers...")

aspect_path = os.path.join(temp_dir, "aspect.tif")
pcurv_path  = os.path.join(temp_dir, "pcurvature.tif")

continuous_layers = [
    (os.path.join(temp_dir, "slope.tif"),      os.path.join(temp_dir, "slope_aligned.tif")),
    (aspect_path,                               os.path.join(temp_dir, "aspect_aligned.tif")),
    (pcurv_path,                                os.path.join(temp_dir, "pcurv_aligned.tif")),
]
for in_path, out_path in continuous_layers:
    align_and_clip(in_path, out_path,
                   extent_str=EXTENT_STR, epsg=EPSG,
                   cell_size=CELL_SIZE, categorical=False)
    print(f"  Aligned: {os.path.basename(out_path)}")

# =============================================================================
# 2.7  LOAD OUTPUTS INTO QGIS
# =============================================================================
print("\nLoading step02 outputs into QGIS...")
load_raster(tri_path,                                 "tri")
load_raster(twi_path,                                 "twi")
load_raster(sti_path,                                 "sti")
load_raster(dist_path,                                "dist_stream")
load_raster(os.path.join(temp_dir,"slope_aligned.tif"), "slope_aligned")
load_raster(os.path.join(temp_dir,"pcurv_aligned.tif"), "pcurv_aligned")

print("\n=== Step 02 complete ===")
print("    Next: run step03_prepare_covariates.py")
