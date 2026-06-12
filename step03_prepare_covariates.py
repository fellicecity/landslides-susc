# -*- coding: utf-8 -*-
"""
step03_prepare_covariates.py — Landslide Susceptibility Mapping
================================================================
Aligns all categorical and environmental layers to the locked DEM extent,
then runs a pre-reclassification shape check to confirm every layer shares
the same pixel grid before the weighted sum in step05.

    3.1  Align categorical layers (land cover, soil drainage, rock)
         → Nearest Neighbour resampling preserves integer class codes
    3.2  Align continuous environmental layer (rainfall)
         → Bilinear resampling for smooth interpolation
    3.3  Pre-reclassification shape check — all 10 layers must match

Run order: step01 → step02 → step03 → step04 → step05
Requires:  step01 and step02 outputs in temp_dir / output_dir
"""

import os
from osgeo import gdal, osr
from qgis.core import QgsRasterLayer, QgsProject

from config import CONFIG
from utils import align_and_clip, load_raster

# =============================================================================
# PATHS AND LOCKED GRID PARAMETERS
# =============================================================================
base       = CONFIG["base_dir"]
input_dir  = os.path.join(base, "input")
output_dir = os.path.join(base, "output")
temp_dir   = os.path.join(base, "temp")

dem_path = os.path.join(input_dir, CONFIG["dem"])

# Input categorical and environmental layers
landcover_path  = os.path.join(input_dir, CONFIG["landcover"])
soil_drain_path = os.path.join(input_dir, CONFIG["soil_drain"])
rock_path       = os.path.join(input_dir, CONFIG["rock"])
rainfall_path   = os.path.join(input_dir, CONFIG["rainfall"])

# Re-derive locked grid from DEM
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
# 3.1  ALIGN CATEGORICAL LAYERS  (Nearest Neighbour)
# =============================================================================
# Nearest Neighbour preserves integer class codes.
# Bilinear interpolation would corrupt them by creating fractional values.
print("Step 3.1: Aligning categorical layers (Nearest Neighbour)...")

landcov_res    = os.path.join(temp_dir, "landcov_aligned.tif")
soil_drain_res = os.path.join(temp_dir, "soil_drain_aligned.tif")
rock_aligned   = os.path.join(temp_dir, "rock_aligned.tif")

align_and_clip(landcover_path,  landcov_res,
               extent_str=EXTENT_STR, epsg=EPSG,
               cell_size=CELL_SIZE, categorical=True)
print("  Land cover aligned.")

align_and_clip(soil_drain_path, soil_drain_res,
               extent_str=EXTENT_STR, epsg=EPSG,
               cell_size=CELL_SIZE, categorical=True)
print("  Soil drainage aligned.")

align_and_clip(rock_path,       rock_aligned,
               extent_str=EXTENT_STR, epsg=EPSG,
               cell_size=CELL_SIZE, categorical=True)
print("  Rock / lithology aligned.")

# =============================================================================
# 3.2  ALIGN RAINFALL (Bilinear — continuous surface)
# =============================================================================
# Rainfall was kriged from HIRDS point data at 2000 m and resampled to 1 m.
# Bilinear interpolation is appropriate for continuous numerical values.
print("\nStep 3.2: Aligning rainfall layer (Bilinear)...")

rainfall_aligned = os.path.join(temp_dir, "rainfall_aligned.tif")
align_and_clip(rainfall_path, rainfall_aligned,
               extent_str=EXTENT_STR, epsg=EPSG,
               cell_size=CELL_SIZE, categorical=False)
print("  Rainfall aligned.")

# =============================================================================
# 3.3  PRE-RECLASSIFICATION SHAPE CHECK
# =============================================================================
# All 10 input layers must share exactly the same (rows, cols) grid before
# the weighted sum. Any mismatch here will cause element-wise addition errors
# downstream and is caught immediately rather than propagating silently.
print("\nStep 3.3: Running pre-reclassification shape check...")

slope_path = os.path.join(temp_dir, "slope_aligned.tif")

_ref_ds       = gdal.Open(slope_path)
EXPECTED_ROWS = _ref_ds.RasterYSize
EXPECTED_COLS = _ref_ds.RasterXSize
_ref_ds       = None
print(f"  Reference grid: {EXPECTED_ROWS} rows × {EXPECTED_COLS} cols")

check_paths = {
    "slope"     : slope_path,
    "pcurv"     : os.path.join(temp_dir, "pcurv_aligned.tif"),
    "tri"       : os.path.join(temp_dir, "tri.tif"),
    "twi"       : os.path.join(temp_dir, "twi.tif"),
    "sti"       : os.path.join(temp_dir, "sti.tif"),
    "dist"      : os.path.join(temp_dir, "dist_stream.tif"),
    "landcover" : landcov_res,
    "soil_drain": soil_drain_res,
    "rainfall"  : rainfall_aligned,
    "rock"      : rock_aligned,
}

print(f"\n  {'Layer':<12} {'Rows':>6} {'Cols':>6}  Status")
print("  " + "-" * 33)
mismatches = []
for name, path in check_paths.items():
    ds   = gdal.Open(path)
    rows = ds.RasterYSize
    cols = ds.RasterXSize
    ds   = None
    ok   = rows == EXPECTED_ROWS and cols == EXPECTED_COLS
    if not ok:
        mismatches.append(name)
    print(f"  {name:<12} {rows:>6} {cols:>6}  {'OK' if ok else '*** MISMATCH ***'}")

if mismatches:
    raise RuntimeError(
        f"\nShape mismatch in: {mismatches}\n"
        f"Expected: ({EXPECTED_ROWS}, {EXPECTED_COLS})\n"
        f"Re-run align_and_clip for the affected layers."
    )
print("\n  All 10 layers aligned — safe to proceed to reclassification.")

# =============================================================================
# 3.4  LOAD OUTPUTS INTO QGIS
# =============================================================================
print("\nLoading step03 outputs into QGIS...")
load_raster(landcov_res,      "landcover_aligned")
load_raster(soil_drain_res,   "soil_drain_aligned")
load_raster(rock_aligned,     "rock_aligned")
load_raster(rainfall_aligned, "rainfall_aligned")

print("\n=== Step 03 complete ===")
print("    Next: run step04_reclassify.py")
