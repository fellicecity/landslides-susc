# -*- coding: utf-8 -*-
"""
step05_weighted_sum.py — Landslide Susceptibility Mapping
==========================================================
Computes the final Landslide Susceptibility Index (LSI) from the ten
reclassified input layers. Can be run as a STANDALONE script in a fresh
QGIS session provided all _r.tif files exist in output_dir.

    5.1  File existence check
    5.2  NaN audit (one layer at a time)
    5.3  Weighted sum (memory-safe, one layer at a time)
    5.4  Normalise to 0–100
    5.5  Distribution analysis + auto-select classification method
    5.6  Classify into 5 susceptibility zones
    5.7  Area summary with topographic surface correction
    5.8  Export zone areas to CSV
    5.9  Attach Raster Attribute Table (RAT)
    5.10 Style and load final outputs into QGIS

STANDALONE USE
--------------
Update output_dir and temp_dir below if running after a QGIS restart.
All _r.tif files must exist in output_dir before running.

Run order: step01 → step02 → step03 → step04 → step05
"""

import os
import gc
import csv
import numpy as np
from osgeo import gdal
from scipy import stats
from scipy.cluster.vq import kmeans
import matplotlib.pyplot as plt
from qgis.core import (QgsRasterLayer, QgsProject,
                        QgsPalettedRasterRenderer)
from qgis.PyQt.QtGui import QColor

from config import CONFIG

# =============================================================================
# PATHS  — update these if running standalone after a QGIS restart
# =============================================================================
base       = CONFIG["base_dir"]
output_dir = os.path.join(base, "output")
temp_dir   = os.path.join(base, "temp")

# =============================================================================
# 5.1  FILE EXISTENCE CHECK
# =============================================================================
layer_files = {
    "slope":      "slope_r.tif",
    "rainfall":   "rain_r.tif",
    "landcover":  "landcover_r.tif",
    "dist":       "dist_r.tif",
    "rock":       "rock_r.tif",
    "sti_twi":    "sti_twi_combined.tif",
    "soil_drain": "soil_drain_r.tif",
    "aspect":     "dwi_r.tif",
    "tri":        "tri_r.tif",
    "curvature":  "curv_r.tif",
}

WEIGHTS = CONFIG["weights"]

print("Step 5.1: Checking required files...")
missing = []
for name, filename in layer_files.items():
    path   = os.path.join(output_dir, filename)
    status = "OK" if os.path.exists(path) else "MISSING"
    print(f"  {status} : {filename}")
    if status == "MISSING":
        missing.append(filename)

if missing:
    raise FileNotFoundError(
        f"\n{len(missing)} file(s) missing — re-run steps 01–04 first.\n"
        f"Missing: {missing}"
    )
print("  All files found.\n")

# Normalise weights if not exactly 1.0
total_weight = sum(WEIGHTS.values())
if not np.isclose(total_weight, 1.0):
    print(f"  Normalising weights (sum={total_weight:.4f})...")
    WEIGHTS = {k: v / total_weight for k, v in WEIGHTS.items()}

# =============================================================================
# 5.2  REFERENCE GRID from slope_r.tif
# =============================================================================
ref_ds       = gdal.Open(os.path.join(output_dir, "slope_r.tif"))
cols         = ref_ds.RasterXSize
rows         = ref_ds.RasterYSize
gt           = ref_ds.GetGeoTransform()
proj         = ref_ds.GetProjection()
ref_ds       = None
CELL_SIZE    = abs(gt[1])
nrows, ncols = rows, cols
print(f"Master shape  : ({nrows}, {ncols})  |  Cell size: {CELL_SIZE} m\n")

# =============================================================================
# 5.3  NaN AUDIT
# =============================================================================
print("Step 5.3: NaN audit...")
for name, filename in layer_files.items():
    path = os.path.join(output_dir, filename)
    ds   = gdal.Open(path)
    band = ds.GetRasterBand(1)
    arr  = band.ReadAsArray().astype(np.float32)
    nd   = band.GetNoDataValue()
    ds   = None
    if nd is not None:
        arr[arr == nd] = np.nan
    total_c = arr.size
    n_nan   = int(np.sum(np.isnan(arr)))
    print(f"  {name:<14} valid={total_c-n_nan:>10,}  NaN={n_nan:>10,} "
          f"({100*n_nan/total_c:.2f}%)")
    del arr; gc.collect()

# =============================================================================
# 5.4  WEIGHTED SUM  (memory-safe: one layer at a time)
# =============================================================================
# Each valid pixel accumulates weight proportional to its data coverage.
# Dividing by weight_tracker at the end ensures pixels missing one layer
# still receive a valid (though differently weighted) score, rather than
# being set to NoData.
print("\nStep 5.4: Weighted sum...")

lsi            = np.zeros((nrows, ncols), dtype=np.float32)
weight_tracker = np.zeros((nrows, ncols), dtype=np.float32)

for name, filename in layer_files.items():
    path   = os.path.join(output_dir, filename)
    weight = WEIGHTS[name]

    ds   = gdal.Open(path)
    band = ds.GetRasterBand(1)
    arr  = band.ReadAsArray().astype(np.float32)
    nd   = band.GetNoDataValue()
    ds   = None
    if nd is not None:
        arr[arr == nd] = np.nan

    r   = min(arr.shape[0], nrows)
    c   = min(arr.shape[1], ncols)
    arr = arr[:r, :c]

    valid = ~np.isnan(arr)
    lsi[:r, :c][valid]            += arr[valid] * weight
    weight_tracker[:r, :c][valid] += weight

    print(f"  Added {name:<14} weight={weight:.4f}  valid={np.sum(valid):,}")
    del arr, valid; gc.collect()

covered       = weight_tracker > 0
lsi[covered] /= weight_tracker[covered]
lsi[~covered] = np.nan
del weight_tracker; gc.collect()

print(f"\n  Raw LSI range : {np.nanmin(lsi):.4f} – {np.nanmax(lsi):.4f}")
print(f"  Valid pixels  : {np.sum(~np.isnan(lsi)):,}")

# Save raw LSI
lsi_raw_path = os.path.join(output_dir, "LSI_Raw.tif")
arr_out      = lsi.copy()
arr_out[np.isnan(arr_out)] = -9999
driver = gdal.GetDriverByName("GTiff")
ds     = driver.Create(lsi_raw_path, ncols, nrows, 1, gdal.GDT_Float32)
ds.SetGeoTransform(gt); ds.SetProjection(proj)
ds.GetRasterBand(1).WriteArray(arr_out)
ds.GetRasterBand(1).SetNoDataValue(-9999)
ds.FlushCache(); ds = None
print(f"  Raw LSI saved → {lsi_raw_path}")

# =============================================================================
# 5.5  NORMALISE TO 0–100
# =============================================================================
lsi_min  = np.nanmin(lsi)
lsi_max  = np.nanmax(lsi)
valid    = ~np.isnan(lsi)
lsi_norm = np.full_like(lsi, np.nan)
lsi_norm[valid] = ((lsi[valid] - lsi_min) / (lsi_max - lsi_min)) * 100

lsi_norm_path = os.path.join(output_dir, "LSI_Normalised_0_100.tif")
arr_out       = lsi_norm.copy()
arr_out[np.isnan(arr_out)] = -9999
ds = driver.Create(lsi_norm_path, ncols, nrows, 1, gdal.GDT_Float32)
ds.SetGeoTransform(gt); ds.SetProjection(proj)
ds.GetRasterBand(1).WriteArray(arr_out)
ds.GetRasterBand(1).SetNoDataValue(-9999)
ds.FlushCache(); ds = None
print(f"  Normalised LSI saved → {lsi_norm_path}")
del arr_out; gc.collect()

# =============================================================================
# 5.6  DISTRIBUTION ANALYSIS + AUTO-SELECT CLASSIFICATION METHOD
# =============================================================================
# Method is auto-selected from distribution statistics:
#   Equal Interval  → approximately normal (skew < 0.5, Q-Q R² > 0.98)
#   Jenks           → moderately skewed (skew < 1.0)
#   Quantile        → highly skewed (skew ≥ 1.0)
print("\nStep 5.6: Distribution analysis...")

v_raw    = lsi[~np.isnan(lsi)].ravel()
skewness = stats.skew(v_raw)
kurtosis = stats.kurtosis(v_raw)
mean     = np.mean(v_raw)
median   = np.median(v_raw)

print("\n" + "=" * 45)
print("  LSI DISTRIBUTION SUMMARY — LANDSLIDE")
print("=" * 45)
print(f"  Min      : {v_raw.min():.4f}")
print(f"  Max      : {v_raw.max():.4f}")
print(f"  Mean     : {mean:.4f}")
print(f"  Median   : {median:.4f}")
print(f"  Std Dev  : {np.std(v_raw):.4f}")
print(f"  Skewness : {skewness:.4f}")
print(f"  Kurtosis : {kurtosis:.4f}")
print("=" * 45)

# Histogram + Q-Q plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle("Landslide Susceptibility Index — Distribution Analysis",
             fontsize=13, fontweight="bold")

ax1 = axes[0]
ax1.hist(v_raw, bins=100, color="darkorange", edgecolor="none",
         density=True, alpha=0.7, label="LSI Distribution")
x_line = np.linspace(v_raw.min(), v_raw.max(), 300)
ax1.plot(x_line, stats.norm.pdf(x_line, mean, np.std(v_raw)),
         "b--", linewidth=2, label="Normal Curve")
ax1.axvline(mean,   color="blue",  linestyle="-",  linewidth=1.5,
            label=f"Mean={mean:.3f}")
ax1.axvline(median, color="green", linestyle="--", linewidth=1.5,
            label=f"Median={median:.3f}")
ax1.set_xlabel("LSI Score (raw weighted sum)")
ax1.set_ylabel("Density")
ax1.set_title(f"Histogram  |  Skewness={skewness:.3f}")
ax1.legend(fontsize=8)

ax2 = axes[1]
sample_qq = v_raw[np.random.default_rng(42).choice(
    len(v_raw), size=min(10_000, len(v_raw)), replace=False)]
(osm, osr_vals), (slope_qq, intercept_qq, r) = stats.probplot(sample_qq, dist="norm")
ax2.scatter(osm, osr_vals, s=1, alpha=0.3, color="darkorange", label="Data")
ax2.plot(osm, slope_qq * np.array(osm) + intercept_qq,
         "b--", linewidth=2, label=f"Normal line (R²={r**2:.3f})")
ax2.set_xlabel("Theoretical Quantiles")
ax2.set_ylabel("Sample Quantiles")
ax2.set_title("Q-Q Plot  |  Points on line = Normal")
ax2.legend(fontsize=8)
plt.tight_layout()
plt.show()

if abs(skewness) < 0.5 and r**2 > 0.98:
    METHOD = "equal_interval"
    reason = f"Skewness={skewness:.3f}, R²={r**2:.3f} → approximately normal."
elif abs(skewness) < 1.0:
    METHOD = "jenks"
    reason = f"Skewness={skewness:.3f} (moderate) → natural clusters likely."
else:
    METHOD = "quantile"
    reason = f"Skewness={skewness:.3f} (high) → quantile forces balance."

print(f"\n  AUTO-SELECTED : {METHOD.upper()}")
print(f"  Reason        : {reason}")

# =============================================================================
# 5.7  CLASSIFY INTO 5 ZONES on the 0–100 normalised scale
# =============================================================================
v_norm = lsi_norm[~np.isnan(lsi_norm)].ravel()

if METHOD == "equal_interval":
    breaks = [20.0, 40.0, 60.0, 80.0]
elif METHOD == "jenks":
    samp     = v_norm if len(v_norm) <= 100_000 else \
               np.random.default_rng(42).choice(v_norm, 100_000, replace=False)
    init     = np.percentile(samp, [10, 30, 50, 70, 90]).reshape(-1, 1)
    cents, _ = kmeans(samp.reshape(-1, 1).astype(np.float32),
                      init.astype(np.float32))
    cents    = np.sort(cents.ravel())
    breaks   = list((cents[:-1] + cents[1:]) / 2)
elif METHOD == "quantile":
    breaks = list(np.percentile(v_norm, [20, 40, 60, 80]))

print(f"\n  {METHOD} breaks (0–100): {[f'{b:.2f}' for b in breaks]}")

lsi_zones             = np.full_like(lsi_norm, np.nan, dtype=np.float32)
valid_mask            = ~np.isnan(lsi_norm)
lsi_zones[valid_mask] = (np.digitize(v_norm, breaks) + 1).astype(np.float32)

lsi_zones_path = os.path.join(output_dir, "LSI_Final_Zones.tif")
arr_out        = lsi_zones.copy()
arr_out[np.isnan(arr_out)] = -9999
ds = driver.Create(lsi_zones_path, ncols, nrows, 1, gdal.GDT_Float32)
ds.SetGeoTransform(gt); ds.SetProjection(proj)
ds.GetRasterBand(1).WriteArray(arr_out)
ds.GetRasterBand(1).SetNoDataValue(-9999)
ds.FlushCache(); ds = None

labels      = {1: "Very Low", 2: "Low", 3: "Moderate", 4: "High", 5: "Very High"}
total       = int(np.sum(valid_mask))
full_breaks = [0.0] + breaks + [100.0]

print(f"\n  {'Zone':<5} {'Class':<12} {'Pixels':>10} {'%':>7}  Range (0–100)")
print("  " + "-" * 55)
for z, label in labels.items():
    count = int(np.sum(lsi_zones == z))
    pct   = count / total * 100 if total > 0 else 0
    print(f"  {z}   {label:<12} {count:>10,}  {pct:>6.1f}%  "
          f"{full_breaks[z-1]:.2f} – {full_breaks[z]:.2f}")
print(f"  {'Total':<17} {total:>10,}  100.00%")

# =============================================================================
# 5.8  AREA SUMMARY WITH TOPOGRAPHIC SURFACE CORRECTION
# =============================================================================
# Flat pixel area underestimates true surface area on steep slopes.
# Correction: surface_area = flat_area / cos(slope_degrees)
print("\nStep 5.8: Area summary with topographic correction...")

slope_snap_path = os.path.join(temp_dir, "slope_snapped_to_lsi.tif")

_ds   = gdal.Open(lsi_zones_path)
_gt   = _ds.GetGeoTransform()
_proj = _ds.GetProjection()
_cols = _ds.RasterXSize
_rows = _ds.RasterYSize
_ds   = None
_xmin = _gt[0]; _ymax = _gt[3]
_xmax = _xmin + _cols * _gt[1]
_ymin = _ymax + _rows * _gt[5]

gdal.Warp(slope_snap_path,
          os.path.join(output_dir, "slope_r.tif"),
          format="GTiff",
          outputBounds=(_xmin, _ymin, _xmax, _ymax),
          width=_cols, height=_rows,
          dstSRS=_proj,
          resampleAlg=gdal.GRA_Bilinear,
          srcNodata=-9999, dstNodata=-9999)

def _read(path):
    ds     = gdal.Open(path)
    arr    = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    gt2    = ds.GetGeoTransform()
    ds     = None
    if nodata is not None:
        arr[arr == nodata] = np.nan
    arr[np.abs(arr) >= 1e30] = np.nan
    return arr, gt2

zones_arr, zone_gt = _read(lsi_zones_path)
slope_deg, _       = _read(slope_snap_path)

flat_area_m2    = abs(zone_gt[1]) * abs(zone_gt[5])
cos_slope       = np.cos(np.deg2rad(slope_deg))
cos_slope[cos_slope <= 0] = np.nan
surface_area_m2 = flat_area_m2 / cos_slope

results           = {}
total_surface_km2 = 0.0
total_flat_km2    = 0.0

for zone_val, label in labels.items():
    zone_mask   = np.abs(zones_arr - zone_val) < 0.5
    valid_z     = zone_mask & ~np.isnan(surface_area_m2)
    pixel_count = int(np.sum(zone_mask))
    flat_km2    = (pixel_count * flat_area_m2) / 1e6
    surface_km2 = float(np.nansum(surface_area_m2[valid_z])) / 1e6
    results[zone_val] = {"label": label, "pixels": pixel_count,
                         "flat_km2": flat_km2, "surface_km2": surface_km2}
    total_surface_km2 += surface_km2
    total_flat_km2    += flat_km2

print("\n" + "=" * 70)
print(f"  {'Zone':<5} {'Class':<12} {'Pixels':>10} "
      f"{'Flat (km²)':>12} {'Surface (km²)':>14} {'% Total':>9}")
print("=" * 70)
for zone_val, res in results.items():
    pct = res["surface_km2"] / total_surface_km2 * 100 if total_surface_km2 > 0 else 0
    print(f"  {zone_val:<5} {res['label']:<12} {res['pixels']:>10,} "
          f"{res['flat_km2']:>12.4f} {res['surface_km2']:>14.4f} {pct:>8.1f}%")
print("=" * 70)
print(f"  {'TOTAL':<17} {sum(r['pixels'] for r in results.values()):>10,} "
      f"{total_flat_km2:>12.4f} {total_surface_km2:>14.4f} {'100.0%':>9}")

# =============================================================================
# 5.9  EXPORT CSV
# =============================================================================
csv_path = os.path.join(output_dir, "LSI_Zone_Areas.csv")
with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Zone", "Class", "Pixels", "Flat_km2",
                     "Surface_km2", "Pct_Total", "Method"])
    for zone_val, res in results.items():
        pct = res["surface_km2"] / total_surface_km2 * 100 if total_surface_km2 > 0 else 0
        writer.writerow([zone_val, res["label"], res["pixels"],
                         round(res["flat_km2"],    4),
                         round(res["surface_km2"], 4),
                         round(pct, 2), METHOD])
print(f"\n  CSV saved → {csv_path}")

# =============================================================================
# 5.10  RASTER ATTRIBUTE TABLE (RAT)
# =============================================================================
def attach_rat(raster_path, zone_labels):
    """Attach zone labels, counts, and percentages as a RAT to a GeoTIFF."""
    ds   = gdal.Open(raster_path, gdal.GA_Update)
    band = ds.GetRasterBand(1)
    arr  = band.ReadAsArray().astype(np.float32)
    nd   = band.GetNoDataValue()
    if nd is not None:
        arr[arr == nd] = np.nan
    total_valid = int(np.sum(~np.isnan(arr)))
    counts      = {z: int(np.sum(np.abs(arr - z) < 0.5)) for z in zone_labels}

    rat = gdal.RasterAttributeTable()
    rat.CreateColumn("Value",      gdal.GFT_Integer, gdal.GFU_MinMax)
    rat.CreateColumn("Class",      gdal.GFT_String,  gdal.GFU_Name)
    rat.CreateColumn("Count",      gdal.GFT_Integer, gdal.GFU_PixelCount)
    rat.CreateColumn("Percentage", gdal.GFT_Real,    gdal.GFU_Generic)

    for row_i, (zone_val, label) in enumerate(sorted(zone_labels.items())):
        count = counts[zone_val]
        pct   = round(100 * count / total_valid, 2) if total_valid > 0 else 0.0
        rat.SetValueAsInt   (row_i, 0, zone_val)
        rat.SetValueAsString(row_i, 1, label)
        rat.SetValueAsInt   (row_i, 2, count)
        rat.SetValueAsDouble(row_i, 3, pct)

    band.SetDefaultRAT(rat)
    ds.FlushCache(); ds = None
    print(f"  RAT attached: {os.path.basename(raster_path)}")

attach_rat(lsi_zones_path, labels)

# =============================================================================
# 5.11  STYLE AND LOAD FINAL OUTPUTS INTO QGIS
# =============================================================================
def load_raster(path, name):
    layer = QgsRasterLayer(path, name)
    if not layer.isValid():
        print(f"  FAILED: {name}"); return None
    QgsProject.instance().addMapLayer(layer)
    print(f"  Loaded: {name}")
    return layer

def style_zones(layer):
    """
    Apply standard 5-class susceptibility colour ramp:
    blue (Very Low) → green → yellow → orange → red (Very High).
    """
    classes = [
        QgsPalettedRasterRenderer.Class(1, QColor("#2b83ba"), "Very Low"),
        QgsPalettedRasterRenderer.Class(2, QColor("#abdda4"), "Low"),
        QgsPalettedRasterRenderer.Class(3, QColor("#ffffbf"), "Moderate"),
        QgsPalettedRasterRenderer.Class(4, QColor("#fdae61"), "High"),
        QgsPalettedRasterRenderer.Class(5, QColor("#d7191c"), "Very High"),
    ]
    layer.setRenderer(QgsPalettedRasterRenderer(
        layer.dataProvider(), 1, classes))
    layer.triggerRepaint()

print("\nLoading final outputs into QGIS...")
load_raster(lsi_raw_path,   "LSI Raw")
load_raster(lsi_norm_path,  "LSI Normalised 0-100")
zones_layer = load_raster(lsi_zones_path, "LSI Final Zones")
if zones_layer:
    style_zones(zones_layer)

from qgis.utils import iface
iface.mapCanvas().refresh()

print("\n=== Step 05 complete ===")
print("=== Landslide Susceptibility Mapping — workflow finished ===")
