# -*- coding: utf-8 -*-
"""
utils.py — Landslide Susceptibility Mapping
============================================
Shared helper functions used across all step scripts.
Import in each step with:
    from utils import (read_raster, write_raster, align_and_clip,
                       verify_raster, load_raster, ...)

Functions
---------
GDAL I/O
    read_raster         Read a GeoTIFF to a NumPy float32 array
    write_raster        Write a NumPy array to a GeoTIFF
    get_cell_size       Return pixel width from a raster's GeoTransform
    get_projection      Extract WKT projection from a raster or vector file
    verify_raster       Print min/max and confirm a file was written

Alignment
    align_and_lidar     Warp any raster to exactly match the reference DEM
    align_and_clip      Warp to a pre-defined locked extent string

Reclassification
    memory_safe_reclass Fixed-threshold reclassification (slope only)
    percentile_reclass  Equal-frequency 5-class reclassification
    jenks_reclass       Natural Breaks via 1-D k-means (scipy)
    dist_reclass        Distance-to-stream proximity inversion
    curv_reclass        Curvature zero-split reclassification
    apply_lookup        Categorical lookup table reclassification

QGIS
    load_raster         Add a raster layer to the QGIS canvas
"""

import os
import numpy as np
from osgeo import gdal, ogr, osr
import processing

gdal.UseExceptions()


# =============================================================================
# GDAL I/O HELPERS
# =============================================================================

def read_raster(path):
    """
    Read a single-band GeoTIFF into a NumPy float32 array.

    Replaces the NoData sentinel with NaN so all downstream NumPy operations
    handle missing pixels transparently.

    Parameters
    ----------
    path : str  Absolute path to the input GeoTIFF.

    Returns
    -------
    arr    : np.ndarray (rows, cols) float32 — pixel values; NoData → NaN.
    nodata : float or None — original NoData sentinel from band metadata.
    gt     : tuple — GDAL GeoTransform (x_min, cell_w, 0, y_max, 0, -cell_h).
    proj   : str   — WKT projection string.
    """
    ds     = gdal.Open(path)
    band   = ds.GetRasterBand(1)
    arr    = band.ReadAsArray().astype(np.float32)
    nodata = band.GetNoDataValue()
    gt     = ds.GetGeoTransform()
    proj   = ds.GetProjection()
    ds     = None
    if nodata is not None:
        arr[arr == nodata] = np.nan
    return arr, nodata, gt, proj


def write_raster(path, arr, gt, proj, nodata=-9999, dtype=gdal.GDT_Float32):
    """
    Write a NumPy array to a GeoTIFF on disk.

    Converts NaN back to the NoData sentinel before writing. Uses
    os.replace() to avoid Windows file-lock errors when the output
    path is already open in QGIS.

    Parameters
    ----------
    path   : str            Absolute output path (overwritten if exists).
    arr    : np.ndarray     Float32 array; NaN pixels → NoData sentinel.
    gt     : tuple          GDAL GeoTransform.
    proj   : str            WKT projection string.
    nodata : float          Sentinel for NaN pixels. Default -9999.
    dtype  : int            GDAL data type. Default gdal.GDT_Float32.
    """
    arr_out = arr.copy()
    arr_out[np.isnan(arr_out)] = nodata
    rows, cols = arr_out.shape

    driver = gdal.GetDriverByName("GTiff")
    ds     = driver.Create(path, cols, rows, 1, dtype)
    ds.SetGeoTransform(gt)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.WriteArray(arr_out)
    band.SetNoDataValue(nodata)
    ds.FlushCache()
    ds = None
    os.replace(path, path)
    print(f"  Saved: {os.path.basename(path)}")


def get_cell_size(path):
    """
    Return the pixel width of a raster in its native CRS units.

    Parameters
    ----------
    path : str  Path to any GDAL-readable raster.

    Returns
    -------
    float  Pixel width from GeoTransform index [1].
    """
    ds        = gdal.Open(path)
    gt        = ds.GetGeoTransform()
    cell_size = gt[1]
    ds        = None
    return cell_size


def get_projection(path):
    """
    Extract the WKT projection string from a raster or vector file.

    Tries GDAL first (rasters), then OGR (shapefiles, GeoPackage).

    Parameters
    ----------
    path : str  Path to a GDAL- or OGR-readable spatial file.

    Returns
    -------
    str  WKT projection string.

    Raises
    ------
    ValueError  If the projection cannot be read.
    """
    ds = gdal.Open(path)
    if ds:
        proj = ds.GetProjection()
        ds   = None
        return proj
    ds = ogr.Open(path)
    if ds:
        proj = ds.GetLayer().GetSpatialRef().ExportToWkt()
        ds   = None
        return proj
    raise ValueError(f"Cannot read projection from {path}")


def verify_raster(path, name):
    """
    Confirm a raster was written and print its min/max values.

    Parameters
    ----------
    path : str  Path to the raster to verify.
    name : str  Label printed to the console.

    Returns
    -------
    str  path, returned for chaining.

    Raises
    ------
    FileNotFoundError  If the file does not exist.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"{name} not saved — check path: {path}")
    ds   = gdal.Open(path)
    band = ds.GetRasterBand(1)
    band.ComputeStatistics(False)
    vmin = band.GetMinimum()
    vmax = band.GetMaximum()
    ds   = None
    print(f"  {name}: min={vmin:.4f}  max={vmax:.4f}")
    return path


# =============================================================================
# ALIGNMENT HELPERS
# =============================================================================

def align_and_lidar(input_path, reference_path, output_path, categorical=True):
    """
    Reproject and resample a raster to exactly match the reference LiDAR DEM.

    Reads target CRS, origin, pixel size, and extent from reference_path,
    then calls gdal:warpreproject. Raises RuntimeError on shape mismatch.

    Parameters
    ----------
    input_path     : str   Raster to align.
    reference_path : str   Master reference raster (filled DEM).
    output_path    : str   Destination path for the aligned output.
    categorical    : bool  True → Nearest Neighbour; False → Bilinear.

    Returns
    -------
    str  output_path.
    """
    ref_ds   = gdal.Open(reference_path)
    gt       = ref_ds.GetGeoTransform()
    ncols    = ref_ds.RasterXSize
    nrows    = ref_ds.RasterYSize
    proj_wkt = ref_ds.GetProjection()
    ref_ds   = None

    xmin = gt[0];  ymax = gt[3]
    xmax = xmin + ncols * gt[1]
    ymin = ymax + nrows * gt[5]

    srs  = osr.SpatialReference()
    srs.ImportFromWkt(proj_wkt)
    epsg = srs.GetAttrValue("AUTHORITY", 1)

    extent_str      = f"{xmin},{xmax},{ymin},{ymax} [EPSG:{epsg}]"
    resample_method = 0 if categorical else 1

    print(f"  -> Aligning {os.path.basename(input_path)}")
    processing.run("gdal:warpreproject", {
        "INPUT":             input_path,
        "SOURCE_CRS":        None,
        "TARGET_CRS":        f"EPSG:{epsg}",
        "RESAMPLING":        resample_method,
        "TARGET_RESOLUTION": gt[1],
        "TARGET_EXTENT":     extent_str,
        "TARGET_EXTENT_CRS": f"EPSG:{epsg}",
        "NODATA":            -9999,
        "OUTPUT":            output_path,
    })

    out_ds = gdal.Open(output_path)
    if out_ds is None:
        raise RuntimeError(f"align_and_lidar failed — output not created: {output_path}")
    out_rows = out_ds.RasterYSize
    out_cols = out_ds.RasterXSize
    out_ds   = None

    if (out_rows, out_cols) != (nrows, ncols):
        raise RuntimeError(
            f"Shape mismatch for {os.path.basename(input_path)}:\n"
            f"  Expected ({nrows}, {ncols}), got ({out_rows}, {out_cols})"
        )
    print(f"     OK → ({out_rows}, {out_cols})")
    return output_path


def align_and_clip(in_path, out_path, extent_str, epsg,
                   cell_size=1, categorical=False):
    """
    Warp a raster to a pre-defined locked extent and resolution.

    Lighter alternative to align_and_lidar when the target grid parameters
    are already known as variables (EXTENT_STR, EPSG, CELL_SIZE).

    Parameters
    ----------
    in_path    : str    Input raster path.
    out_path   : str    Output raster path.
    extent_str : str    QGIS extent string: "xmin,xmax,ymin,ymax [EPSG:NNNN]"
    epsg       : str    EPSG code string (e.g. "2193").
    cell_size  : float  Target pixel size. Default 1 (1 m LiDAR).
    categorical: bool   True → Nearest Neighbour; False → Bilinear.
    """
    processing.run("gdal:warpreproject", {
        "INPUT"            : in_path,
        "TARGET_CRS"       : f"EPSG:{epsg}",
        "TARGET_EXTENT"    : extent_str,
        "TARGET_EXTENT_CRS": f"EPSG:{epsg}",
        "TARGET_RESOLUTION": cell_size,
        "RESAMPLING"       : 0 if categorical else 1,
        "NODATA"           : -9999,
        "OUTPUT"           : out_path,
    })


# =============================================================================
# RECLASSIFICATION HELPERS
# =============================================================================

def memory_safe_reclass(arr, bins):
    """
    Fixed-threshold reclassification into 1–5 classes.

    Reserved for slope only — geomorphic thresholds are literature-established
    (Yalcin 2008) and should not be data-driven.

    Parameters
    ----------
    arr  : np.ndarray  Input raster (float32, NaN for NoData).
    bins : array-like  Four monotonically increasing break values.

    Returns
    -------
    np.ndarray  Float32 classes 1–5; NaN preserved.
    """
    reclassed = np.digitize(arr, bins) + 1
    reclassed = reclassed.astype(np.float32)
    reclassed[np.isnan(arr)] = np.nan
    return reclassed


def percentile_reclass(arr, ascending=True, label=""):
    """
    Equal-frequency (quantile) reclassification into 5 susceptibility classes.

    Breaks at the 20th, 40th, 60th, and 80th percentiles guarantee ~20% of
    valid pixels per class regardless of raw value range.

    Parameters
    ----------
    arr       : np.ndarray  Input raster (float32, NaN for NoData).
    ascending : bool        True  → higher value = higher class (e.g. slope).
                            False → lower  value = higher class (e.g. DWI).
    label     : str         Console QA label. "" to suppress.

    Returns
    -------
    np.ndarray  Float32 classes 1–5; NaN preserved.
    """
    result     = np.full_like(arr, np.nan, dtype=np.float32)
    valid_mask = ~np.isnan(arr)
    valid_vals = arr[valid_mask]

    breaks = np.percentile(valid_vals, [20, 40, 60, 80])
    result[valid_mask] = np.digitize(valid_vals, breaks) + 1

    if not ascending:
        result[valid_mask] = 6 - result[valid_mask]
    result[~valid_mask] = np.nan

    if label:
        valid_out = result[valid_mask]
        total     = len(valid_out)
        pcts      = [100 * np.sum(valid_out == c) / total for c in range(1, 6)]
        print(f"  {label:<22} breaks: {breaks[0]:.3f} / {breaks[1]:.3f} / "
              f"{breaks[2]:.3f} / {breaks[3]:.3f}")
        print(f"  {'':22} C1:{pcts[0]:.1f}% C2:{pcts[1]:.1f}% "
              f"C3:{pcts[2]:.1f}% C4:{pcts[3]:.1f}% C5:{pcts[4]:.1f}%")
    return result


def jenks_reclass(arr, ascending=True, label="", log_transform=False):
    """
    Natural Breaks reclassification into 5 classes via 1-D k-means.

    Uses scipy.cluster.vq.kmeans (bundled with QGIS). Subsamples 100k pixels
    with seed=42 for reproducibility. log_transform=True applies log10 before
    clustering — mandatory for power-law distributions (STI, flow accumulation).

    Parameters
    ----------
    arr           : np.ndarray  Input raster (float32, NaN for NoData).
    ascending     : bool        True  → higher value = higher class.
                                False → lower  value = higher class.
    label         : str         Console QA label. "" to suppress.
    log_transform : bool        Apply log10 before clustering. Default False.

    Returns
    -------
    np.ndarray  Float32 classes 1–5; NaN preserved.
    """
    from scipy.cluster.vq import kmeans

    result     = np.full_like(arr, np.nan, dtype=np.float32)
    valid_mask = ~np.isnan(arr)
    valid_vals = arr[valid_mask].copy().astype(np.float64)

    if log_transform:
        valid_vals = np.where(valid_vals > 0, valid_vals, np.nan)
        valid_vals = np.log10(valid_vals)
        valid_vals = valid_vals[~np.isnan(valid_vals)]

    if len(valid_vals) > 100_000:
        rng    = np.random.default_rng(42)
        sample = rng.choice(valid_vals, size=100_000, replace=False)
    else:
        sample = valid_vals.copy()

    init_centroids = np.percentile(sample, [10, 30, 50, 70, 90]).reshape(-1, 1)
    centroids, _   = kmeans(sample.reshape(-1, 1).astype(np.float32),
                            init_centroids.astype(np.float32))
    centroids      = np.sort(centroids.ravel())
    inner_breaks   = (centroids[:-1] + centroids[1:]) / 2

    if log_transform:
        orig_valid    = arr[valid_mask]
        classify_vals = np.log10(np.where(orig_valid > 0, orig_valid, np.nan))
        log_nan_mask  = np.isnan(classify_vals)
        classified    = np.where(
            log_nan_mask, np.nan,
            (np.digitize(classify_vals, inner_breaks) + 1).astype(np.float32)
        )
    else:
        classified = (np.digitize(valid_vals, inner_breaks) + 1).astype(np.float32)

    if not ascending:
        classified = 6 - classified

    result[valid_mask]  = classified
    result[~valid_mask] = np.nan

    if label:
        valid_out = result[~np.isnan(result)].ravel()
        total     = len(valid_out)
        pcts      = [100 * np.sum(valid_out == c) / total for c in range(1, 6)]
        print(f"  {label:<22} breaks: {' / '.join(f'{b:.3f}' for b in inner_breaks)}")
        print(f"  {'':22} C1:{pcts[0]:.1f}% C2:{pcts[1]:.1f}% "
              f"C3:{pcts[2]:.1f}% C4:{pcts[3]:.1f}% C5:{pcts[4]:.1f}%")
    return result


def dist_reclass(arr, label="Dist to Stream"):
    """
    Reclassify cost-path distance to stream into 5 susceptibility classes.

    Applies equal-frequency breaks to positive values, then inverts so
    pixels closest to streams receive class 5 (Very High susceptibility).

    Parameters
    ----------
    arr   : np.ndarray  Cost-path distance raster (float32, NaN for NoData).
    label : str         Console QA label.

    Returns
    -------
    np.ndarray  Float32 — class 1 = far, class 5 = near. NaN preserved.
    """
    result     = np.full_like(arr, np.nan, dtype=np.float32)
    valid_mask = ~np.isnan(arr)
    pos_mask   = valid_mask & (arr >= 0)
    pos_vals   = arr[pos_mask]

    breaks  = np.percentile(pos_vals, [20, 40, 60, 80])
    raw     = np.digitize(pos_vals, breaks) + 1
    flipped = 6 - raw
    result[pos_mask] = flipped

    if label:
        valid_out = result[valid_mask]
        total     = len(valid_out)
        pcts = [100 * np.sum(valid_out == c) / total for c in range(1, 6)]
        print(f"  {label}")
        print(f"  breaks: {breaks[0]:.1f}m / {breaks[1]:.1f}m / "
              f"{breaks[2]:.1f}m / {breaks[3]:.1f}m")
        print(f"  C1:{pcts[0]:.1f}% C2:{pcts[1]:.1f}% C3:{pcts[2]:.1f}% "
              f"C4:{pcts[3]:.1f}% C5:{pcts[4]:.1f}%")
    return result


def curv_reclass(arr, label="Curvature (zero-split)"):
    """
    Reclassify profile curvature into 5 susceptibility classes.

    Splits on zero to avoid duplicate breaks on flat terrain:
        Negative (concave) → classes 4–5   (water convergence = higher risk)
        Zero (flat)        → class  3       (neutral)
        Positive (convex)  → classes 1–2   (water divergence = lower risk)
    Convention: GRASS r.slope.aspect (negative = concave).

    Parameters
    ----------
    arr   : np.ndarray  Profile curvature raster (float32, NaN for NoData).
    label : str         Console QA label.

    Returns
    -------
    np.ndarray  Float32 — class 1 = convex, class 5 = concave. NaN preserved.
    """
    result     = np.full_like(arr, np.nan, dtype=np.float32)
    valid_mask = ~np.isnan(arr)
    neg_vals   = arr[valid_mask & (arr < 0)]
    pos_vals   = arr[valid_mask & (arr > 0)]

    if len(neg_vals) > 0:
        neg_break = np.percentile(neg_vals, 50)
        result[valid_mask & (arr < neg_break)]               = 5
        result[valid_mask & (arr >= neg_break) & (arr < 0)]  = 4
    result[valid_mask & (arr == 0)] = 3
    if len(pos_vals) > 0:
        pos_break = np.percentile(pos_vals, 50)
        result[valid_mask & (arr > 0) & (arr <= pos_break)] = 2
        result[valid_mask & (arr > pos_break)]               = 1

    if label:
        valid_out = result[valid_mask]
        total     = len(valid_out)
        pcts      = [100 * np.sum(valid_out == c) / total for c in range(1, 6)]
        print(f"  {label}")
        print(f"  Neg break: {np.percentile(neg_vals, 50):.4f} | "
              f"Pos break: {np.percentile(pos_vals, 50):.4f}")
        print(f"  C1:{pcts[0]:.1f}% C2:{pcts[1]:.1f}% C3:{pcts[2]:.1f}% "
              f"C4:{pcts[3]:.1f}% C5:{pcts[4]:.1f}%")
    return result


def apply_lookup(arr, mapping):
    """
    Reclassify a categorical raster using an expert-assigned lookup table.

    Any pixel value not present in the mapping becomes NaN (unclassified),
    rather than silently receiving a wrong score.

    Parameters
    ----------
    arr     : np.ndarray  Categorical raster (float32).
    mapping : dict        {source_value: susceptibility_score 1–5}
                          Defined in config.py.

    Returns
    -------
    np.ndarray  Float32 scores 1–5; NaN for unmapped values.
    """
    result = np.full_like(arr, np.nan, dtype=np.float32)
    for src_val, dst_val in mapping.items():
        result[arr == src_val] = dst_val
    return result


# =============================================================================
# QGIS HELPER
# =============================================================================

def load_raster(path, name):
    """
    Load a raster file into the active QGIS project and map canvas.

    Parameters
    ----------
    path : str  Absolute path to the GeoTIFF.
    name : str  Display name in the QGIS Layers panel.

    Returns
    -------
    QgsRasterLayer or None.
    """
    layer = QgsRasterLayer(path, name)
    if not layer.isValid():
        print(f"  FAILED to load: {name}")
        return None
    QgsProject.instance().addMapLayer(layer)
    print(f"  Loaded: {name}")
    return layer
