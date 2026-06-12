# -*- coding: utf-8 -*-
"""
config.py — Landslide Susceptibility Mapping
=============================================
Central configuration block. Edit ONLY this file when adapting the
workflow to a new study area, changing weights, or updating lookup tables.

All other scripts import CONFIG from here:
    from config import CONFIG
"""

CONFIG = {
    # -------------------------------------------------------------------------
    # PATHS  — update base_dir for your machine
    # -------------------------------------------------------------------------
    "base_dir": r"E:\Fellice\Mod2Landslide",

    # Input file names — place all files in base_dir/input/
    "dem":        "ext_DEM_filled.tif",
    "landcover":  "ext_landcov.tif",
    "soil_drain": "ext_soildrain.tif",
    "rock":       "ext_rock.tif",
    "rainfall":   "ext_rainfall.tif",

    # -------------------------------------------------------------------------
    # HYDROLOGY
    # -------------------------------------------------------------------------
    # Flow accumulation threshold for stream network delineation.
    # Cells with accumulation > this value become stream pixels.
    "stream_threshold": 100_000,

    # -------------------------------------------------------------------------
    # AHP WEIGHTS  — must sum to 1.0
    # Derived from Analytic Hierarchy Process (AHP) pairwise comparisons.
    # -------------------------------------------------------------------------
    "weights": {
        "slope":      0.1748,
        "rainfall":   0.0547,
        "landcover":  0.2171,
        "dist":       0.1221,
        "rock":       0.1030,
        "sti_twi":    0.1035,
        "soil_drain": 0.1220,
        "aspect":     0.0451,
        "tri":        0.0317,
        "curvature":  0.0260,
    },

    # -------------------------------------------------------------------------
    # LOOKUP TABLES  — {raster_value: susceptibility_score 1–5}
    # -------------------------------------------------------------------------

    # Soil drainage: higher score = poorer drainage = greater susceptibility
    "soil_drain_map": {
        3: 5,   # Very poorly drained
        2: 4,   # Poorly drained
        1: 2,   # Imperfectly drained
        5: 3,   # Moderately well drained
        4: 1,   # Well drained
    },

    # Land cover (LCDB v6): based on root cohesion, surface roughness,
    # and rainfall interception capacity
    "landcover_map": {
        1:  4,
        2:  2,
        3:  2,
        4:  2,
        5:  2,
        6:  1,
        7:  3,
        8:  4,
        9:  4,
        10: 5,
        11: 5,
        12: 1,
        13: 1,
        14: 3,
        15: 1,
        16: 1,
        17: 3,
        18: 3,
    },

    # Rock / lithology: higher score = weaker, more erodible lithology
    "rock_map": {
        1: 1,
        2: 2,
        3: 3,
        4: 4,
        5: 5,
    },
}

# Sanity check — weights must sum to 1.0
assert abs(sum(CONFIG["weights"].values()) - 1.0) < 1e-6, \
    "CONFIG weights must sum to 1.0"
