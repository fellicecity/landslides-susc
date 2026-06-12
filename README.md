# Landslide Susceptibility Mapping — QGIS Python Workflow

A reproducible, open-source workflow for AHP-based landslide susceptibility
mapping using freely available geospatial data and the QGIS Python Console.
Produces a classified 1 m resolution Landslide Susceptibility Index (LSI)
raster by combining ten terrain, environmental, and anthropogenic conditioning
factors.

Developed for **Hawke's Bay, New Zealand** (NZTM2000 / EPSG:2193),
anchored by the Cyclone Gabrielle (2023) landslide inventory.

---

## Repository Structure

```
landslide-susceptibility-mapping/
├── config.py                       ← All weights, paths, and lookup tables
├── utils.py                        ← Shared helper functions
├── step01_terrain_derivatives.py   ← Slope, DWI, Curvature, Flow accumulation
├── step02_advanced_indices.py      ← TRI, TWI, STI, Distance to stream
├── step03_prepare_covariates.py    ← Categorical layer alignment + shape check
├── step04_reclassify.py            ← Reclassify all layers to 1–5
├── step05_weighted_sum.py          ← Weighted sum → classify → export
├── README.md
├── LICENSE
└── .gitignore
```

---

## Conditioning Factors and AHP Weights

| Factor | Input | Weight |
|---|---|---|
| Land Cover (LCDB v6) | `ext_landcov.tif` | 0.2171 |
| Slope | Derived from LiDAR DEM | 0.1748 |
| Distance to Stream | GRASS r.cost | 0.1221 |
| Soil Drainage | `ext_soildrain.tif` | 0.1220 |
| STI + TWI Combined | Derived from LiDAR DEM | 0.1035 |
| Rock / Lithology | `ext_rock.tif` | 0.1030 |
| Rainfall (HIRDS) | `ext_rainfall.tif` | 0.0547 |
| Directional Wetness Index | Derived from aspect | 0.0451 |
| Terrain Ruggedness Index | Derived from LiDAR DEM | 0.0317 |
| Profile Curvature | Derived from LiDAR DEM | 0.0260 |

---

## How to Run

### Requirements

| Software | Notes |
|---|---|
| QGIS ≥ 3.16 LTR | Run from Plugins → Python Console → Show Editor |
| GRASS GIS | Bundled with QGIS |
| Python, NumPy, SciPy, Matplotlib, GDAL | All bundled with QGIS — no pip install needed |

### Setup

1. Clone or download this repository
2. Open `config.py` and update `base_dir` to your local path
3. Place all input `.tif` files in `base_dir\input\`
4. Open QGIS → `Plugins` → `Python Console` → click **Show Editor**

### Run order

```
step01_terrain_derivatives.py
    ↓
step02_advanced_indices.py
    ↓
step03_prepare_covariates.py
    ↓
step04_reclassify.py
    ↓
step05_weighted_sum.py
```

> **Tip:** `step05_weighted_sum.py` can run as a **standalone** script in a
> fresh QGIS session if all `_r.tif` files already exist in `output_dir`.

---

## Input Files

Place all in `base_dir\input\`:

| File | Description |
|---|---|
| `ext_DEM_filled.tif` | Sink-filled LiDAR DEM (1 m, EPSG:2193) |
| `ext_landcov.tif` | LCDB v6 land cover class codes |
| `ext_soildrain.tif` | Soil drainage class raster |
| `ext_rock.tif` | Rock / lithology class raster |
| `ext_rainfall.tif` | Mean annual rainfall (kriged from HIRDS) |

---

## Output Files

All written to `base_dir\output\`:

| File | Description |
|---|---|
| `slope_r.tif` … `curv_r.tif` | Factor rasters reclassified to 1–5 |
| `sti_twi_combined.tif` | Equal-weight STI + TWI average |
| `LSI_Raw.tif` | Raw weighted sum |
| `LSI_Normalised_0_100.tif` | Normalised 0–100 |
| `LSI_Final_Zones.tif` | 5-class susceptibility map + embedded RAT |
| `LSI_Zone_Areas.csv` | Zone areas (flat + topographic correction) |

---

## Adapting to a New Study Area

Edit **only** `config.py` — update `base_dir`, `weights`, and the three
lookup tables. No changes to any step script required.

---

## References

- Yalcin (2008). GIS-based landslide susceptibility mapping. *CATENA*, 72(1).
- Pourghasemi et al. (2012). Landslide susceptibility using SVM. *IJDRS*, 3(2).
- Tien Bui et al. (2016). Landslide susceptibility mapping. *Geomorphology*, 256.
- Massey et al. (2025). Cyclone Gabrielle landslide inventory. GNS Science.

---

## License

MIT License — see `LICENSE` file.
