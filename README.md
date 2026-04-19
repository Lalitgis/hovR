# hovR <img src="inst/logo.png" align="right" height="120" alt=""/>

> End-to-end R tools for drone hyperspectral time-series analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg)](https://cran.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.2.0-green.svg)]()

---

## What it solves

| Without hovR | With hovR |
|---|---|
| Draw 500 plot polygons by hand in QGIS | `segment_plots(raster)` — done in seconds |
| Click panels in every flight image, note DN values, run spreadsheet | `detect_panels()` + `calibrate_reflectance()` — one call per flight |
| Manually check each image for blur / cloud / saturation | `flight_qc(stack)` → HTML report with accept/reject per flight |
| Analyse each flight date independently | `integrate_season()` — seasonal AUC per pixel across all flights |
| Try flight dates one by one to find the best predictor | `optimal_window()` — tests every date and window automatically |

---

## Installation

```r
# Install from GitHub
remotes::install_github("dronehsi/hovR")
```

**Dependencies** (installed automatically): `terra`, `sf`, `ggplot2`, `dplyr`, `cli`, `lubridate`

---

## Modules

### 1 — Flight quality assessment

Run before any analysis. Checks every flight for blur, saturation,
irradiance instability, and no-data fraction. Outputs a structured
accept/caution/reject decision and an HTML report.

```r
qc <- flight_qc(my_flight_stack)
flight_qc_report(qc, output_file = "qc_season_2024.html")
```

### 2 — Calibration panel detection

Locates grey calibration panels automatically in flight images using
brightness thresholding and spectral flatness scoring. Applies
empirical-line reflectance calibration without manual clicking.

```r
panels <- detect_panels(raw_raster, n_panels = 2,
                         panel_size_m2 = c(0.25, 0.25))
refl   <- calibrate_reflectance(raw_raster, panels,
            known_reflectance = c(panel_1 = 0.50, panel_2 = 0.20))
```

### 3 — Plot auto-segmentation

Detects trial plot boundaries directly from the drone image.
Handles regular grids, staggered rows, and irregular layouts.
Exports an `sf` polygon layer ready for QGIS or modelling.

```r
plots <- segment_plots(refl, wavelengths = wl,
                        method = "combined",
                        expected_rows = 20, expected_cols = 10)
plot_segmentation_qc(plots, refl)          # visual QC
sf::st_write(plots, "trial_plots.gpkg")    # export to QGIS
```

### 4 — Temporal integration

The core module. Treats a growing season of drone flights as one
coherent temporal dataset. Computes seasonal VI integrals (AUC),
temporal derivatives, and optimal prediction windows.

```r
# Stack all calibrated flights
fs <- flight_stack(raster_list, dates = flight_dates,
                    wavelengths = wl)

# Tag phenological stages from field notebook
fs <- phenostage_tag(fs, bbch_obs = data.frame(
  date = c("2024-03-20", "2024-05-01", "2024-06-10"),
  bbch = c("21", "55", "83")
))

# Compute vegetation indices
vi <- compute_vi(fs, indices = c("NDVI", "NDRE", "LCI", "MCARI"))

# Seasonal integral (AUC) — one number per pixel capturing the whole story
auc <- integrate_season(vi, index = "NDVI")

# Rate of change between flights — detect stress onset
rates <- vi_derivative(vi, index = "NDRE", smooth = TRUE)

# Find which dates/windows best predict yield
best <- optimal_window(vi, ground_truth = yield_df,
                        trait = "yield_t_ha", plots = plots)
```

### 5 — Ground truth fusion

Spatially extracts VI statistics within each plot polygon and matches
them to field-measured traits by date and plot ID.

```r
plot_data <- extract_plot(vi, plots,
                           indices = c("NDVI", "NDRE"),
                           stats   = c("mean", "sd"))

fused     <- fuse_ground(plot_data, ground_truth = chl_df,
                          traits = "spad_chl", max_lag_days = 3)

plot_vi_traits(fused, vi_col = "mean", trait_col = "spad_chl",
               colour_by = "date", fit = "lm")
```

---

## Complete workflow

```
Raw drone images
      │
      ▼
flight_qc()          ← check blur, saturation, irradiance per flight
      │
      ▼
detect_panels()      ← find calibration panels automatically
calibrate_reflectance() ← convert DN → reflectance
      │
      ▼
segment_plots()      ← auto-detect plot boundaries (no QGIS needed)
      │
      ▼
flight_stack()       ← assemble chronological raster stack
compute_vi()         ← NDVI, NDRE, LCI, MCARI, ... (wavelength-agnostic)
phenostage_tag()     ← attach BBCH growth stages
      │
      ▼
integrate_season()   ← seasonal AUC per pixel
vi_derivative()      ← rate of change between flights
optimal_window()     ← which dates/windows predict traits best?
      │
      ▼
extract_plot()       ← per-plot VI statistics
fuse_ground()        ← match with field measurements
      │
      ▼
lm() / plsr() / tidymodels  ← trait models
```

---

## Built-in vegetation indices

| Index | Full name | Key application |
|---|---|---|
| NDVI | Normalised Difference Vegetation Index | General greenness |
| NDRE | Normalised Difference Red-Edge | Chlorophyll in dense canopies |
| VARI | Visible Atmospherically Resistant Index | RGB-based greenness |
| MCARI | Modified Chlorophyll Absorption RI | Chlorophyll content |
| mNDVI705 | Modified NDVI 705/750 nm | Canopy chlorophyll |
| LCI | Leaf Chlorophyll Index | Leaf chlorophyll |
| GNDVI | Green NDVI | Chlorophyll, LAI |
| CRI1 | Carotenoid Reflectance Index 1 | Carotenoid content |
| WBI | Water Band Index | Canopy water content |
| PSRI | Plant Senescence RI | Senescence / ripening |

Add your own with `define_index()`.

---

## Research gaps filled

- **Seasonal VI integration (AUC)** — no existing R package provides this; shown to outperform single-date analysis for yield prediction (Killian et al., 2025)
- **Automated calibration panel detection** — currently done by hand per flight; costs a full day per season
- **Plot auto-segmentation** — currently costs days of manual QGIS work for large breeding trials
- **Automated flight QC** — currently invisible; researchers do not know which flights are borderline until analysis fails

---

## License

MIT © 2026 hovR Authors
