library(testthat)
library(terra)
library(sf)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_field_raster <- function(nrow = 80, ncol = 120, n_bands = 5,
                               add_panels = FALSE) {
  wl <- seq(450, 850, length.out = n_bands)
  r  <- terra::rast(nrows = nrow, ncols = ncol, nlyr = n_bands)
  terra::ext(r) <- c(0, 12, 0, 8)        # 12m wide, 8m tall
  terra::crs(r) <- "EPSG:32632"

  # Simulate a 4-row x 6-column trial: plots are bright, alleys are dark
  m <- matrix(0.05, nrow = nrow, ncol = ncol)  # dark background (alleys)

  # Create plot regions
  row_bounds <- round(seq(0, nrow, length.out = 6))  # 5 plot rows
  col_bounds <- round(seq(0, ncol, length.out = 8))  # 7 plot cols

  for (ri in seq_len(length(row_bounds) - 1)) {
    for (ci in seq_len(length(col_bounds) - 1)) {
      r_from <- row_bounds[ri] + 2; r_to <- row_bounds[ri + 1] - 2
      c_from <- col_bounds[ci] + 2; c_to <- col_bounds[ci + 1] - 2
      if (r_from < r_to && c_from < c_to)
        m[r_from:r_to, c_from:c_to] <- 0.4 + runif(1, 0, 0.1)
    }
  }

  # If panels requested, add two bright uniform regions
  if (add_panels) {
    # Panel 1 (top-left corner): high, spectrally flat
    m[2:8, 2:10] <- 0.50
    # Panel 2 (bottom-right corner): moderate, spectrally flat
    m[(nrow - 8):(nrow - 2), (ncol - 10):(ncol - 2)] <- 0.20
  }

  for (b in seq_len(n_bands)) {
    noise <- matrix(rnorm(nrow * ncol, 0, 0.01), nrow, ncol)
    band_scale <- if (wl[b] > 700) 1.5 else 1.0  # NIR brighter
    terra::values(r[[b]]) <- as.vector(t(pmax(0, pmin(1, m * band_scale + noise))))
  }
  list(raster = r, wavelengths = wl)
}

# ---------------------------------------------------------------------------
# Tests: segment_plots()
# ---------------------------------------------------------------------------

test_that("segment_plots() returns an sf object", {
  d <- make_field_raster()
  polys <- segment_plots(
    d$raster, wavelengths = d$wavelengths,
    method = "profile", min_plot_area_m2 = 0.1
  )
  expect_s3_class(polys, "sf")
})

test_that("segment_plots() returns polygons with required columns", {
  d <- make_field_raster()
  polys <- segment_plots(
    d$raster, wavelengths = d$wavelengths,
    method = "profile", min_plot_area_m2 = 0.1
  )
  expect_true(all(c("plot_id", "row_idx", "col_idx",
                     "area_m2", "quality") %in% names(polys)))
})

test_that("segment_plots() plot_ids are unique", {
  d <- make_field_raster()
  polys <- segment_plots(
    d$raster, wavelengths = d$wavelengths,
    method = "profile", min_plot_area_m2 = 0.1
  )
  expect_equal(length(unique(polys$plot_id)), nrow(polys))
})

test_that("segment_plots() quality scores are in [0, 1]", {
  d <- make_field_raster()
  polys <- segment_plots(
    d$raster, wavelengths = d$wavelengths,
    method = "profile", min_plot_area_m2 = 0.1
  )
  q <- polys$quality[!is.na(polys$quality)]
  expect_true(all(q >= 0 & q <= 1))
})

test_that("segment_plots() works with RGB raster (no wavelengths)", {
  r <- terra::rast(nrows = 60, ncols = 80, nlyr = 3)
  terra::ext(r) <- c(0, 8, 0, 6)
  terra::crs(r) <- "EPSG:32632"
  # Green channel higher in plot areas
  m <- matrix(0.1, 60, 80)
  m[5:25, 5:35] <- 0.5; m[5:25, 45:75] <- 0.5
  m[35:55, 5:35] <- 0.5; m[35:55, 45:75] <- 0.5
  for (b in 1:3) terra::values(r[[b]]) <- as.vector(t(m)) + rnorm(60*80, 0, 0.02)
  terra::values(r[[2]]) <- terra::values(r[[2]]) * 1.3  # green boost

  polys <- segment_plots(r, method = "profile", min_plot_area_m2 = 0.5)
  expect_s3_class(polys, "sf")
  expect_gte(nrow(polys), 1)
})

test_that("segment_plots() errors on missing wavelengths for >3 bands", {
  d <- make_field_raster(n_bands = 10)
  expect_error(segment_plots(d$raster), regexp = "wavelengths")
})

test_that("segment_plots() custom id_prefix is applied", {
  d <- make_field_raster()
  polys <- segment_plots(
    d$raster, wavelengths = d$wavelengths,
    method = "profile", min_plot_area_m2 = 0.1, id_prefix = "PLOT"
  )
  expect_true(all(grepl("^PLOT", polys$plot_id)))
})

test_that("refine_plot_boundaries() returns sf with shift_m column", {
  d <- make_field_raster()
  polys <- segment_plots(
    d$raster, wavelengths = d$wavelengths,
    method = "profile", min_plot_area_m2 = 0.1
  )
  refined <- refine_plot_boundaries(polys, d$raster[[1]], max_shift_m = 0.2)
  expect_s3_class(refined, "sf")
  expect_true("shift_m" %in% names(refined))
  expect_equal(nrow(refined), nrow(polys))
})

test_that("plot_segmentation_qc() returns a ggplot", {
  d <- make_field_raster()
  polys <- segment_plots(
    d$raster, wavelengths = d$wavelengths,
    method = "profile", min_plot_area_m2 = 0.1
  )
  p <- plot_segmentation_qc(polys, d$raster, show_ids = FALSE)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# Tests: detect_panels()
# ---------------------------------------------------------------------------

test_that("detect_panels() returns sf with required columns", {
  d <- make_field_raster(add_panels = TRUE)
  panels <- detect_panels(
    d$raster, wavelengths = d$wavelengths,
    n_panels = 1, panel_size_m2 = c(0.01, 2),
    brightness_quantile = 0.90
  )
  expect_s3_class(panels, "sf")
  expect_true(all(c("panel_id", "brightness", "spectral_cv",
                     "confidence", "area_m2") %in% names(panels)))
})

test_that("detect_panels() returns correct number of panels", {
  d <- make_field_raster(add_panels = TRUE)
  panels <- detect_panels(
    d$raster, wavelengths = d$wavelengths,
    n_panels = 1, panel_size_m2 = c(0.01, 2),
    brightness_quantile = 0.90
  )
  expect_lte(nrow(panels), 1)
})

test_that("detect_panels() panel_ids are correctly prefixed", {
  d <- make_field_raster(add_panels = TRUE)
  panels <- detect_panels(
    d$raster, wavelengths = d$wavelengths,
    n_panels = 1, panel_size_m2 = c(0.01, 2),
    brightness_quantile = 0.90
  )
  if (nrow(panels) > 0)
    expect_true(all(grepl("^panel_", panels$panel_id)))
})

test_that("calibrate_reflectance() returns SpatRaster with values in [0,1]", {
  d <- make_field_raster(add_panels = TRUE)
  panels <- detect_panels(
    d$raster, wavelengths = d$wavelengths,
    n_panels = 1, panel_size_m2 = c(0.01, 2),
    brightness_quantile = 0.90
  )
  skip_if(nrow(panels) == 0, "No panels detected in synthetic raster")

  calib <- calibrate_reflectance(
    raster            = d$raster,
    panels            = panels,
    known_reflectance = c(panel_1 = 0.50)
  )
  expect_s4_class(calib, "SpatRaster")
  expect_equal(terra::nlyr(calib), terra::nlyr(d$raster))
  vals <- terra::values(calib, na.rm = TRUE)
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("calibrate_reflectance() errors on unmatched panel IDs", {
  d <- make_field_raster(add_panels = TRUE)
  panels <- detect_panels(
    d$raster, wavelengths = d$wavelengths,
    n_panels = 1, panel_size_m2 = c(0.01, 2),
    brightness_quantile = 0.90
  )
  skip_if(nrow(panels) == 0)
  expect_error(
    calibrate_reflectance(d$raster, panels,
                          known_reflectance = c(wrong_id = 0.5)),
    regexp = "matching panel"
  )
})

test_that("panel_qc_plot() returns a list with map and spectra", {
  d <- make_field_raster(add_panels = TRUE)
  panels <- detect_panels(
    d$raster, wavelengths = d$wavelengths,
    n_panels = 1, panel_size_m2 = c(0.01, 2),
    brightness_quantile = 0.90
  )
  skip_if(nrow(panels) == 0)
  plots <- panel_qc_plot(d$raster, panels, wavelengths = d$wavelengths)
  expect_true(is.list(plots))
  expect_s3_class(plots$map,     "ggplot")
  expect_s3_class(plots$spectra, "ggplot")
})

# ---------------------------------------------------------------------------
# Tests: flight_qc()
# ---------------------------------------------------------------------------

make_clean_raster <- function() {
  d <- make_field_raster(nrow = 50, ncol = 60, n_bands = 4)
  d$raster
}

make_blurry_raster <- function() {
  r <- make_clean_raster()
  # Constant value = zero Laplacian variance
  terra::values(r) <- 0.5
  r
}

make_nodata_raster <- function() {
  r <- make_clean_raster()
  # Set 50% of pixels to NA
  v <- terra::values(r)
  v[sample(nrow(v), round(nrow(v) * 0.5)), ] <- NA
  terra::values(r) <- v
  r
}

test_that("flight_qc() returns FlightQC object from SpatRaster", {
  r  <- make_clean_raster()
  qc <- flight_qc(r, flight_date = "2024-04-10", expected_bands = 4)
  expect_s3_class(qc, "FlightQC")
  expect_true(is.data.frame(qc$summary))
  expect_equal(nrow(qc$summary), 1)
})

test_that("flight_qc() summary has all required columns", {
  r  <- make_clean_raster()
  qc <- flight_qc(r)
  expected_cols <- c("label", "n_bands", "nodata_pct", "mean_brightness",
                     "irradiance_cv", "blur_laplacian", "flag_blur",
                     "flag_nodata", "flag_irradiance", "flag_saturation",
                     "overall_score", "recommendation")
  expect_true(all(expected_cols %in% names(qc$summary)))
})

test_that("flight_qc() recommendation is one of Accept/Caution/Reject", {
  r  <- make_clean_raster()
  qc <- flight_qc(r)
  expect_true(qc$summary$recommendation %in% c("Accept", "Caution", "Reject"))
})

test_that("flight_qc() flags blurry raster", {
  r  <- make_blurry_raster()
  qc <- flight_qc(r, blur_threshold = 100)
  expect_true(qc$summary$flag_blur)
})

test_that("flight_qc() flags high no-data fraction", {
  r  <- make_nodata_raster()
  qc <- flight_qc(r, nodata_threshold = 0.10)
  expect_true(qc$summary$flag_nodata)
})

test_that("flight_qc() flags wrong band count", {
  r  <- make_clean_raster()  # 4 bands
  qc <- flight_qc(r, expected_bands = 10)
  expect_true(qc$summary$flag_bands)
})

test_that("flight_qc() processes FlightStack with multiple flights", {
  # Build a small FlightStack
  wl <- seq(450, 850, length.out = 4)
  r1 <- make_clean_raster(); r2 <- make_clean_raster()

  fs <- flight_stack(
    rasters     = list(r1, r2),
    dates       = c("2024-04-01", "2024-05-01"),
    wavelengths = wl
  )
  qc <- flight_qc(fs)
  expect_s3_class(qc, "FlightQC")
  expect_equal(nrow(qc$summary), 2)
})

test_that("flight_qc() overall_score is between 0 and 100", {
  r  <- make_clean_raster()
  qc <- flight_qc(r)
  expect_gte(qc$summary$overall_score, 0)
  expect_lte(qc$summary$overall_score, 100)
})

test_that("flight_qc_report() creates an HTML file", {
  r   <- make_clean_raster()
  qc  <- flight_qc(r)
  tmp <- tempfile(fileext = ".html")
  out <- flight_qc_report(qc, output_file = tmp, open = FALSE)
  expect_true(file.exists(tmp))
  expect_equal(out, tmp)
  content <- readLines(tmp, warn = FALSE)
  expect_true(any(grepl("hovR", content)))
  expect_true(any(grepl("<table", content)))
  unlink(tmp)
})

test_that("flight_qc_report() HTML contains accept/reject keywords", {
  r   <- make_clean_raster()
  qc  <- flight_qc(r)
  tmp <- tempfile(fileext = ".html")
  flight_qc_report(qc, output_file = tmp, open = FALSE)
  content <- paste(readLines(tmp, warn = FALSE), collapse = " ")
  expect_true(grepl("Accept|Caution|Reject", content))
  unlink(tmp)
})
