library(testthat)
library(terra)

# ---- Helpers: synthetic FlightStack -----------------------------------------

make_stack <- function(n_flights = 4, nrow = 20, ncol = 20, n_bands = 10) {
  wl <- seq(450, 900, length.out = n_bands)
  dates <- as.Date("2024-03-01") + seq(0, (n_flights - 1) * 21, by = 21)

  rasters <- lapply(seq_len(n_flights), function(i) {
    # Simulate increasing greenness over season then decline
    phase <- sin(pi * i / n_flights)
    r <- terra::rast(nrows = nrow, ncols = ncol, nlyr = n_bands)
    terra::ext(r) <- c(0, 1, 0, 1)
    terra::crs(r) <- "EPSG:4326"
    # Set values: NIR high, Red moderate, with seasonal variation
    for (b in seq_len(n_bands)) {
      base_val <- if (wl[b] > 700) 0.5 + 0.3 * phase else 0.1 + 0.05 * phase
      terra::values(r[[b]]) <- base_val + rnorm(nrow * ncol, 0, 0.02)
    }
    r
  })

  flight_stack(
    rasters     = rasters,
    dates       = dates,
    wavelengths = wl,
    bbch        = c("21", "41", "65", "85")
  )
}

# ---- FlightStack tests -------------------------------------------------------

test_that("flight_stack() creates object with correct structure", {
  fs <- make_stack()

  expect_s3_class(fs, "FlightStack")
  expect_equal(length(fs$rasters), 4)
  expect_equal(length(fs$dates), 4)
  expect_equal(length(fs$wavelengths), 10)
  expect_equal(fs$bbch, c("21", "41", "65", "85"))
})

test_that("flight_stack() sorts dates chronologically", {
  wl <- seq(450, 900, length.out = 5)
  r  <- terra::rast(nrows = 10, ncols = 10, nlyr = 5)
  terra::ext(r) <- c(0, 1, 0, 1)
  terra::crs(r) <- "EPSG:4326"

  # Provide out-of-order dates
  fs <- flight_stack(
    rasters     = list(r, r, r),
    dates       = c("2024-06-01", "2024-03-01", "2024-09-01"),
    wavelengths = wl
  )
  expect_true(all(diff(as.integer(fs$dates)) > 0))
})

test_that("flight_stack() validates input lengths", {
  wl <- seq(450, 900, length.out = 5)
  r  <- terra::rast(nrows = 5, ncols = 5, nlyr = 5)
  terra::ext(r) <- c(0, 1, 0, 1)
  terra::crs(r) <- "EPSG:4326"

  expect_error(
    flight_stack(rasters = list(r, r), dates = "2024-01-01", wavelengths = wl),
    regexp = "length"
  )
  expect_error(
    flight_stack(rasters = list(r), dates = "2024-01-01", wavelengths = 1:3),
    regexp = "wavelengths"
  )
})

test_that("[ subsetting works", {
  fs <- make_stack(n_flights = 4)
  sub <- fs[1:2]
  expect_equal(length(sub), 2)
  expect_equal(sub$dates, fs$dates[1:2])
})

# ---- compute_vi tests -------------------------------------------------------

test_that("compute_vi() returns VIStack with correct structure", {
  fs <- make_stack()
  vi <- compute_vi(fs, indices = c("NDVI", "NDRE"))

  expect_s3_class(vi, "VIStack")
  expect_named(vi, c("NDVI", "NDRE"))
  expect_equal(length(vi$NDVI), 4)   # 4 flight dates
})

test_that("NDVI values are in [-1, 1] after clamping", {
  fs   <- make_stack()
  vi   <- compute_vi(fs, indices = "NDVI", clamp = TRUE)
  vals <- unlist(lapply(vi$NDVI, function(r) terra::values(r, na.rm = TRUE)))
  expect_true(all(vals >= -1 & vals <= 1))
})

test_that("list_vi() returns a data frame with expected columns", {
  df <- list_vi()
  expect_s3_class(df, "data.frame")
  expect_true(all(c("name", "description", "wavelengths", "tolerance") %in% names(df)))
  expect_gte(nrow(df), 10)
})

test_that("custom index is accepted and computed", {
  fs <- make_stack()
  custom <- list(
    bands = c(NIR = 800, Red = 670),
    tol   = 30,
    fun   = function(b) (b$NIR - b$Red) / (b$NIR + b$Red + 0.5),  # SAVI-like
    desc  = "Custom SAVI-like index"
  )
  vi <- compute_vi(fs, indices = "CUSTOM", custom_index = custom)
  expect_true("CUSTOM" %in% names(vi))
})

test_that("compute_vi() errors on unknown index name", {
  fs <- make_stack()
  expect_error(compute_vi(fs, indices = "FAKE_INDEX"), regexp = "Unknown")
})

# ---- integrate_season tests -------------------------------------------------

test_that("integrate_season() returns a single-band SpatRaster", {
  fs  <- make_stack()
  vi  <- compute_vi(fs, indices = "NDVI")
  auc <- integrate_season(vi, index = "NDVI")

  expect_s4_class(auc, "SpatRaster")
  expect_equal(terra::nlyr(auc), 1)
})

test_that("integrate_season() produces positive values for always-positive VI", {
  fs  <- make_stack()
  vi  <- compute_vi(fs, indices = "NDVI")
  auc <- integrate_season(vi, index = "NDVI")
  # NDVI of healthy veg > 0, so AUC should be positive
  expect_true(mean(terra::values(auc, na.rm = TRUE)) > 0)
})

test_that("integrate_season() date window subsetting reduces AUC", {
  fs      <- make_stack(n_flights = 4)
  vi      <- compute_vi(fs, indices = "NDVI")
  auc_all <- integrate_season(vi, index = "NDVI")
  auc_sub <- integrate_season(vi, index = "NDVI",
                               date_from = fs$dates[2], date_to = fs$dates[3])

  mean_all <- mean(terra::values(auc_all, na.rm = TRUE))
  mean_sub <- mean(terra::values(auc_sub, na.rm = TRUE))
  expect_lt(mean_sub, mean_all)
})

test_that("integrate_season() step and trapezoid methods give different results", {
  fs   <- make_stack()
  vi   <- compute_vi(fs, indices = "NDVI")
  trap <- integrate_season(vi, index = "NDVI", method = "trapezoid")
  step <- integrate_season(vi, index = "NDVI", method = "step")

  trap_mean <- mean(terra::values(trap, na.rm = TRUE))
  step_mean <- mean(terra::values(step, na.rm = TRUE))
  expect_false(isTRUE(all.equal(trap_mean, step_mean)))
})

test_that("integrate_season() errors when fewer than 2 flights in window", {
  fs <- make_stack(n_flights = 2)
  vi <- compute_vi(fs, indices = "NDVI")
  expect_error(
    integrate_season(vi, index = "NDVI",
                     date_from = fs$dates[1], date_to = fs$dates[1]),
    regexp = "2 flights"
  )
})

# ---- vi_derivative tests ----------------------------------------------------

test_that("vi_derivative() returns n-1 rasters", {
  fs   <- make_stack(n_flights = 4)
  vi   <- compute_vi(fs, indices = "NDVI")
  deriv <- vi_derivative(vi, index = "NDVI")

  expect_equal(length(deriv), 3)  # 4 flights → 3 intervals
})

test_that("vi_derivative() names follow date_to_date format", {
  fs    <- make_stack(n_flights = 3)
  vi    <- compute_vi(fs, indices = "NDVI")
  deriv <- vi_derivative(vi, index = "NDVI")

  expect_match(names(deriv)[1], "\\d{4}-\\d{2}-\\d{2}_to_\\d{4}-\\d{2}-\\d{2}")
})

test_that("vi_derivative() errors with fewer than 2 flights", {
  fs <- make_stack(n_flights = 1)
  vi <- compute_vi(fs, indices = "NDVI")
  expect_error(vi_derivative(vi, index = "NDVI"), regexp = "2 flights")
})

# ---- phenostage_tag tests ---------------------------------------------------

test_that("phenostage_tag() populates bbch slot", {
  fs_no_bbch <- make_stack()
  fs_no_bbch$bbch <- NULL

  obs <- data.frame(
    date = fs_no_bbch$dates[c(1, 3)],
    bbch = c("21", "65")
  )
  tagged <- phenostage_tag(fs_no_bbch, bbch_obs = obs, method = "nearest")
  expect_equal(length(tagged$bbch), length(fs_no_bbch$dates))
  expect_false(any(is.na(tagged$bbch)))
})

test_that("phenostage_tag() linear interpolation produces intermediate stages", {
  fs  <- make_stack(n_flights = 4)
  obs <- data.frame(
    date = as.Date(c("2024-03-01", "2024-09-01")),
    bbch = c("10", "90")
  )
  tagged <- phenostage_tag(fs, bbch_obs = obs, method = "linear")
  bbch_num <- as.integer(tagged$bbch)
  # Middle flights should have intermediate BBCH values
  expect_true(bbch_num[2] > 10 && bbch_num[2] < 90)
  expect_true(bbch_num[3] > 10 && bbch_num[3] < 90)
})
