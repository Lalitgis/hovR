# =========================================================
# hovR Example: UAV Multispectral Time-Series Analysis
# Parrot Sequoia 4-band sensor (Green, Red, RedEdge, NIR)
# =========================================================

library(hovR)
library(terra)
library(ggplot2)

# ---- 1. Load data -------------------------------------------------------
# Replace img_dir with the folder containing your .tif flight images.
# Files should be named with the acquisition date as YYYYMMDD_*.tif
# (e.g. "20180410_sequoia.tif").

img_dir <- "path/to/your/flight/images"

files <- sort(list.files(img_dir, pattern = "\\.tif$", full.names = TRUE))
dates <- as.Date(substr(basename(files), 1, 8), format = "%Y%m%d")

# Parrot Sequoia band order: Green, Red, RedEdge, NIR
band_names  <- c("Green", "Red", "RedEdge", "NIR")
wavelengths <- c(550, 660, 735, 790)

rasters <- lapply(files, function(f) {
  r <- rast(f)[[1:4]]   # keep first 4 bands only
  names(r) <- band_names
  r
})

# ---- 2. Build flight stack ----------------------------------------------
stack <- flight_stack(
  rasters     = rasters,
  dates       = dates,
  wavelengths = wavelengths
)
print(stack)

# ---- 3. Quality control -------------------------------------------------
qc <- flight_qc(stack, expected_bands = 4)
print(qc)
flight_qc_report(qc, output_file = "qc_season.html", open = FALSE)

# ---- 4. Calibration panels (if available) --------------------------------
# panels <- detect_panels(rasters[[1]], wavelengths = wavelengths, n_panels = 2)
# rasters <- lapply(rasters, calibrate_reflectance,
#                   panels = panels,
#                   known_reflectance = c(panel_1 = 0.50, panel_2 = 0.20))

# ---- 5. Compute vegetation indices --------------------------------------
vi <- compute_vi(stack, indices = c("NDVI", "NDRE", "GNDVI"))

# ---- 6. Time-series (field mean) ----------------------------------------
mean_vi <- function(x) mean(terra::values(x), na.rm = TRUE)

ts_data <- data.frame(
  date  = dates,
  NDVI  = sapply(vi$NDVI,  mean_vi),
  NDRE  = sapply(vi$NDRE,  mean_vi),
  GNDVI = sapply(vi$GNDVI, mean_vi)
)

ts_long <- reshape(ts_data,
  varying   = c("NDVI", "NDRE", "GNDVI"),
  v.names   = "value",
  timevar   = "index",
  times     = c("NDVI", "NDRE", "GNDVI"),
  direction = "long"
)

ggplot(ts_long, aes(date, value, colour = index)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  labs(title = "Spring Wheat Growth Curve", x = "Date", y = "Mean VI") +
  theme_minimal()

# ---- 7. Seasonal integration (AUC) --------------------------------------
season_ndvi <- integrate_season(vi, index = "NDVI", method = "trapezoid")
terra::plot(season_ndvi, main = "Seasonal NDVI (Productivity Map)",
            col = hcl.colors(50, "Greens"))
terra::writeRaster(season_ndvi, "wheat_seasonal_ndvi.tif", overwrite = TRUE)

# ---- 8. Growth rate (derivative) ----------------------------------------
rates      <- vi_derivative(vi, index = "NDVI", smooth = TRUE, units = "days")
mean_rates <- sapply(rates, mean_vi)

rate_table <- data.frame(interval = names(mean_rates), rate = mean_rates,
                          row.names = NULL)
print(rate_table)
rate_table[which.max(mean_rates), ]   # fastest green-up
rate_table[which.min(mean_rates), ]   # peak senescence
