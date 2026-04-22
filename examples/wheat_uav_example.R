# =========================================================
# hovR Example: UAV Multispectral Time-Series Analysis
# =========================================================

# ---- Install & load -----------------------------------------------------
# install.packages(c("terra", "ggplot2"))
# devtools::install_github("Lalitgis/hovR", upgrade = "never")

library(hovR)
library(terra)
library(ggplot2)

# ---- 1. Load data -------------------------------------------------------

# Folder containing UAV images
img_dir <- "C:/Users/Hp/Downloads/UAV-based multispectral images of spring wheat/images/PC"

# List all .tif files
files <- sort(list.files(img_dir, pattern = "\\.tif$", full.names = TRUE))

# Extract dates from filenames (YYYYMMDD format)
dates <- as.Date(substr(basename(files), 1, 8), format = "%Y%m%d")

# ---- 2. Read rasters (keep 4 bands only) --------------------------------
# Parrot Sequoia bands: Green, Red, RedEdge, NIR

band_names  <- c("Green", "Red", "RedEdge", "NIR")
wavelengths <- c(550, 660, 735, 790)

rasters <- lapply(files, function(f) {
  r <- rast(f)[[1:4]]       # drop band 5 (empty)
  names(r) <- band_names
  r
})

# ---- 3. Build flight stack ----------------------------------------------

stack <- flight_stack(
  rasters     = rasters,
  dates       = dates,
  wavelengths = wavelengths
)

print(stack)
# Expected: Flights ~19 | Bands: 4

# ---- 4. Quality control -------------------------------------------------

qc <- flight_qc(stack, expected_bands = 4)
print(qc)

# ---- 5. Compute vegetation indices -------------------------------------

vi <- compute_vi(
  stack,
  indices = c("NDVI", "NDRE", "GNDVI")
)

# ---- 6. Visualize NDVI (selected dates) --------------------------------

sample_idx <- c(2, 10, 17)

par(mfrow = c(1, 3), mar = c(3, 3, 3, 2))

for (i in sample_idx) {
  plot(
    vi$NDVI[[i]],
    main  = paste("NDVI:", dates[i]),
    col   = hcl.colors(50, "RdYlGn"),
    range = c(0, 0.9)
  )
}

par(mfrow = c(1, 1))

# ---- 7. Time-series (field mean) ---------------------------------------

mean_index <- function(x) mean(values(x), na.rm = TRUE)

ts_data <- data.frame(
  date  = dates,
  NDVI  = sapply(vi$NDVI,  mean_index),
  NDRE  = sapply(vi$NDRE,  mean_index),
  GNDVI = sapply(vi$GNDVI, mean_index)
)

# Convert to long format
ts_long <- data.frame(
  date  = rep(ts_data$date, 3),
  value = c(ts_data$NDVI, ts_data$NDRE, ts_data$GNDVI),
  index = rep(c("NDVI", "NDRE", "GNDVI"), each = nrow(ts_data))
)

# Plot
ggplot(ts_long, aes(date, value, colour = index)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  labs(
    title = "Spring Wheat Growth Curve (2018)",
    x = "Date",
    y = "Mean Vegetation Index"
  ) +
  theme_minimal()

# ---- 8. Seasonal integration -------------------------------------------

season_ndvi <- integrate_season(
  vi,
  index  = "NDVI",
  method = "trapezoid"
)

plot(
  season_ndvi,
  main = "Seasonal NDVI (Productivity Map)",
  col  = hcl.colors(50, "Greens")
)

# Save output
writeRaster(
  season_ndvi,
  "wheat_seasonal_ndvi.tif",
  overwrite = TRUE
)

# ---- 9. Growth rate (derivative) ---------------------------------------

rates <- vi_derivative(
  vi,
  index  = "NDVI",
  smooth = TRUE,
  units  = "days"
)

mean_rates <- sapply(rates, mean_index)

rate_table <- data.frame(
  interval = names(mean_rates),
  rate     = mean_rates
)

print(rate_table)

# Identify key periods
rate_table[which.max(mean_rates), ]  # fastest growth
rate_table[which.min(mean_rates), ]  # strongest decline