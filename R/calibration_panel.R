#' Automated calibration panel detection and reflectance calibration
#'
#' @description
#' Before each drone flight, grey calibration panels are placed on the ground.
#' Converting raw DN values to surface reflectance requires detecting those
#' panels in the image, reading their known reflectance values, and fitting a
#' calibration curve. This module automates the entire process.
#'
#' Currently this is done by hand - opening each flight's image, manually
#' clicking on the panels, noting the pixel values, and running a spreadsheet.
#' For a 10-flight season that is a full day of work. These functions
#' reduce it to a single function call per flight.
#'
#' @name calibration_panel
NULL

# ---------------------------------------------------------------------------
# detect_panels()
# ---------------------------------------------------------------------------

#' Detect calibration panels in a drone image
#'
#' @description
#' Locates grey calibration panels (e.g. Spectralon, MicaSense panels, or
#' DIY grey cards) in a drone orthomosaic or raw flight image. Uses brightness
#' thresholding and shape filtering to find candidate regions, then ranks
#' them by their spectral flatness (grey panels have near-constant reflectance
#' across bands).
#'
#' @param raster A \code{terra::SpatRaster}. An RGB orthomosaic (3 bands) or
#'   a hyperspectral raster with \code{wavelengths} provided.
#' @param wavelengths Numeric vector of band wavelengths (nm). Required for
#'   hyperspectral input.
#' @param n_panels Integer. Expected number of panels in the scene. The
#'   function returns exactly this many candidate regions, ranked by
#'   confidence. Default: \code{1}.
#' @param panel_size_m2 Numeric vector of length 1 or 2. Expected panel area
#'   in square metres. If length 2, treated as (min, max). Default: \code{c(0.1, 2)}.
#' @param brightness_quantile Numeric 0-1. Panels are expected to be among the
#'   brightest objects in the scene. Pixels below this brightness quantile are
#'   excluded before panel search. Default: \code{0.92}.
#' @param spectral_flatness_threshold Numeric. Maximum allowed coefficient of
#'   variation (CV) across bands within a candidate panel region. Lower = stricter
#'   grey requirement. Default: \code{0.08} (8% CV).
#' @param search_area An optional \code{sf} polygon restricting the search area.
#'   Useful when panels are always placed in a corner of the field.
#'
#' @return An \code{sf} data frame with one row per detected panel:
#'   \describe{
#'     \item{\code{panel_id}}{Character label (\code{"panel_1"} etc.)}
#'     \item{\code{mean_dn}}{Per-band mean DN values (as a list-column)}
#'     \item{\code{brightness}}{Mean brightness across all bands}
#'     \item{\code{spectral_cv}}{Coefficient of variation across bands (grey flatness)}
#'     \item{\code{area_m2}}{Detected panel area}
#'     \item{\code{confidence}}{Numeric 0-1: detection confidence}
#'     \item{\code{geometry}}{Panel polygon}
#'   }
#'
#' @examples
#' \dontrun{
#' ortho   <- terra::rast("flight_01_ortho.tif")
#' panels  <- detect_panels(ortho, n_panels = 3, panel_size_m2 = c(0.25, 0.25))
#' print(panels[, c("panel_id", "brightness", "spectral_cv", "confidence")])
#' }
#'
#' @seealso \code{\link{calibrate_reflectance}}, \code{\link{panel_qc_plot}}
#' @export
detect_panels <- function(raster,
                          wavelengths              = NULL,
                          n_panels                 = 1L,
                          panel_size_m2            = c(0.1, 2),
                          brightness_quantile      = 0.92,
                          spectral_flatness_threshold = 0.08,
                          search_area              = NULL) {

  if (!inherits(raster, "SpatRaster"))
    cli::cli_abort("{.arg raster} must be a terra::SpatRaster.")

  if (length(panel_size_m2) == 1)
    panel_size_m2 <- c(panel_size_m2 * 0.5, panel_size_m2 * 1.5)

  cli::cli_h1("Calibration Panel Detection")

  # Restrict to search area if provided
  if (!is.null(search_area)) {
    raster <- terra::mask(terra::crop(raster, terra::vect(search_area)),
                          terra::vect(search_area))
  }

  n_bands <- terra::nlyr(raster)

  # Step 1: Compute brightness (mean across bands)
  cli::cli_inform("Step 1/4: Computing brightness map...")
  brightness <- terra::app(raster, fun = mean, na.rm = TRUE)

  # Step 2: Threshold to bright candidate pixels
  cli::cli_inform("Step 2/4: Isolating bright candidate pixels...")
  bq  <- terra::global(brightness, fun = function(x) {
    stats::quantile(x, brightness_quantile, na.rm = TRUE)
  })[[1]]
  bright_mask <- terra::ifel(brightness >= bq, 1L, NA)

  # Step 3: Clump connected bright regions -> candidate blobs
  cli::cli_inform("Step 3/4: Clustering candidate regions...")
  clumps <- terra::patches(bright_mask, directions = 8, allowGaps = FALSE)

  # Convert clumps to polygons
  clump_polys <- terra::as.polygons(clumps, dissolve = FALSE)
  clump_sf    <- sf::st_as_sf(clump_polys)
  clump_sf    <- clump_sf[sf::st_is_valid(clump_sf), ]
  clump_sf$area_m2 <- as.numeric(sf::st_area(clump_sf))

  # Filter by size
  clump_sf <- clump_sf[
    clump_sf$area_m2 >= panel_size_m2[1] &
    clump_sf$area_m2 <= panel_size_m2[2], ]

  if (nrow(clump_sf) == 0)
    cli::cli_abort(c(
      "No candidate regions found in expected size range \\
       ({panel_size_m2[1]} - {panel_size_m2[2]} m2).",
      "i" = "Try adjusting {.arg panel_size_m2} or {.arg brightness_quantile}."
    ))

  # Step 4: Score by spectral flatness (grey panels are spectrally flat)
  cli::cli_inform("Step 4/4: Scoring {nrow(clump_sf)} candidates by spectral flatness...")

  scores <- lapply(seq_len(nrow(clump_sf)), function(i) {
    poly  <- clump_sf[i, ]
    extr  <- terra::extract(raster, terra::vect(poly),
                             fun = mean, na.rm = TRUE)
    means <- as.numeric(extr[1, -1])

    if (all(is.na(means)) || length(means) == 0)
      return(list(mean_dn = means, brightness = NA_real_,
                  spectral_cv = Inf, confidence = 0))

    brightness_val <- mean(means, na.rm = TRUE)
    cv_val <- if (brightness_val > 0) {
      stats::sd(means, na.rm = TRUE) / brightness_val
    } else Inf

    conf <- max(0, 1 - cv_val / spectral_flatness_threshold)
    list(mean_dn     = means,
         brightness  = brightness_val,
         spectral_cv = round(cv_val, 4),
         confidence  = round(min(conf, 1), 3))
  })

  clump_sf$mean_dn    <- lapply(scores, `[[`, "mean_dn")
  clump_sf$brightness <- vapply(scores, `[[`, "brightness",   FUN.VALUE = numeric(1))
  clump_sf$spectral_cv<- vapply(scores, `[[`, "spectral_cv",  FUN.VALUE = numeric(1))
  clump_sf$confidence <- vapply(scores, `[[`, "confidence",   FUN.VALUE = numeric(1))

  # Keep only flat-spectrum candidates and rank
  candidates <- clump_sf[clump_sf$spectral_cv <= spectral_flatness_threshold, ]
  candidates <- candidates[order(-candidates$confidence), ]

  if (nrow(candidates) == 0) {
    cli::cli_warn(c(
      "No candidates passed spectral flatness threshold ({spectral_flatness_threshold}).",
      "i" = "Returning top {n_panels} bright region(s) regardless of flatness."
    ))
    candidates <- clump_sf[order(-clump_sf$confidence), ]
  }

  out <- head(candidates, n_panels)
  out$panel_id <- paste0("panel_", seq_len(nrow(out)))
  drop_col <- names(out)[1]   # clump ID column from terra
  out <- out[, setdiff(names(out), drop_col)]
  keep <- c("panel_id", "mean_dn", "brightness", "spectral_cv",
            "area_m2", "confidence", "geometry")
  out <- out[, intersect(keep, names(out))]

  cli::cli_inform(c(
    "v" = "{nrow(out)} panel(s) detected.",
    "i" = "Run {.fn calibrate_reflectance} to apply the calibration curve."
  ))
  out
}

# ---------------------------------------------------------------------------
# calibrate_reflectance()
# ---------------------------------------------------------------------------

#' Apply empirical line reflectance calibration using detected panels
#'
#' @description
#' Fits a per-band linear calibration from DN to reflectance using the
#' detected panel(s) as ground-control points. With a single panel, a
#' gain-only (through-origin) calibration is applied. With two or more
#' panels of different known reflectances, a full gain+offset (empirical
#' line) fit is used.
#'
#' @param raster A \code{terra::SpatRaster} to calibrate (the raw DN image).
#' @param panels An \code{sf} object from \code{\link{detect_panels}}.
#' @param known_reflectance Named numeric vector of known panel reflectance
#'   values (0-1 fraction). Names must match \code{panel_id} values in
#'   \code{panels}. Example: \code{c(panel_1 = 0.50, panel_2 = 0.20)}.
#'   If a single unnamed value is provided, it is applied to \code{panel_1}.
#' @param per_band Logical. If \code{TRUE} (default), fits separate calibration
#'   coefficients for each spectral band. If \code{FALSE}, fits a single
#'   global gain/offset across all bands.
#' @param clamp_output Logical. If \code{TRUE} (default), clamps output
#'   reflectance to \code{[0, 1]}.
#'
#' @return A \code{terra::SpatRaster} with values in units of reflectance (0-1).
#'   Calibration coefficients are stored as raster metadata.
#'
#' @examples
#' \dontrun{
#' panels <- detect_panels(ortho, n_panels = 2)
#' refl   <- calibrate_reflectance(
#'   raster           = ortho,
#'   panels           = panels,
#'   known_reflectance = c(panel_1 = 0.50, panel_2 = 0.20)
#' )
#' terra::plot(refl[[4]], main = "Band 4 reflectance")
#' }
#'
#' @export
calibrate_reflectance <- function(raster,
                                  panels,
                                  known_reflectance,
                                  per_band     = TRUE,
                                  clamp_output = TRUE) {

  if (!inherits(raster, "SpatRaster"))
    cli::cli_abort("{.arg raster} must be a terra::SpatRaster.")
  if (!inherits(panels, "sf"))
    cli::cli_abort("{.arg panels} must be an sf object from detect_panels().")

  # Resolve known reflectance
  if (is.null(names(known_reflectance)))
    known_reflectance <- stats::setNames(known_reflectance,
                                         paste0("panel_", seq_along(known_reflectance)))

  common_ids <- intersect(panels$panel_id, names(known_reflectance))
  if (length(common_ids) == 0)
    cli::cli_abort(
      "No matching panel IDs between {.arg panels} ({panels$panel_id}) \\
       and {.arg known_reflectance} ({names(known_reflectance)}).")

  cli::cli_h1("Reflectance Calibration")
  cli::cli_inform("Using {length(common_ids)} panel(s): {common_ids}")

  n_bands <- terra::nlyr(raster)
  method  <- if (length(common_ids) >= 2) "empirical_line" else "gain_only"
  cli::cli_inform("Calibration method: {method}")

  # Extract per-panel mean DN for each band
  panel_dns <- lapply(common_ids, function(pid) {
    p_row <- panels[panels$panel_id == pid, ]
    extr  <- terra::extract(raster, terra::vect(p_row), fun = mean, na.rm = TRUE)
    as.numeric(extr[1, -1])
  })
  names(panel_dns) <- common_ids

  known_refl_matched <- known_reflectance[common_ids]

  # Fit calibration per band (or globally)
  band_range <- if (per_band) seq_len(n_bands) else 1L

  gains   <- numeric(n_bands)
  offsets <- numeric(n_bands)

  for (b in band_range) {
    dns   <- vapply(common_ids, function(pid) panel_dns[[pid]][b], numeric(1))
    refls <- known_refl_matched

    if (method == "gain_only") {
      g <- refls[[1]] / max(dns[[1]], 1e-9)
      gains[b]   <- g
      offsets[b] <- 0
    } else {
      # Linear fit: reflectance = gain * DN + offset
      fit       <- stats::lm(refls ~ dns)
      gains[b]  <- stats::coef(fit)[["dns"]]
      offsets[b]<- stats::coef(fit)[["(Intercept)"]]
    }

    if (!per_band) {
      gains   <- rep(gains[1],   n_bands)
      offsets <- rep(offsets[1], n_bands)
      break
    }
  }

  cli::cli_inform(
    "Gain range: [{round(min(gains),4)}, {round(max(gains),4)}]  \\
     Offset range: [{round(min(offsets),4)}, {round(max(offsets),4)}]"
  )

  # Apply calibration band-by-band
  calib <- terra::rast(lapply(seq_len(n_bands), function(b) {
    band <- raster[[b]] * gains[b] + offsets[b]
    if (clamp_output) band <- terra::clamp(band, lower = 0, upper = 1)
    band
  }))
  names(calib) <- names(raster)

  # Store coefficients as metadata
  attr(calib, "hovR_meta") <- list(
    calibration_method  = method,
    n_panels            = as.character(length(common_ids)),
    panel_ids           = paste(common_ids, collapse = ","),
    known_reflectances  = paste(round(known_refl_matched, 3), collapse = ",")
  )

  cli::cli_inform(c("v" = "Calibration applied. Output in reflectance units (0-1)."))
  calib
}

# ---------------------------------------------------------------------------
# panel_qc_plot()
# ---------------------------------------------------------------------------

#' Visual QC for panel detection
#'
#' @description
#' Plots detected panels overlaid on a greyscale brightness map of the image,
#' and shows per-panel spectral profiles vs the expected flat line.
#'
#' @param raster A \code{terra::SpatRaster}.
#' @param panels An \code{sf} object from \code{\link{detect_panels}}.
#' @param wavelengths Numeric vector of band wavelengths (nm). If \code{NULL},
#'   band indices are used on the x-axis.
#'
#' @return A named list of two \code{ggplot} objects:
#'   \code{$map} (spatial panel locations) and
#'   \code{$spectra} (spectral profiles of detected panels).
#'
#' @export
panel_qc_plot <- function(raster, panels, wavelengths = NULL) {

  n_bands <- terra::nlyr(raster)
  x_vals  <- if (!is.null(wavelengths) && length(wavelengths) == n_bands)
               wavelengths else seq_len(n_bands)
  x_lab   <- if (!is.null(wavelengths)) "Wavelength (nm)" else "Band index"

  # Brightness map
  bright <- terra::app(raster, mean, na.rm = TRUE)
  bdf    <- as.data.frame(bright, xy = TRUE, na.rm = TRUE)
  names(bdf)[3] <- "brightness"

  map_plot <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = bdf,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$brightness)) +
    ggplot2::scale_fill_gradient(low = "black", high = "white",
                                  name = "Brightness") +
    ggplot2::geom_sf(data = panels, fill = NA,
                     colour = "#FF4500", linewidth = 1) +
    ggplot2::geom_sf_label(data = panels,
      ggplot2::aes(label = .data$panel_id),
      size = 3, colour = "#FF4500", fill = "white", alpha = 0.8) +
    ggplot2::labs(title = "Detected calibration panels", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::coord_sf()

  # Spectral profiles
  spec_rows <- lapply(seq_len(nrow(panels)), function(i) {
    dns <- panels$mean_dn[[i]]
    if (is.null(dns)) return(NULL)
    data.frame(panel_id  = panels$panel_id[i],
               x         = x_vals[seq_along(dns)],
               mean_dn   = dns)
  })
  spec_df <- do.call(rbind, Filter(Negate(is.null), spec_rows))

  spec_plot <- ggplot2::ggplot(spec_df,
    ggplot2::aes(x = .data$x, y = .data$mean_dn,
                 colour = .data$panel_id, group = .data$panel_id)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(title  = "Panel spectral profiles (should be flat)",
                  x      = x_lab,
                  y      = "Mean DN",
                  colour = "Panel") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  list(map = map_plot, spectra = spec_plot)
}
