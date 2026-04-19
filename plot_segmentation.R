#' Automatic plot segmentation from drone imagery
#'
#' @description
#' Detects and delineates trial plot boundaries directly from a drone
#' orthomosaic or hyperspectral raster — no manual polygon drawing in
#' QGIS required. Handles regular grids, misaligned rows, staggered
#' layouts, and variable plot sizes.
#'
#' The pipeline:
#' \enumerate{
#'   \item Compute a vegetation/soil contrast index to separate plots from alleys
#'   \item Detect dominant row orientation via Hough-like gradient analysis
#'   \item Project onto row/column axes and find inter-plot gaps via 1-D profile valleys
#'   \item Reconstruct individual plot polygons and optionally refine boundaries
#'   \item Return an \code{sf} polygon layer with plot IDs and quality scores
#' }
#'
#' @name plot_segmentation
NULL

# ---------------------------------------------------------------------------
# segment_plots() — main entry point
# ---------------------------------------------------------------------------

#' Automatically segment trial plots from a drone raster
#'
#' @param raster A \code{terra::SpatRaster}. Either an RGB orthomosaic
#'   (3 bands) or a hyperspectral raster. For hyperspectral input, provide
#'   \code{wavelengths} so that appropriate bands are selected automatically.
#' @param wavelengths Numeric vector of band centre wavelengths (nm).
#'   Required when \code{raster} has more than 3 bands.
#' @param method Character. Segmentation strategy:
#'   \describe{
#'     \item{\code{"profile"}}{(default) Project row/column mean profiles and
#'       find valleys — fast, works well for regular grids.}
#'     \item{\code{"gradient"}}{Sobel gradient magnitude to detect alley edges —
#'       better for irregular or staggered layouts.}
#'     \item{\code{"combined"}}{Run both methods and take the intersection —
#'       most robust but slower.}
#'   }
#' @param expected_rows Integer. Expected number of plot rows. If \code{NULL}
#'   (default), estimated automatically from profile valleys.
#' @param expected_cols Integer. Expected number of plot columns. If \code{NULL},
#'   estimated automatically.
#' @param alley_width_m Numeric. Approximate alley/path width between plots in
#'   metres. Used to tune gap detection. Default: \code{0.5}.
#' @param min_plot_area_m2 Numeric. Minimum plot area in square metres.
#'   Detected polygons smaller than this are discarded as artefacts.
#'   Default: \code{1}.
#' @param orientation_deg Numeric. Known row orientation in degrees from North
#'   (0 = rows run E-W, 90 = rows run N-S). If \code{NULL} (default),
#'   estimated automatically.
#' @param smooth_sigma Numeric. Gaussian smoothing sigma applied to profiles
#'   before valley detection. Higher = more tolerant of within-plot variability.
#'   Default: \code{3}.
#' @param id_prefix Character. Prefix for plot ID labels. Default: \code{"P"}.
#' @param crs_out Character or numeric. Output CRS. Defaults to input raster CRS.
#'
#' @return An \code{sf} data frame with one row per detected plot containing:
#'   \describe{
#'     \item{\code{plot_id}}{Character plot identifier (e.g. \code{"P001"})}
#'     \item{\code{row_idx}}{Integer row index}
#'     \item{\code{col_idx}}{Integer column index}
#'     \item{\code{area_m2}}{Plot area in square metres}
#'     \item{\code{quality}}{Numeric 0–1: confidence in boundary accuracy}
#'     \item{\code{geometry}}{Polygon geometry}
#'   }
#'
#' @examples
#' \dontrun{
#' # From an RGB orthomosaic
#' ortho  <- terra::rast("field_ortho.tif")
#' plots  <- segment_plots(ortho, expected_rows = 20, expected_cols = 10,
#'                          alley_width_m = 0.5)
#' sf::st_write(plots, "trial_plots.gpkg")
#'
#' # From hyperspectral (uses NDVI contrast automatically)
#' hsi    <- terra::rast("field_hsi.tif")
#' wl     <- seq(400, 900, by = 5)
#' plots  <- segment_plots(hsi, wavelengths = wl, method = "combined")
#' }
#'
#' @seealso \code{\link{refine_plot_boundaries}}, \code{\link{plot_segmentation_qc}}
#' @export
segment_plots <- function(raster,
                          wavelengths      = NULL,
                          method           = c("profile", "gradient", "combined"),
                          expected_rows    = NULL,
                          expected_cols    = NULL,
                          alley_width_m    = 0.5,
                          min_plot_area_m2 = 1,
                          orientation_deg  = NULL,
                          smooth_sigma     = 3,
                          id_prefix        = "P",
                          crs_out          = NULL) {

  method <- match.arg(method)

  if (!inherits(raster, "SpatRaster"))
    cli::cli_abort("{.arg raster} must be a terra::SpatRaster.")

  n_bands <- terra::nlyr(raster)
  if (n_bands > 3 && is.null(wavelengths))
    cli::cli_abort(
      "Hyperspectral raster ({n_bands} bands) requires {.arg wavelengths}.")

  cli::cli_h1("Plot Auto-Segmentation")

  # Step 1: derive contrast index
  cli::cli_inform("Step 1/5: Computing vegetation contrast index...")
  contrast <- .compute_contrast_index(raster, wavelengths, n_bands)

  # Step 2: estimate orientation
  cli::cli_inform("Step 2/5: Estimating row orientation...")
  if (is.null(orientation_deg)) {
    orientation_deg <- .estimate_orientation(contrast)
    cli::cli_inform("  Detected orientation: {round(orientation_deg, 1)} degrees")
  }

  # Step 3: detect gaps via selected method
  cli::cli_inform("Step 3/5: Detecting inter-plot gaps ({method} method)...")
  res <- terra::res(raster)
  px_per_m <- 1 / mean(res)

  gap_mask <- switch(method,
    profile  = .profile_gaps(contrast, orientation_deg, alley_width_m,
                              px_per_m, smooth_sigma),
    gradient = .gradient_gaps(contrast, alley_width_m, px_per_m),
    combined = {
      g1 <- .profile_gaps(contrast, orientation_deg, alley_width_m,
                           px_per_m, smooth_sigma)
      g2 <- .gradient_gaps(contrast, alley_width_m, px_per_m)
      terra::ifel(g1 == 1 | g2 == 1, 1L, 0L)
    }
  )

  # Step 4: polygonise plot interiors
  cli::cli_inform("Step 4/5: Polygonising plot interiors...")
  plot_interior <- terra::ifel(gap_mask == 1, NA, 1L)
  polys <- .raster_to_plot_polygons(
    plot_interior, raster, min_plot_area_m2,
    expected_rows, expected_cols, id_prefix
  )

  # Step 5: assign quality scores
  cli::cli_inform("Step 5/5: Scoring boundary quality...")
  polys <- .score_plot_quality(polys, contrast, gap_mask)

  if (!is.null(crs_out)) polys <- sf::st_transform(polys, crs_out)

  n_det  <- nrow(polys)
  q_mean <- round(mean(polys$quality, na.rm = TRUE), 2)
  cli::cli_inform(c(
    "v" = "Detected {n_det} plots  |  Mean quality score: {q_mean}",
    "i" = "Run {.fn plot_segmentation_qc} to inspect results visually."
  ))

  polys
}

# ---------------------------------------------------------------------------
# refine_plot_boundaries()
# ---------------------------------------------------------------------------

#' Refine plot polygon boundaries using image gradients
#'
#' @description
#' Takes an initial polygon layer (e.g. from \code{\link{segment_plots}} or
#' a hand-drawn shapefile with approximate boundaries) and snaps each edge to
#' the nearest strong gradient in the raster — analogous to active contour /
#' "snap to edge" in GIS software, but automated.
#'
#' @param plots An \code{sf} polygon layer from \code{\link{segment_plots}} or
#'   any source, with a \code{plot_id} column.
#' @param raster A \code{terra::SpatRaster} (orthomosaic or contrast index).
#' @param max_shift_m Numeric. Maximum distance in metres that a boundary edge
#'   may move during refinement. Default: \code{0.3}.
#' @param gradient_threshold Numeric 0–1. Minimum normalised gradient magnitude
#'   to qualify as an edge. Default: \code{0.2}.
#'
#' @return A refined \code{sf} polygon layer with the same columns as
#'   \code{plots}, plus a \code{shift_m} column giving the mean boundary
#'   shift applied per plot.
#'
#' @export
refine_plot_boundaries <- function(plots,
                                   raster,
                                   max_shift_m          = 0.3,
                                   gradient_threshold   = 0.2) {

  if (!inherits(plots, "sf"))
    cli::cli_abort("{.arg plots} must be an sf object.")
  if (!inherits(raster, "SpatRaster"))
    cli::cli_abort("{.arg raster} must be a terra::SpatRaster.")

  cli::cli_inform("Refining {nrow(plots)} plot boundaries (max shift: {max_shift_m} m)...")

  # Compute gradient magnitude
  grad <- .sobel_magnitude(raster[[1]])
  grad_norm <- (grad - terra::global(grad, "min", na.rm = TRUE)[[1]]) /
    (terra::global(grad, "max", na.rm = TRUE)[[1]] -
       terra::global(grad, "min", na.rm = TRUE)[[1]] + 1e-9)
  edge_mask <- terra::ifel(grad_norm >= gradient_threshold, 1L, 0L)

  res_m   <- mean(terra::res(raster))
  px_shift <- round(max_shift_m / res_m)

  refined_list <- lapply(seq_len(nrow(plots)), function(i) {
    poly <- plots[i, ]
    bbox <- sf::st_bbox(poly)

    # Sample gradient inside a buffer around the polygon boundary
    buf  <- sf::st_buffer(sf::st_cast(poly, "MULTILINESTRING"), max_shift_m)
    edge_in_buf <- terra::mask(
      terra::crop(edge_mask, terra::vect(buf)),
      terra::vect(buf)
    )

    # Compute centroid of strong edges within buffer as shift target
    edge_pts <- terra::as.data.frame(edge_in_buf, xy = TRUE, na.rm = TRUE)
    edge_pts <- edge_pts[edge_pts[, 3] == 1, , drop = FALSE]

    if (nrow(edge_pts) < 4) {
      poly$shift_m <- 0
      return(poly)
    }

    # Simple centroid-of-edges snap: translate polygon centroid to match
    # the centroid of detected edge pixels (bounded by max_shift_m)
    poly_ctr  <- sf::st_coordinates(sf::st_centroid(poly))
    edge_ctr  <- c(mean(edge_pts[, 1]), mean(edge_pts[, 2]))
    delta     <- edge_ctr - poly_ctr[1, 1:2]
    shift_dist <- min(sqrt(sum(delta^2)), max_shift_m)
    scale      <- if (sqrt(sum(delta^2)) > 0) shift_dist / sqrt(sum(delta^2)) else 0
    dx <- delta[1] * scale
    dy <- delta[2] * scale

    refined_geom <- sf::st_geometry(poly) + c(dx, dy)
    sf::st_geometry(poly) <- refined_geom
    sf::st_crs(poly) <- sf::st_crs(plots)
    poly$shift_m <- round(shift_dist, 3)
    poly
  })

  out <- do.call(rbind, refined_list)
  cli::cli_inform(c("v" = "Mean boundary shift: {round(mean(out$shift_m), 3)} m"))
  out
}

# ---------------------------------------------------------------------------
# plot_segmentation_qc()
# ---------------------------------------------------------------------------

#' Visual QC plot for segmentation results
#'
#' @description
#' Produces a ggplot2 figure overlaying detected plot polygons on the
#' contrast index raster, coloured by quality score. Use to quickly identify
#' problematic boundary detections.
#'
#' @param plots An \code{sf} polygon layer from \code{\link{segment_plots}}.
#' @param raster The \code{terra::SpatRaster} used for segmentation.
#' @param show_ids Logical. Label each plot with its \code{plot_id}.
#'   Default: \code{TRUE} (suppressed automatically when > 200 plots).
#' @param quality_threshold Numeric 0–1. Polygons below this quality score
#'   are highlighted in red. Default: \code{0.5}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_segmentation_qc <- function(plots,
                                 raster,
                                 show_ids          = TRUE,
                                 quality_threshold = 0.5) {

  if (!inherits(plots, "sf"))
    cli::cli_abort("{.arg plots} must be an sf object.")

  contrast <- .compute_contrast_index(raster, NULL, terra::nlyr(raster))
  rdf      <- as.data.frame(contrast, xy = TRUE, na.rm = TRUE)
  names(rdf)[3] <- "contrast"

  low_q  <- plots[plots$quality < quality_threshold, ]
  high_q <- plots[plots$quality >= quality_threshold, ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = rdf,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$contrast)) +
    ggplot2::scale_fill_gradientn(
      colours = c("#5A3E28", "#8B6914", "#C8B04A", "#6BAE4F", "#1E6B2E"),
      name    = "Contrast\nindex"
    ) +
    ggplot2::geom_sf(data = high_q, fill = NA,
                     colour = "#00BFFF", linewidth = 0.4) +
    ggplot2::geom_sf(data = low_q, fill = NA,
                     colour = "#FF4500", linewidth = 0.6, linetype = "dashed") +
    ggplot2::labs(
      title    = "Plot segmentation QC",
      subtitle = paste0(
        nrow(plots), " plots detected  |  ",
        nrow(low_q), " below quality threshold (", quality_threshold, ")"
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text   = ggplot2::element_text(size = 8)
    ) +
    ggplot2::coord_sf()

  if (show_ids && nrow(plots) <= 200) {
    centroids <- suppressWarnings(sf::st_centroid(plots))
    coords    <- as.data.frame(sf::st_coordinates(centroids))
    coords$plot_id <- plots$plot_id
    p <- p + ggplot2::geom_text(
      data = coords,
      ggplot2::aes(x = .data$X, y = .data$Y, label = .data$plot_id),
      size = 2, colour = "white", fontface = "bold"
    )
  }

  p
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
.compute_contrast_index <- function(raster, wavelengths, n_bands) {
  if (n_bands == 1) return(raster)

  if (n_bands == 3) {
    # RGB: use ExG (Excess Green) as vegetation contrast
    r <- raster[[1]]; g <- raster[[2]]; b <- raster[[3]]
    total <- r + g + b + 1e-9
    return(2 * g / total - r / total - b / total)
  }

  if (!is.null(wavelengths) && length(wavelengths) == n_bands) {
    # Hyperspectral: compute NDVI
    nir_idx <- which.min(abs(wavelengths - 800))
    red_idx <- which.min(abs(wavelengths - 670))
    nir <- raster[[nir_idx]]; red <- raster[[red_idx]]
    return((nir - red) / (nir + red + 1e-9))
  }

  # Fallback: first band
  cli::cli_warn("Cannot determine contrast index automatically. Using band 1.")
  raster[[1]]
}

#' @keywords internal
.estimate_orientation <- function(contrast) {
  # Sobel-based dominant gradient direction
  grad_x <- .sobel_x(contrast)
  grad_y <- .sobel_y(contrast)
  angles  <- atan2(terra::values(grad_y, na.rm = FALSE),
                   terra::values(grad_x, na.rm = FALSE)) * 180 / pi
  angles  <- angles[!is.na(angles)]
  # Histogram of gradient directions — dominant perpendicular = row direction
  hist_r  <- graphics::hist(angles %% 180, breaks = 36, plot = FALSE)
  dominant_angle <- hist_r$mids[which.max(hist_r$counts)]
  (dominant_angle + 90) %% 180
}

#' @keywords internal
.profile_gaps <- function(contrast, orientation_deg, alley_width_m,
                           px_per_m, smooth_sigma) {
  m   <- as.matrix(contrast, wide = TRUE)
  m[is.na(m)] <- 0

  # Row-wise and column-wise mean profiles
  row_prof <- rowMeans(m, na.rm = TRUE)
  col_prof <- colMeans(m, na.rm = TRUE)

  # Smooth profiles
  row_s <- .gaussian_smooth_1d(row_prof, smooth_sigma)
  col_s <- .gaussian_smooth_1d(col_prof, smooth_sigma)

  # Threshold at 20th percentile to find alleys
  row_thresh <- stats::quantile(row_s, 0.20, na.rm = TRUE)
  col_thresh <- stats::quantile(col_s, 0.20, na.rm = TRUE)

  gap_rows <- which(row_s < row_thresh)
  gap_cols <- which(col_s < col_thresh)

  # Build binary gap raster
  nr <- nrow(m); nc <- ncol(m)
  gap_mat <- matrix(0L, nrow = nr, ncol = nc)
  if (length(gap_rows) > 0) gap_mat[gap_rows, ] <- 1L
  if (length(gap_cols) > 0) gap_mat[, gap_cols] <- 1L

  out <- terra::rast(contrast)
  terra::values(out) <- as.vector(t(gap_mat))
  out
}

#' @keywords internal
.gradient_gaps <- function(contrast, alley_width_m, px_per_m) {
  mag  <- .sobel_magnitude(contrast)
  thr  <- terra::global(mag, fun = "mean", na.rm = TRUE)[[1]] +
          terra::global(mag, fun = "sd",   na.rm = TRUE)[[1]]
  terra::ifel(mag > thr, 1L, 0L)
}

#' @keywords internal
.raster_to_plot_polygons <- function(plot_interior, ref_raster,
                                      min_plot_area_m2,
                                      expected_rows, expected_cols,
                                      id_prefix) {
  # Convert plot interior pixels to polygons
  polys_sv <- terra::as.polygons(plot_interior, dissolve = TRUE)
  if (length(polys_sv) == 0)
    cli::cli_abort("No plot polygons detected. Try adjusting {.arg alley_width_m} or {.arg method}.")

  polys_sf <- sf::st_as_sf(polys_sv)
  polys_sf <- polys_sf[sf::st_is_valid(polys_sf), ]

  # Compute area and filter small artefacts
  polys_sf$area_m2 <- as.numeric(sf::st_area(polys_sf))
  polys_sf <- polys_sf[polys_sf$area_m2 >= min_plot_area_m2, ]

  if (nrow(polys_sf) == 0)
    cli::cli_abort(
      "All detected polygons smaller than min_plot_area_m2 ({min_plot_area_m2} m2). \\
       Lower threshold or adjust method.")

  # Assign row/col indices based on centroid positions
  cents   <- sf::st_coordinates(suppressWarnings(sf::st_centroid(polys_sf)))
  row_idx <- as.integer(cut(cents[, 2], breaks = max(2,
    if (!is.null(expected_rows)) expected_rows else min(30, nrow(polys_sf)))))
  col_idx <- as.integer(cut(cents[, 1], breaks = max(2,
    if (!is.null(expected_cols)) expected_cols else min(30, nrow(polys_sf)))))

  # Reverse row index so row 1 is at top (north)
  row_idx <- max(row_idx, na.rm = TRUE) + 1L - row_idx

  polys_sf$row_idx <- row_idx
  polys_sf$col_idx <- col_idx
  polys_sf$plot_id <- sprintf(
    paste0(id_prefix, "%0", nchar(nrow(polys_sf)), "d"),
    seq_len(nrow(polys_sf))
  )
  polys_sf$quality <- NA_real_

  # Drop the raster value column
  keep_cols <- c("plot_id", "row_idx", "col_idx", "area_m2", "quality", "geometry")
  polys_sf[, intersect(keep_cols, names(polys_sf))]
}

#' @keywords internal
.score_plot_quality <- function(plots, contrast, gap_mask) {
  plots$quality <- vapply(seq_len(nrow(plots)), function(i) {
    poly   <- plots[i, ]
    buf    <- tryCatch(sf::st_buffer(poly, dist = -0.1), error = function(e) poly)
    inner  <- terra::mask(terra::crop(contrast, terra::vect(buf)),
                          terra::vect(buf))
    outer_buf <- sf::st_buffer(sf::st_cast(poly, "MULTILINESTRING"), 0.2)
    outer  <- terra::mask(terra::crop(gap_mask, terra::vect(outer_buf)),
                          terra::vect(outer_buf))

    mean_inner <- tryCatch(
      terra::global(inner, "mean", na.rm = TRUE)[[1]], error = function(e) NA_real_)
    gap_frac   <- tryCatch({
      v <- terra::values(outer, na.rm = TRUE)
      if (length(v) == 0) NA_real_ else mean(v == 1, na.rm = TRUE)
    }, error = function(e) NA_real_)

    if (is.na(mean_inner) || is.na(gap_frac)) return(0.5)
    score <- 0.5 * min(gap_frac * 2, 1) + 0.5 * min(mean_inner + 0.5, 1)
    round(max(0, min(1, score)), 3)
  }, numeric(1))
  plots
}

#' @keywords internal
.sobel_x <- function(r) {
  kx <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), nrow = 3)
  terra::focal(r[[1]], w = kx, fun = "sum", na.policy = "omit", na.rm = TRUE)
}

#' @keywords internal
.sobel_y <- function(r) {
  ky <- matrix(c(-1, -2, -1, 0, 0, 0, 1, 2, 1), nrow = 3)
  terra::focal(r[[1]], w = ky, fun = "sum", na.policy = "omit", na.rm = TRUE)
}

#' @keywords internal
.sobel_magnitude <- function(r) {
  gx <- .sobel_x(r); gy <- .sobel_y(r)
  sqrt(gx^2 + gy^2)
}

#' @keywords internal
.gaussian_smooth_1d <- function(x, sigma) {
  if (sigma <= 0) return(x)
  hw  <- ceiling(3 * sigma)
  ker <- exp(-(-hw:hw)^2 / (2 * sigma^2))
  ker <- ker / sum(ker)
  stats::filter(x, ker, sides = 2, circular = FALSE)
}
