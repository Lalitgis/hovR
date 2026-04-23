#' Generate a trial plot grid from user-specified geometry
#'
#' @description
#' When a user knows their trial layout (number of rows, columns, plot dimensions,
#' alley widths) this function generates an sf polygon layer matching the field
#' layout without relying on auto-segmentation. The user supplies an anchor point
#' (typically the NW corner of the first plot) either programmatically or by
#' clicking on a displayed raster.
#'
#' This is often more reliable than auto-segmentation for structured breeding
#' trials where the layout is precisely known. The resulting polygons can be
#' used directly with \code{\link{extract_plot}} and \code{\link{fuse_ground}}.
#'
#' @param anchor_x,anchor_y Numeric. Coordinates of the north-west corner of
#'   the first plot (row 1, column 1) in the raster's CRS. If NULL and
#'   \code{raster} is provided, an interactive click is used.
#' @param n_rows Integer. Number of plot rows in the trial.
#' @param n_cols Integer. Number of plot columns in the trial.
#' @param plot_length Numeric. Plot length along the row direction (metres).
#' @param plot_width Numeric. Plot width along the column direction (metres).
#' @param alley_row Numeric. Gap between plots in the same row (metres).
#'   Default 0 (plots abut horizontally).
#' @param alley_col Numeric. Gap between rows (metres). Default 0.
#' @param orientation_deg Numeric. Rotation angle in degrees (0 = plots aligned
#'   east-west). Use negative values for clockwise rotation. Default 0.
#' @param crs Character. CRS for the output polygons. If NULL and \code{raster}
#'   is supplied, uses the raster's CRS. Otherwise defaults to EPSG:4326.
#' @param id_prefix Character. Prefix for generated plot IDs. Default "P".
#' @param id_order Character. How to number plots: \code{"row_major"} (default)
#'   numbers left-to-right, row by row. \code{"col_major"} numbers top-to-bottom,
#'   column by column. \code{"snake"} alternates direction each row.
#' @param raster Optional \code{terra::SpatRaster}. If provided and
#'   \code{anchor_x}/\code{anchor_y} are NULL, the raster is displayed and the
#'   user clicks to set the anchor point.
#' @param starting_corner Character. Which corner of the first plot the anchor
#'   represents. One of \code{"NW"} (default), \code{"NE"}, \code{"SW"}, or
#'   \code{"SE"}.
#'
#' @return An \code{sf} polygon layer with columns \code{plot_id}, \code{row},
#'   \code{col}, \code{area_m2}, \code{quality} (set to 1.0 for compatibility
#'   with \code{\link{extract_plot}}).
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' # Programmatic anchor (known corner coordinates)
#' plots <- generate_plot_grid(
#'   anchor_x        = 434825.0,
#'   anchor_y        = 6668735.0,
#'   n_rows          = 14,
#'   n_cols          = 20,
#'   plot_length     = 2.5,
#'   plot_width      = 1.0,
#'   alley_row       = 0.3,
#'   alley_col       = 0.3,
#'   orientation_deg = 0,
#'   crs             = "EPSG:32722",
#'   id_prefix       = "W"
#' )
#'
#' # Interactive mode - click the NW corner on the map
#' plots <- generate_plot_grid(
#'   raster      = peak_raster_utm,
#'   n_rows      = 14,
#'   n_cols      = 20,
#'   plot_length = 2.5,
#'   plot_width  = 1.0,
#'   alley_row   = 0.3,
#'   alley_col   = 0.3
#' )
#' }
#'
#' @seealso \code{\link{segment_plots}}, \code{\link{extract_plot}}
#' @export
generate_plot_grid <- function(n_rows,
                                n_cols,
                                plot_length,
                                plot_width,
                                anchor_x        = NULL,
                                anchor_y        = NULL,
                                alley_row       = 0,
                                alley_col       = 0,
                                orientation_deg = 0,
                                crs             = NULL,
                                id_prefix       = "P",
                                id_order        = c("row_major", "col_major", "snake"),
                                raster          = NULL,
                                starting_corner = c("NW", "NE", "SW", "SE")) {

  id_order        <- match.arg(id_order)
  starting_corner <- match.arg(starting_corner)

  # ---- Validate inputs ----
  if (!is.numeric(n_rows) || n_rows <= 0 || n_rows != as.integer(n_rows))
    cli::cli_abort("{.arg n_rows} must be a positive integer.")
  if (!is.numeric(n_cols) || n_cols <= 0 || n_cols != as.integer(n_cols))
    cli::cli_abort("{.arg n_cols} must be a positive integer.")
  if (plot_length <= 0 || plot_width <= 0)
    cli::cli_abort("Plot dimensions must be positive.")
  if (alley_row < 0 || alley_col < 0)
    cli::cli_abort("Alley widths cannot be negative.")

  n_rows <- as.integer(n_rows)
  n_cols <- as.integer(n_cols)

  # ---- Resolve anchor ----
  if (is.null(anchor_x) || is.null(anchor_y)) {
    if (is.null(raster))
      cli::cli_abort("Provide either {.arg anchor_x}/{.arg anchor_y} or a {.arg raster} for clicking.")

    anchor <- .pick_anchor_interactive(raster, starting_corner)
    anchor_x <- anchor$x
    anchor_y <- anchor$y
  }

  # ---- Determine CRS ----
  if (is.null(crs)) {
    if (!is.null(raster)) {
      crs <- terra::crs(raster)
    } else {
      crs <- "EPSG:4326"
      cli::cli_warn("No CRS supplied - defaulting to EPSG:4326.")
    }
  }

  # ---- Compute plot step sizes ----
  step_x <- plot_length + alley_row    # spacing between plot origins along row
  step_y <- plot_width + alley_col     # spacing between rows

  cli::cli_h1("Generating plot grid")
  cli::cli_inform(c(
    "*" = "Layout      : {n_rows} rows x {n_cols} cols = {n_rows * n_cols} plots",
    "*" = "Plot size   : {plot_length} m x {plot_width} m",
    "*" = "Alleys      : row {alley_row} m, column {alley_col} m",
    "*" = "Orientation : {orientation_deg} deg",
    "*" = "Anchor ({starting_corner}) : ({round(anchor_x, 2)}, {round(anchor_y, 2)})"
  ))

  # ---- Build grid (relative to anchor, pre-rotation) ----
  # We build with anchor treated as NW corner, then transform for other corners
  polys_list <- vector("list", n_rows * n_cols)
  idx <- 1L

  for (r in seq_len(n_rows)) {
    for (c in seq_len(n_cols)) {
      # Top-left corner of this plot (relative to anchor)
      x0_rel <- (c - 1) * step_x
      y0_rel <- -(r - 1) * step_y    # negative because y decreases as we go south

      # Four corners (closing polygon)
      corners <- rbind(
        c(x0_rel,              y0_rel),
        c(x0_rel + plot_length, y0_rel),
        c(x0_rel + plot_length, y0_rel - plot_width),
        c(x0_rel,              y0_rel - plot_width),
        c(x0_rel,              y0_rel)
      )

      # Apply rotation if requested
      if (orientation_deg != 0) {
        theta <- orientation_deg * pi / 180
        rot_matrix <- matrix(c(cos(theta), -sin(theta),
                                sin(theta),  cos(theta)),
                              nrow = 2, byrow = TRUE)
        corners <- corners %*% rot_matrix
      }

      # Apply starting corner offset if not NW
      corners <- .apply_corner_offset(
        corners, starting_corner, plot_length, plot_width,
        n_rows, n_cols, step_x, step_y
      )

      # Translate by anchor
      corners[, 1] <- corners[, 1] + anchor_x
      corners[, 2] <- corners[, 2] + anchor_y

      polys_list[[idx]] <- sf::st_polygon(list(corners))
      idx <- idx + 1L
    }
  }

  # ---- Build sf object ----
  plots_sf <- sf::st_sf(
    plot_id = .make_plot_ids(n_rows, n_cols, id_prefix, id_order),
    row     = rep(seq_len(n_rows), each = n_cols),
    col     = rep(seq_len(n_cols), times = n_rows),
    geometry = sf::st_sfc(polys_list, crs = crs)
  )

  plots_sf$area_m2 <- as.numeric(sf::st_area(plots_sf))
  plots_sf$quality <- 1.0

  plots_sf <- plots_sf[, c("plot_id", "row", "col", "area_m2", "quality")]

  cli::cli_inform(c(
    "v" = "Generated {nrow(plots_sf)} plots",
    "*" = "Mean area: {round(mean(plots_sf$area_m2), 3)} m^2"
  ))

  plots_sf
}


#' Preview plot grid overlaid on a raster before finalising
#'
#' @description
#' Helper to visually check whether the generated plot grid aligns correctly
#' with field plots in a raster. Use this before running \code{extract_plot()}.
#'
#' @param plots sf polygon layer (output of \code{\link{generate_plot_grid}})
#' @param raster terra::SpatRaster - the image to overlay on
#' @param band Integer. Which band of \code{raster} to display. Default 1.
#' @param show_ids Logical. Whether to label plots with their IDs.
#'   Default FALSE for grids larger than 50 plots.
#' @param border_colour Character. Polygon border colour. Default "cyan".
#'
#' @return Invisibly returns the plots layer.
#' @export
preview_plot_grid <- function(plots, raster,
                               band = 1,
                               show_ids = NULL,
                               border_colour = "cyan") {

  if (is.null(show_ids)) show_ids <- nrow(plots) <= 50

  terra::plot(
    raster[[band]],
    main = paste("Plot grid preview -", nrow(plots), "plots"),
    col  = hcl.colors(50, "RdYlGn")
  )

  plot(sf::st_geometry(plots),
       add    = TRUE,
       border = border_colour,
       lwd    = 0.8)

  if (show_ids) {
    centroids <- suppressWarnings(sf::st_centroid(plots))
    coords    <- sf::st_coordinates(centroids)
    graphics::text(coords[, 1], coords[, 2],
                    labels = plots$plot_id,
                    cex = 0.5, col = "black")
  }

  invisible(plots)
}


# ---- Internal helpers ----

#' @keywords internal
.pick_anchor_interactive <- function(raster, corner_label) {
  cli::cli_inform(c(
    "i" = "Interactive mode: click the {corner_label} corner of plot (row 1, col 1).",
    "*" = "Tip: zoom in first for accuracy."
  ))

  terra::plot(raster[[1]],
              main = paste("Click the", corner_label,
                           "corner of plot 1 (row 1, col 1)"),
              col  = hcl.colors(50, "RdYlGn"))

  loc <- graphics::locator(n = 1)
  if (is.null(loc))
    cli::cli_abort("No click detected. Provide anchor_x and anchor_y instead.")

  cli::cli_inform("Anchor set: ({round(loc$x, 2)}, {round(loc$y, 2)})")
  list(x = loc$x, y = loc$y)
}

#' @keywords internal
.apply_corner_offset <- function(corners, starting_corner,
                                  plot_length, plot_width,
                                  n_rows, n_cols, step_x, step_y) {
  # Default build assumes anchor = NW corner of plot (1,1)
  # If user specifies a different corner, shift the entire grid
  switch(starting_corner,
    NW = corners,
    NE = { corners[, 1] <- corners[, 1] - (n_cols - 1) * step_x - plot_length; corners },
    SW = { corners[, 2] <- corners[, 2] + (n_rows - 1) * step_y + plot_width; corners },
    SE = {
      corners[, 1] <- corners[, 1] - (n_cols - 1) * step_x - plot_length
      corners[, 2] <- corners[, 2] + (n_rows - 1) * step_y + plot_width
      corners
    }
  )
}

#' @keywords internal
.make_plot_ids <- function(n_rows, n_cols, prefix, order) {
  total <- n_rows * n_cols
  ids   <- character(total)
  idx   <- 1L

  if (order == "row_major") {
    for (r in seq_len(n_rows)) {
      for (c in seq_len(n_cols)) {
        ids[idx] <- sprintf("%s%03d", prefix, idx)
        idx <- idx + 1L
      }
    }
  } else if (order == "col_major") {
    rc_map <- matrix(0L, nrow = n_rows, ncol = n_cols)
    n <- 1L
    for (c in seq_len(n_cols)) {
      for (r in seq_len(n_rows)) {
        rc_map[r, c] <- n
        n <- n + 1L
      }
    }
    for (r in seq_len(n_rows)) {
      for (c in seq_len(n_cols)) {
        ids[idx] <- sprintf("%s%03d", prefix, rc_map[r, c])
        idx <- idx + 1L
      }
    }
  } else if (order == "snake") {
    for (r in seq_len(n_rows)) {
      cols_in_order <- if (r %% 2 == 1) seq_len(n_cols) else rev(seq_len(n_cols))
      for (c in cols_in_order) {
        ids[idx] <- sprintf("%s%03d", prefix, idx)
        idx <- idx + 1L
      }
    }
  }
  ids
}
