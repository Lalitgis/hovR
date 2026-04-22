#' Ground truth fusion for drone hyperspectral time-series
#'
#' @description
#' Spatially matches field-measured trait data (from destructive sampling,
#' hand-held meters, or lab analysis) to hyperspectral raster pixels using
#' plot polygon boundaries. Produces tidy paired data frames ready for
#' modelling with \code{lm()}, \code{pls::plsr()}, or \code{tidymodels}.

# ---- extract_plot -----------------------------------------------------------

#' Extract per-plot VI statistics from a VIStack
#'
#' @description
#' Computes summary statistics (mean, sd, median, CV) of VI values within each
#' plot polygon across all flights. Returns a tidy long-format data frame with
#' one row per plot x flight x index combination.
#'
#' @param vi_stack A \code{VIStack} object.
#' @param plots An \code{sf} polygon layer or \code{terra::SpatVector} with at
#'   least a \code{plot_id} column. All other columns are retained in the output.
#' @param indices Character vector. Which indices to extract. Default: all.
#' @param stats Character vector of statistics to compute per plot.
#'   Options: \code{"mean"} (default), \code{"sd"}, \code{"median"}, \code{"cv"}.
#' @param min_pixels Integer. Minimum number of valid (non-NA) pixels required
#'   within a plot polygon. Plots with fewer pixels are flagged as \code{NA}
#'   with a warning. Default: 5.
#'
#' @return A tidy data frame with columns:
#'   \code{plot_id}, \code{date}, \code{bbch} (if tagged), \code{index},
#'   and one column per requested statistic.
#'
#' @export
extract_plot <- function(vi_stack,
                         plots,
                         indices    = NULL,
                         stats      = "mean",
                         min_pixels = 5L) {

  if (is.null(indices)) indices <- names(vi_stack)
  stats <- match.arg(stats, c("mean", "sd", "median", "cv"), several.ok = TRUE)

  plots_sv <- if (inherits(plots, "sf")) terra::vect(plots) else plots
  if (!"plot_id" %in% names(plots_sv))
    cli::cli_abort("{.arg plots} must contain a 'plot_id' column.")

  all_dates <- as.Date(attr(vi_stack, "dates") %||% names(vi_stack[[indices[1]]]))
  rows <- list()

  cli::cli_inform(
    "Extracting plot values for {length(indices)} index/indices across {length(all_dates)} flights..."
  )

  for (idx in indices) {
    for (dt in as.character(all_dates)) {
      r    <- vi_stack[[idx]][[dt]]
      extr <- terra::extract(r, plots_sv, fun = NULL, na.rm = TRUE, ID = TRUE)

      # Aggregate per plot
      for (pid in unique(plots_sv$plot_id)) {
        plot_rows <- extr[extr$ID == which(plots_sv$plot_id == pid), 2]
        valid     <- plot_rows[!is.na(plot_rows)]

        row <- list(
          plot_id = pid,
          date    = as.Date(dt),
          index   = idx
        )

        if (length(valid) < min_pixels) {
          cli::cli_warn(
            "Plot {pid}, {dt}: only {length(valid)} valid pixels (<{min_pixels}). Setting NA.")
          for (s in stats) row[[s]] <- NA_real_
        } else {
          if ("mean"   %in% stats) row$mean   <- mean(valid)
          if ("sd"     %in% stats) row$sd     <- stats::sd(valid)
          if ("median" %in% stats) row$median <- stats::median(valid)
          if ("cv"     %in% stats) row$cv     <- stats::sd(valid) / mean(valid) * 100
        }
        rows[[length(rows) + 1]] <- as.data.frame(row)
      }
    }
  }

  result <- do.call(rbind, rows)
  result$date <- as.Date(result$date)
  result
}

# ---- fuse_ground ------------------------------------------------------------

#' Fuse extracted VI data with field-measured trait data
#'
#' @description
#' Joins the output of \code{extract_plot()} with ground-truth measurements.
#' Supports exact date matching and lagged matching (e.g. when field data
#' was collected several days before or after the drone flight).
#'
#' @param vi_data A data frame from \code{\link{extract_plot}}.
#' @param ground_truth A data frame with at minimum:
#'   \describe{
#'     \item{\code{plot_id}}{Matches \code{plot_id} in \code{vi_data}}
#'     \item{\code{date}}{Date of field measurement (as \code{Date} or string)}
#'     \item{\code{...}}{One or more trait columns (yield, LAI, nitrogen, etc.)}
#'   }
#' @param traits Character vector. Names of trait columns to retain.
#'   If \code{NULL}, all non-join columns in \code{ground_truth} are kept.
#' @param max_lag_days Integer. Maximum allowable days between flight date and
#'   field measurement date. Pairs with larger gaps are dropped. Default: 7.
#' @param allow_missing Logical. If \code{TRUE}, plots with no matching
#'   ground-truth measurement are retained with \code{NA} trait values.
#'   Default: \code{FALSE}.
#'
#' @return A tidy data frame combining VI statistics and trait values, with an
#'   added \code{date_lag_days} column showing the flight-to-measurement gap.
#'
#' @export
fuse_ground <- function(vi_data,
                        ground_truth,
                        traits       = NULL,
                        max_lag_days = 7L,
                        allow_missing = FALSE) {

  ground_truth$date <- as.Date(ground_truth$date)
  vi_data$date      <- as.Date(vi_data$date)

  if (is.null(traits)) {
    traits <- setdiff(names(ground_truth), c("plot_id", "date"))
  }

  cli::cli_inform(
    "Fusing {nrow(vi_data)} VI observations with {nrow(ground_truth)} \\
     ground truth records. Max lag: {max_lag_days} days.")

  # For each VI row, find the closest ground-truth date per plot
  fused <- lapply(seq_len(nrow(vi_data)), function(i) {
    row     <- vi_data[i, ]
    gt_plot <- ground_truth[ground_truth$plot_id == row$plot_id, ]
    if (nrow(gt_plot) == 0) {
      if (!allow_missing) return(NULL)
      result <- row
      result[traits] <- NA_real_
      result$date_lag_days <- NA_integer_
      return(result)
    }

    lags    <- abs(as.integer(gt_plot$date - row$date))
    best    <- which.min(lags)

    if (lags[best] > max_lag_days) return(NULL)

    result               <- row
    result[traits]       <- gt_plot[best, traits]
    result$date_lag_days <- lags[best]
    result
  })

  out <- do.call(rbind, Filter(Negate(is.null), fused))

  n_matched  <- sum(!is.na(out$date_lag_days))
  n_dropped  <- nrow(vi_data) - nrow(out)
  cli::cli_inform(
    "{n_matched} matches retained. {n_dropped} rows dropped \\
     (lag > {max_lag_days} days or no ground truth).")
  out
}

# ---- plot_vi_traits ---------------------------------------------------------

#' Scatter plot: per-plot VI versus field-measured trait
#'
#' @description
#' Produces a ggplot2 scatter with a LOESS or linear fit, coloured by flight
#' date or BBCH stage. Designed for rapid visual model diagnostics.
#'
#' @param fused_data A data frame from \code{\link{fuse_ground}}.
#' @param vi_col Character. Name of the VI column (e.g. \code{"mean"}).
#' @param trait_col Character. Name of the trait column.
#' @param colour_by One of \code{"date"} or \code{"bbch"}.
#' @param facet_by Optional variable name to facet on (e.g. \code{"index"}).
#' @param fit One of \code{"lm"}, \code{"loess"}, or \code{"none"}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_vi_traits <- function(fused_data,
                           vi_col     = "mean",
                           trait_col,
                           colour_by  = c("date", "bbch"),
                           facet_by   = NULL,
                           fit        = c("lm", "loess", "none")) {

  colour_by <- match.arg(colour_by)
  fit       <- match.arg(fit)

  colour_var <- if (colour_by == "date") "date" else "bbch"

  p <- ggplot2::ggplot(fused_data,
    ggplot2::aes(
      x      = .data[[vi_col]],
      y      = .data[[trait_col]],
      colour = .data[[colour_var]]
    )
  ) +
  ggplot2::geom_point(alpha = 0.7, size = 2) +
  ggplot2::labs(
    x      = vi_col,
    y      = trait_col,
    colour = colour_by
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "right"
  )

  if (fit == "lm") {
    p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "grey30",
                                   linewidth = 0.8)
  } else if (fit == "loess") {
    p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, colour = "grey30",
                                   linewidth = 0.8)
  }

  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
  }

  p
}
