#' Temporal integration of drone hyperspectral VI time-series
#'
#' @description
#' These are the core functions that no existing R package provides.
#' They treat the seasonal sequence of drone flights as a coherent temporal
#' object, enabling:
#'
#' \itemize{
#'   \item Seasonal area-under-the-curve (AUC) integration per pixel/plot
#'   \item Rate-of-change (first derivative) of any VI over time
#'   \item Automatic detection of optimal flight windows for trait prediction
#'   \item BBCH phenological stage tagging and subsetting
#' }
#'
#' Research context: integrating across the full time series correlates more
#' strongly with yield, heading date, and maturity than any single flight date
#' (Killian et al., 2025; Kuenzer et al., 2015).

# ---- integrate_season -------------------------------------------------------

#' Compute pixel-wise seasonal integral (AUC) of a vegetation index
#'
#' @description
#' Computes the area under the VI curve over the growing season using the
#' trapezoidal rule. Handles irregular flight spacing — flights do not need
#' to be evenly spaced. Optionally restricts integration to a phenological
#' window (e.g. only during grain filling).
#'
#' @param vi_stack A \code{VIStack} object (output of \code{\link{compute_vi}}),
#'   or a named list of single-band \code{terra::SpatRaster} objects keyed by
#'   ISO date string.
#' @param index Character. Which index to integrate. Must be present in
#'   \code{vi_stack}. Default: first available index.
#' @param date_from,date_to Optional \code{Date} or character limits to restrict
#'   the integration window. If \code{NULL}, the full season is used.
#' @param bbch_from,bbch_to Optional BBCH stage codes (character). If provided,
#'   the integration window is bounded by the first/last flight tagged with
#'   these BBCH stages. Requires \code{vi_stack} to retain BBCH metadata, i.e.
#'   you must pass the original \code{FlightStack} via the \code{stack} arg.
#' @param stack A \code{FlightStack} object, required only when \code{bbch_from}
#'   or \code{bbch_to} are specified.
#' @param method Integration method. \code{"trapezoid"} (default) or
#'   \code{"step"} (block/left-endpoint rule, simpler but less accurate).
#' @param units Character: \code{"days"} (default) or \code{"doy"}
#'   (day of year). Controls the time axis unit for the integral.
#'
#' @return A single-band \code{terra::SpatRaster} with the seasonal VI integral
#'   (AUC) at each pixel. Units are VI-units × days (or × doy).
#'
#' @examples
#' \dontrun{
#' vi   <- compute_vi(my_stack, indices = c("NDVI", "NDRE"))
#' auc  <- integrate_season(vi, index = "NDVI")
#' auc_gf <- integrate_season(vi, index = "NDVI",
#'                             bbch_from = "71", bbch_to = "89",
#'                             stack = my_stack)
#' terra::plot(auc, main = "Seasonal NDVI integral")
#' }
#'
#' @export
integrate_season <- function(vi_stack,
                             index      = NULL,
                             date_from  = NULL,
                             date_to    = NULL,
                             bbch_from  = NULL,
                             bbch_to    = NULL,
                             stack      = NULL,
                             method     = c("trapezoid", "step"),
                             units      = c("days", "doy")) {

  method <- match.arg(method)
  units  <- match.arg(units)
  index  <- .resolve_index(vi_stack, index)

  rasters_by_date <- vi_stack[[index]]
  all_dates       <- as.Date(names(rasters_by_date))

  # Resolve BBCH-based window limits
  if (!is.null(bbch_from) || !is.null(bbch_to)) {
    if (is.null(stack) || !inherits(stack, "FlightStack"))
      cli::cli_abort("Provide {.arg stack} (a FlightStack) when using bbch_from/bbch_to.")
    if (is.null(stack$bbch))
      cli::cli_abort("FlightStack has no BBCH stage labels. Run {.fn phenostage_tag} first.")
    if (!is.null(bbch_from))
      date_from <- min(stack$dates[stack$bbch >= bbch_from], na.rm = TRUE)
    if (!is.null(bbch_to))
      date_to   <- max(stack$dates[stack$bbch <= bbch_to],   na.rm = TRUE)
    cli::cli_inform("BBCH window resolved: {date_from} to {date_to}")
  }

  # Filter dates
  keep <- rep(TRUE, length(all_dates))
  if (!is.null(date_from)) keep <- keep & (all_dates >= as.Date(date_from))
  if (!is.null(date_to))   keep <- keep & (all_dates <= as.Date(date_to))
  if (sum(keep) < 2)
    cli::cli_abort("At least 2 flights required in window. Only {sum(keep)} found.")

  sel_dates   <- all_dates[keep]
  sel_rasters <- rasters_by_date[keep]

  cli::cli_inform(
    "Integrating {index} over {length(sel_dates)} flights \\
     ({min(sel_dates)} to {max(sel_dates)}) using {method} rule.")

  # Compute time axis in chosen units
  t_axis <- if (units == "doy") {
    as.integer(format(sel_dates, "%j"))
  } else {
    as.integer(sel_dates - sel_dates[1])
  }

  # Pixel-wise trapezoidal or step integration
  ref <- sel_rasters[[1]]
  auc <- terra::rast(ref)
  terra::values(auc) <- 0

  n <- length(sel_rasters)

  if (method == "trapezoid") {
    for (i in seq_len(n - 1)) {
      dt   <- t_axis[i + 1] - t_axis[i]
      trap <- (sel_rasters[[i]] + sel_rasters[[i + 1]]) / 2 * dt
      auc  <- auc + trap
    }
  } else {  # step
    for (i in seq_len(n - 1)) {
      dt   <- t_axis[i + 1] - t_axis[i]
      auc  <- auc + sel_rasters[[i]] * dt
    }
  }

  names(auc) <- paste0(index, "_AUC_", units)
  terra::metags(auc) <- list(
    index     = index,
    method    = method,
    date_from = as.character(min(sel_dates)),
    date_to   = as.character(max(sel_dates)),
    n_flights = as.character(length(sel_dates))
  )
  auc
}

# ---- vi_derivative ----------------------------------------------------------

#' Compute the temporal derivative (rate of change) of a VI time-series
#'
#' @description
#' Computes the first temporal derivative of a vegetation index at each pixel:
#' ΔVI / Δt between consecutive flight dates. Returns a list of rasters, one
#' per inter-flight interval. Useful for detecting the onset of stress,
#' rapid green-up, or senescence.
#'
#' @param vi_stack A \code{VIStack} object.
#' @param index Character. Index to differentiate.
#' @param smooth Logical. If \code{TRUE}, applies a 3-point moving average
#'   across dates before differencing to reduce noise. Default: \code{FALSE}.
#' @param units Character. Time units for the denominator: \code{"days"} (default)
#'   or \code{"weeks"}.
#'
#' @return A named list of \code{terra::SpatRaster} objects. Names are interval
#'   labels of the form \code{"2024-04-01_to_2024-05-15"}.
#'
#' @export
vi_derivative <- function(vi_stack,
                          index  = NULL,
                          smooth = FALSE,
                          units  = c("days", "weeks")) {

  units <- match.arg(units)
  index <- .resolve_index(vi_stack, index)
  rlist <- vi_stack[[index]]
  dates <- as.Date(names(rlist))
  n     <- length(rlist)

  if (n < 2)
    cli::cli_abort("At least 2 flights required to compute derivative.")

  if (smooth) {
    cli::cli_inform("Applying 3-point temporal smoothing before differentiation.")
    rlist <- .smooth_vi_stack(rlist)
  }

  divisor <- if (units == "weeks") 7 else 1

  intervals <- lapply(seq_len(n - 1), function(i) {
    dt    <- as.integer(dates[i + 1] - dates[i]) / divisor
    deriv <- (rlist[[i + 1]] - rlist[[i]]) / dt
    names(deriv) <- paste0(index, "_ddt")
    deriv
  })

  interval_names <- paste0(
    as.character(dates[-n]), "_to_", as.character(dates[-1])
  )
  stats::setNames(intervals, interval_names)
}

# ---- optimal_window ---------------------------------------------------------

#' Find flight dates or windows with highest predictive power for a trait
#'
#' @description
#' Correlates per-plot VI values (across all flight dates) with a ground-truth
#' trait variable (e.g. yield, heading date, nitrogen content). Returns
#' a ranked data frame showing which individual flight dates and which
#' cumulative windows explain the most trait variance.
#'
#' @details
#' This implements the core recommendation from Killian et al. (2025): rather
#' than relying on a single "best" flight date found by trial and error,
#' systematically evaluate all dates and all contiguous window integrals.
#'
#' @param vi_stack A \code{VIStack} object.
#' @param ground_truth A data frame with at least two columns:
#'   \describe{
#'     \item{\code{plot_id}}{Plot identifier matching \code{plots}}
#'     \item{\code{<trait>}}{The trait value to predict}
#'   }
#' @param trait Character. Name of the trait column in \code{ground_truth}.
#' @param plots An \code{sf} or \code{SpatVector} polygon layer with a
#'   \code{plot_id} column. VI values are spatially averaged within each polygon.
#' @param index Character. VI to use. Default: first available.
#' @param method Correlation method: \code{"pearson"} (default), \code{"spearman"}.
#' @param n_top Integer. Number of top-ranked dates/windows to return. Default 10.
#'
#' @return A data frame ranked by \code{abs(r)}, with columns:
#'   \code{type} (single/window), \code{date_from}, \code{date_to},
#'   \code{r}, \code{r2}, \code{p_value}, \code{n_plots}.
#'
#' @export
optimal_window <- function(vi_stack,
                           ground_truth,
                           trait,
                           plots,
                           index  = NULL,
                           method = c("pearson", "spearman"),
                           n_top  = 10) {

  method <- match.arg(method)
  index  <- .resolve_index(vi_stack, index)
  rlist  <- vi_stack[[index]]
  dates  <- as.Date(names(rlist))
  n      <- length(dates)

  if (!trait %in% names(ground_truth))
    cli::cli_abort("Column '{trait}' not found in ground_truth.")

  cli::cli_inform("Extracting plot-mean {index} for {n} flights across {nrow(plots)} plots...")

  # Extract per-plot mean VI for each flight
  plot_vi <- lapply(stats::setNames(nm = as.character(dates)), function(dt) {
    r    <- rlist[[dt]]
    extr <- terra::extract(r, terra::vect(plots), fun = mean, na.rm = TRUE)
    data.frame(plot_id = plots$plot_id, vi_mean = extr[, 2])
  })

  # Single-date correlations
  single_rows <- lapply(as.character(dates), function(dt) {
    df  <- merge(plot_vi[[dt]], ground_truth, by = "plot_id")
    if (nrow(df) < 3) return(NULL)
    ct  <- stats::cor.test(df$vi_mean, df[[trait]], method = method)
    data.frame(type = "single", date_from = as.Date(dt), date_to = as.Date(dt),
               r = ct$estimate, r2 = ct$estimate^2, p_value = ct$p.value,
               n_plots = nrow(df))
  })

  # Window (AUC) correlations — all contiguous pairs
  window_rows <- list()
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      sub_stack <- list(
        rasters     = rlist[i:j],
        dates       = dates[i:j],
        wavelengths = attr(vi_stack, "wavelengths") %||% numeric(0),
        bbch        = NULL,
        meta        = NULL
      )
      # Simple inline AUC for this window
      t_ax  <- as.integer(dates[i:j] - dates[i])
      auc_r <- .inline_trapz(rlist[i:j], t_ax)

      extr <- terra::extract(auc_r, terra::vect(plots), fun = mean, na.rm = TRUE)
      df   <- merge(
        data.frame(plot_id = plots$plot_id, vi_mean = extr[, 2]),
        ground_truth, by = "plot_id"
      )
      if (nrow(df) < 3) next
      ct <- stats::cor.test(df$vi_mean, df[[trait]], method = method)
      window_rows[[length(window_rows) + 1]] <- data.frame(
        type      = "window",
        date_from = dates[i],
        date_to   = dates[j],
        r         = ct$estimate,
        r2        = ct$estimate^2,
        p_value   = ct$p.value,
        n_plots   = nrow(df)
      )
    }
  }

  all_rows <- do.call(rbind, c(Filter(Negate(is.null), single_rows), window_rows))
  all_rows <- all_rows[order(-abs(all_rows$r)), ]
  rownames(all_rows) <- NULL
  head(all_rows, n_top)
}

# ---- phenostage_tag ---------------------------------------------------------

#' Tag FlightStack flights with BBCH phenological stage codes
#'
#' @description
#' Adds BBCH (Biologische Bundesanstalt, Bundessortenamt and CHemische Industrie)
#' growth stage codes to a \code{FlightStack}. Codes can be supplied manually
#' or interpolated from a sparse reference table (e.g. from a field notebook
#' with only a few dated observations).
#'
#' @param stack A \code{FlightStack} object.
#' @param bbch_obs A data frame with columns \code{date} (Date) and
#'   \code{bbch} (character BBCH code). Observations need not coincide with
#'   flight dates — stages are interpolated to all flight dates.
#' @param method Interpolation method: \code{"nearest"} (default — assigns
#'   the nearest observed BBCH code to each flight date) or \code{"linear"}
#'   (linearly interpolates the numeric BBCH code and rounds).
#'
#' @return A \code{FlightStack} with the \code{bbch} slot populated.
#'
#' @examples
#' \dontrun{
#' obs <- data.frame(
#'   date = as.Date(c("2024-03-20", "2024-04-25", "2024-06-01")),
#'   bbch = c("21", "55", "89")
#' )
#' my_stack <- phenostage_tag(my_stack, bbch_obs = obs)
#' summary(my_stack)
#' }
#'
#' @export
phenostage_tag <- function(stack,
                           bbch_obs,
                           method = c("nearest", "linear")) {

  method <- match.arg(method)
  if (!inherits(stack, "FlightStack"))
    cli::cli_abort("{.arg stack} must be a FlightStack.")

  obs_dates <- as.Date(bbch_obs$date)
  obs_bbch  <- as.integer(bbch_obs$bbch)
  fl_dates  <- stack$dates

  tagged <- vapply(fl_dates, function(fd) {
    if (method == "nearest") {
      nearest_idx <- which.min(abs(obs_dates - fd))
      bbch_obs$bbch[nearest_idx]
    } else {
      interp_val <- stats::approx(
        x      = as.integer(obs_dates),
        y      = obs_bbch,
        xout   = as.integer(fd),
        method = "linear",
        rule   = 2
      )$y
      as.character(round(interp_val))
    }
  }, character(1))

  stack$bbch <- tagged
  cli::cli_inform("BBCH stages tagged for {length(fl_dates)} flights:")
  cli::cli_dl(stats::setNames(tagged, as.character(fl_dates)))
  stack
}

# ---- Internal helpers -------------------------------------------------------

#' @keywords internal
.resolve_index <- function(vi_stack, index) {
  available <- names(vi_stack)
  if (is.null(index)) return(available[1])
  if (!index %in% available)
    cli::cli_abort("Index '{index}' not in VIStack. Available: {available}.")
  index
}

#' @keywords internal
.smooth_vi_stack <- function(rlist) {
  n <- length(rlist)
  lapply(seq_len(n), function(i) {
    idxs <- max(1, i-1):min(n, i+1)
    Reduce(`+`, rlist[idxs]) / length(idxs)
  })
}

#' @keywords internal
.inline_trapz <- function(rlist, t_axis) {
  auc <- terra::rast(rlist[[1]]); terra::values(auc) <- 0
  for (i in seq_len(length(rlist) - 1)) {
    dt  <- t_axis[i + 1] - t_axis[i]
    auc <- auc + (rlist[[i]] + rlist[[i + 1]]) / 2 * dt
  }
  auc
}

#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b
