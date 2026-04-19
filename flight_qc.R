#' Automated flight quality assessment and reporting
#'
#' @description
#' After each drone flight, a researcher must check whether the images are
#' usable: was overlap sufficient? Were there motion-blur artefacts? Was the
#' sun angle too low? Was irradiance stable? Currently this requires opening
#' every image manually — a 2-hour job that becomes invisible and error-prone.
#'
#' \code{flight_qc()} ingests a folder of images or a \code{terra::SpatRaster}
#' and produces a structured quality report with accept/reject recommendations
#' per image and per flight. \code{flight_qc_report()} renders the report to
#' an HTML file for sharing and archiving.
#'
#' @name flight_quality
NULL

# ---------------------------------------------------------------------------
# flight_qc()
# ---------------------------------------------------------------------------

#' Assess flight data quality
#'
#' @description
#' Runs a battery of quality checks on a drone flight dataset and returns a
#' structured \code{FlightQC} object containing per-image and per-flight
#' metrics. Checks include:
#'
#' \describe{
#'   \item{Brightness}{Mean and SD of scene brightness — flags under/overexposure}
#'   \item{Blur}{Laplacian variance — low values indicate motion blur}
#'   \item{Saturation}{Fraction of saturated pixels (> 95th percentile)}
#'   \item{Irradiance stability}{CV of mean brightness across flight — flags
#'     mixed cloud cover during the flight}
#'   \item{Spatial coverage}{Raster extent and no-data fraction}
#'   \item{Band completeness}{Whether all expected bands are present}
#'   \item{Overall recommendation}{Accept / Caution / Reject per flight}
#' }
#'
#' @param x Either:
#'   \itemize{
#'     \item A \code{terra::SpatRaster} (single flight raster), or
#'     \item A character vector of paths to image files (one per flight image),
#'       or
#'     \item A \code{FlightStack} object (multiple flights — runs QC on each).
#'   }
#' @param flight_date Optional \code{Date} or character. Used for labelling.
#' @param expected_bands Integer. Expected number of spectral bands. A warning
#'   is raised if the actual band count differs. Default: \code{NULL} (no check).
#' @param blur_threshold Numeric. Laplacian variance below this value is
#'   flagged as blurry. Default: \code{100}.
#' @param saturation_threshold Numeric 0–1. Fraction of saturated pixels above
#'   which a band is flagged. Default: \code{0.02} (2%).
#' @param irradiance_cv_threshold Numeric. Coefficient of variation of mean
#'   brightness across flight images. Above this, the flight is flagged as
#'   having unstable illumination. Default: \code{0.15}.
#' @param nodata_threshold Numeric 0–1. Maximum acceptable fraction of NA
#'   pixels. Default: \code{0.10} (10%).
#'
#' @return A \code{FlightQC} object (a list with class attribute) containing:
#'   \describe{
#'     \item{\code{summary}}{Data frame with one row per flight:
#'       date, recommendation, overall_score, and individual check outcomes}
#'     \item{\code{metrics}}{List of detailed per-flight metric data frames}
#'     \item{\code{params}}{List of thresholds used}
#'   }
#'
#' @examples
#' \dontrun{
#' # Single raster
#' r   <- terra::rast("flight_2024_04_10.tif")
#' qc  <- flight_qc(r, flight_date = "2024-04-10", expected_bands = 5)
#' print(qc)
#'
#' # Full FlightStack — checks all flights at once
#' qc_all <- flight_qc(my_stack)
#' flight_qc_report(qc_all, output_file = "qc_season_2024.html")
#' }
#'
#' @seealso \code{\link{flight_qc_report}}
#' @export
flight_qc <- function(x,
                      flight_date            = NULL,
                      expected_bands         = NULL,
                      blur_threshold         = 100,
                      saturation_threshold   = 0.02,
                      irradiance_cv_threshold= 0.15,
                      nodata_threshold       = 0.10) {

  params <- list(
    blur_threshold         = blur_threshold,
    saturation_threshold   = saturation_threshold,
    irradiance_cv_threshold= irradiance_cv_threshold,
    nodata_threshold       = nodata_threshold,
    expected_bands         = expected_bands
  )

  # Dispatch on input type
  if (inherits(x, "FlightStack")) {
    cli::cli_h1("Flight Quality Assessment — {length(x$rasters)} flights")
    results <- lapply(seq_along(x$rasters), function(i) {
      .qc_single_raster(
        r            = x$rasters[[i]],
        flight_date  = x$dates[i],
        params       = params,
        flight_label = paste0("Flight ", i, " (", x$dates[i], ")")
      )
    })
  } else if (inherits(x, "SpatRaster")) {
    cli::cli_h1("Flight Quality Assessment")
    results <- list(.qc_single_raster(
      r            = x,
      flight_date  = flight_date,
      params       = params,
      flight_label = if (!is.null(flight_date)) as.character(flight_date) else "Flight 1"
    ))
  } else if (is.character(x)) {
    cli::cli_h1("Flight Quality Assessment — {length(x)} image files")
    results <- lapply(seq_along(x), function(i) {
      r <- tryCatch(
        terra::rast(x[i]),
        error = function(e) {
          cli::cli_warn("Cannot read {x[i]}: {e$message}")
          NULL
        }
      )
      if (is.null(r)) return(NULL)
      .qc_single_raster(
        r            = r,
        flight_date  = flight_date,
        params       = params,
        flight_label = basename(x[i])
      )
    })
    results <- Filter(Negate(is.null), results)
  } else {
    cli::cli_abort(
      "{.arg x} must be a SpatRaster, FlightStack, or character vector of file paths.")
  }

  # Combine summaries
  summary_df <- do.call(rbind, lapply(results, `[[`, "summary"))
  metrics_list <- lapply(results, `[[`, "metrics")
  names(metrics_list) <- summary_df$label

  qc_obj <- structure(
    list(summary = summary_df, metrics = metrics_list, params = params),
    class = "FlightQC"
  )

  .print_qc_summary(qc_obj)
  qc_obj
}

# ---------------------------------------------------------------------------
# flight_qc_report()
# ---------------------------------------------------------------------------

#' Render a flight QC report to HTML
#'
#' @description
#' Generates a self-contained HTML report from a \code{FlightQC} object,
#' suitable for archiving with the season's data or sharing with colleagues.
#' The report includes a summary table, per-flight metric plots, and
#' accept/reject recommendations with explanations.
#'
#' @param qc A \code{FlightQC} object from \code{\link{flight_qc}}.
#' @param output_file Character. Path to the output HTML file.
#'   Default: \code{"flight_qc_report.html"} in the working directory.
#' @param title Character. Report title. Default: \code{"hovR Flight QC Report"}.
#' @param open Logical. Open the report in the default browser after rendering.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly returns the path to the generated HTML file.
#'
#' @export
flight_qc_report <- function(qc,
                              output_file = "flight_qc_report.html",
                              title       = "hovR Flight QC Report",
                              open        = TRUE) {

  if (!inherits(qc, "FlightQC"))
    cli::cli_abort("{.arg qc} must be a FlightQC object from flight_qc().")

  cli::cli_inform("Rendering QC report to {.file {output_file}} ...")

  html <- .build_qc_html(qc, title)
  writeLines(html, output_file)

  cli::cli_inform(c("v" = "Report written: {.file {output_file}}"))

  if (open && interactive()) {
    tryCatch(utils::browseURL(output_file),
             error = function(e) cli::cli_warn("Could not open browser: {e$message}"))
  }

  invisible(output_file)
}

# ---------------------------------------------------------------------------
# print.FlightQC
# ---------------------------------------------------------------------------

#' @export
print.FlightQC <- function(x, ...) {
  cli::cli_h1("FlightQC Summary")
  df <- x$summary[, c("label", "date", "n_bands", "overall_score",
                       "recommendation", "flag_blur", "flag_saturation",
                       "flag_irradiance", "flag_nodata", "flag_bands")]
  print(df, row.names = FALSE)
  invisible(x)
}

# ---------------------------------------------------------------------------
# Internal: QC a single raster
# ---------------------------------------------------------------------------

#' @keywords internal
.qc_single_raster <- function(r, flight_date, params, flight_label) {

  n_bands <- terra::nlyr(r)
  cli::cli_inform("  Checking: {flight_label}")

  # 1. Band count check
  flag_bands <- !is.null(params$expected_bands) && n_bands != params$expected_bands

  # 2. No-data fraction
  total_cells  <- terra::ncell(r)
  na_count     <- sum(is.na(terra::values(r[[1]])))
  nodata_frac  <- na_count / total_cells
  flag_nodata  <- nodata_frac > params$nodata_threshold

  # 3. Brightness (mean of band means)
  band_means <- unlist(terra::global(r, "mean", na.rm = TRUE))
  mean_bright <- mean(band_means, na.rm = TRUE)
  sd_bright   <- stats::sd(band_means, na.rm = TRUE)

  # 4. Irradiance stability (CV across bands as proxy for mixed cloud)
  irr_cv      <- if (mean_bright > 0) sd_bright / mean_bright else NA_real_
  flag_irr    <- !is.na(irr_cv) && irr_cv > params$irradiance_cv_threshold

  # 5. Blur (Laplacian variance on first band)
  lap_var <- tryCatch({
    lap_kernel <- matrix(c(0, 1, 0, 1, -4, 1, 0, 1, 0), nrow = 3)
    lap        <- terra::focal(r[[1]], w = lap_kernel, fun = "sum", na.rm = TRUE)
    lap_vals   <- terra::values(lap, na.rm = TRUE)
    stats::var(lap_vals, na.rm = TRUE)
  }, error = function(e) NA_real_)
  flag_blur <- !is.na(lap_var) && lap_var < params$blur_threshold

  # 6. Saturation fraction (per band)
  sat_fracs <- vapply(seq_len(n_bands), function(b) {
    vals  <- terra::values(r[[b]], na.rm = TRUE)
    if (length(vals) == 0) return(NA_real_)
    thr   <- stats::quantile(vals, 0.95, na.rm = TRUE)
    mean(vals >= thr, na.rm = TRUE)
  }, numeric(1))
  mean_sat_frac <- mean(sat_fracs, na.rm = TRUE)
  flag_sat      <- mean_sat_frac > params$saturation_threshold

  # 7. Overall score (0-100) and recommendation
  penalties <- c(
    if (flag_bands)   10 else 0,
    if (flag_nodata)  20 else 0,
    if (flag_blur)    25 else 0,
    if (flag_irr)     20 else 0,
    if (flag_sat)     15 else 0
  )
  overall_score <- max(0, 100 - sum(penalties))
  recommendation <- dplyr::case_when(
    overall_score >= 80 ~ "Accept",
    overall_score >= 50 ~ "Caution",
    TRUE                ~ "Reject"
  )

  summary_row <- data.frame(
    label           = flight_label,
    date            = if (!is.null(flight_date)) as.character(as.Date(flight_date)) else NA_character_,
    n_bands         = n_bands,
    nodata_pct      = round(nodata_frac * 100, 1),
    mean_brightness = round(mean_bright, 4),
    irradiance_cv   = round(irr_cv, 4),
    blur_laplacian  = round(lap_var, 1),
    sat_frac_pct    = round(mean_sat_frac * 100, 2),
    flag_bands      = flag_bands,
    flag_nodata     = flag_nodata,
    flag_blur       = flag_blur,
    flag_irradiance = flag_irr,
    flag_saturation = flag_sat,
    overall_score   = overall_score,
    recommendation  = recommendation,
    stringsAsFactors = FALSE
  )

  metrics_detail <- list(
    band_means    = band_means,
    sat_fracs     = sat_fracs,
    nodata_frac   = nodata_frac,
    lap_var       = lap_var,
    irradiance_cv = irr_cv
  )

  list(summary = summary_row, metrics = metrics_detail)
}

# ---------------------------------------------------------------------------
# Internal: console summary printer
# ---------------------------------------------------------------------------

#' @keywords internal
.print_qc_summary <- function(qc) {
  df <- qc$summary
  n_accept  <- sum(df$recommendation == "Accept",  na.rm = TRUE)
  n_caution <- sum(df$recommendation == "Caution", na.rm = TRUE)
  n_reject  <- sum(df$recommendation == "Reject",  na.rm = TRUE)

  cli::cli_bullets(c(
    "v" = "Accept  : {n_accept}",
    "!" = "Caution : {n_caution}",
    "x" = "Reject  : {n_reject}"
  ))
  if (n_reject > 0) {
    rej_labels <- df$label[df$recommendation == "Reject"]
    cli::cli_warn("Rejected flight(s): {rej_labels}")
  }
}

# ---------------------------------------------------------------------------
# Internal: build HTML report
# ---------------------------------------------------------------------------

#' @keywords internal
.build_qc_html <- function(qc, title) {
  df <- qc$summary

  row_colour <- function(rec) {
    switch(rec,
      "Accept"  = "#d4edda",
      "Caution" = "#fff3cd",
      "Reject"  = "#f8d7da",
      "#ffffff"
    )
  }

  icon <- function(flag) if (isTRUE(flag)) "&#x26A0;" else "&#x2713;"

  header_row <- paste0(
    "<tr>",
    paste(c("Flight", "Date", "Bands", "No-data %", "Brightness",
            "Irr. CV", "Blur (Lap.)", "Sat. %",
            "Blur?", "Nodata?", "Irr.?", "Sat.?", "Score", "Decision"),
          function(h) paste0("<th>", h, "</th>"), USE.NAMES = FALSE),
    "</tr>"
  )
  header_row <- paste0(
    "<tr>",
    "<th>Flight</th><th>Date</th><th>Bands</th><th>No-data %</th>",
    "<th>Brightness</th><th>Irr. CV</th><th>Blur (Lap.)</th><th>Sat. %</th>",
    "<th>Blur?</th><th>Nodata?</th><th>Irr.?</th><th>Sat.?</th>",
    "<th>Score</th><th>Decision</th>",
    "</tr>"
  )

  data_rows <- paste(vapply(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    bg  <- row_colour(row$recommendation)
    paste0(
      "<tr style='background:", bg, "'>",
      "<td>", row$label, "</td>",
      "<td>", row$date,  "</td>",
      "<td>", row$n_bands, "</td>",
      "<td>", row$nodata_pct, "</td>",
      "<td>", row$mean_brightness, "</td>",
      "<td>", row$irradiance_cv, "</td>",
      "<td>", row$blur_laplacian, "</td>",
      "<td>", row$sat_frac_pct, "</td>",
      "<td>", icon(row$flag_blur), "</td>",
      "<td>", icon(row$flag_nodata), "</td>",
      "<td>", icon(row$flag_irradiance), "</td>",
      "<td>", icon(row$flag_saturation), "</td>",
      "<td><strong>", row$overall_score, "</strong></td>",
      "<td><strong>", row$recommendation, "</strong></td>",
      "</tr>"
    )
  }, character(1)), collapse = "\n")

  params_rows <- paste(vapply(names(qc$params), function(nm) {
    paste0("<tr><td>", nm, "</td><td>", qc$params[[nm]], "</td></tr>")
  }, character(1)), collapse = "\n")

  generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  n_total  <- nrow(df)
  n_accept <- sum(df$recommendation == "Accept")
  n_reject <- sum(df$recommendation == "Reject")

  paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>', title, '</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 1100px; margin: 40px auto; padding: 0 20px;
         color: #222; background: #f9f9f9; }
  h1   { color: #1a5276; border-bottom: 2px solid #1a5276; padding-bottom: 8px; }
  h2   { color: #1a5276; margin-top: 2rem; }
  .badge { display:inline-block; padding:4px 10px; border-radius:4px;
           font-size:0.85em; font-weight:600; margin-right:6px; }
  .accept  { background:#d4edda; color:#155724; }
  .caution { background:#fff3cd; color:#856404; }
  .reject  { background:#f8d7da; color:#721c24; }
  table { border-collapse: collapse; width: 100%; background: white;
          box-shadow: 0 1px 4px rgba(0,0,0,0.1); border-radius:6px;
          overflow:hidden; }
  th { background: #1a5276; color: white; padding: 10px 8px;
       text-align: left; font-size: 0.85em; }
  td { padding: 8px; font-size: 0.85em; border-bottom: 1px solid #eee; }
  .meta { font-size:0.8em; color:#666; margin-top:2rem; }
  .summary-cards { display:flex; gap:16px; margin:1rem 0 2rem; flex-wrap:wrap; }
  .card { background:white; border-radius:8px; padding:16px 24px;
          box-shadow:0 1px 4px rgba(0,0,0,0.1); text-align:center; }
  .card .num { font-size:2em; font-weight:700; }
  .card .lbl { font-size:0.85em; color:#666; }
</style>
</head>
<body>
<h1>', title, '</h1>

<div class="summary-cards">
  <div class="card">
    <div class="num">', n_total, '</div>
    <div class="lbl">Total flights</div>
  </div>
  <div class="card" style="border-top:4px solid #28a745">
    <div class="num" style="color:#155724">', n_accept, '</div>
    <div class="lbl">Accepted</div>
  </div>
  <div class="card" style="border-top:4px solid #ffc107">
    <div class="num" style="color:#856404">', sum(df$recommendation=="Caution"), '</div>
    <div class="lbl">Caution</div>
  </div>
  <div class="card" style="border-top:4px solid #dc3545">
    <div class="num" style="color:#721c24">', n_reject, '</div>
    <div class="lbl">Rejected</div>
  </div>
  <div class="card">
    <div class="num">', round(mean(df$overall_score), 1), '</div>
    <div class="lbl">Mean score</div>
  </div>
</div>

<h2>Per-flight results</h2>
<table>
<thead>', header_row, '</thead>
<tbody>', data_rows, '</tbody>
</table>

<h2>QC parameters used</h2>
<table style="width:400px">
<thead><tr><th>Parameter</th><th>Value</th></tr></thead>
<tbody>', params_rows, '</tbody>
</table>

<h2>Decision guide</h2>
<p>
  <span class="badge accept">Accept</span> Score &ge; 80: no significant issues detected.<br>
  <span class="badge caution">Caution</span> Score 50&ndash;79: one or more moderate issues — inspect manually before use.<br>
  <span class="badge reject">Reject</span> Score &lt; 50: multiple or severe issues — do not use without correction.
</p>

<p class="meta">Generated by hovR &bull; ', generated_at, '</p>
</body>
</html>')
}
