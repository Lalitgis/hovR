#' hovR: Drone Hyperspectral Time-Series Analysis
#'
#' @description
#' End-to-end tools for drone-based hyperspectral image and time-series analysis.
#' The package is organised into four modules:
#'
#' \describe{
#'   \item{Temporal integration}{
#'     \code{\link{flight_stack}}, \code{\link{compute_vi}},
#'     \code{\link{integrate_season}}, \code{\link{vi_derivative}},
#'     \code{\link{optimal_window}}, \code{\link{phenostage_tag}}
#'   }
#'   \item{Plot auto-segmentation}{
#'     \code{\link{segment_plots}}, \code{\link{refine_plot_boundaries}},
#'     \code{\link{plot_segmentation_qc}}
#'   }
#'   \item{Calibration panel detection}{
#'     \code{\link{detect_panels}}, \code{\link{calibrate_reflectance}},
#'     \code{\link{panel_qc_plot}}
#'   }
#'   \item{Flight quality assessment}{
#'     \code{\link{flight_qc}}, \code{\link{flight_qc_report}}
#'   }
#'   \item{Ground truth fusion}{
#'     \code{\link{extract_plot}}, \code{\link{fuse_ground}},
#'     \code{\link{plot_vi_traits}}
#'   }
#' }
#'
#' @section Typical workflow:
#' \preformatted{
#' library(hovR)
#'
#' # 1. Check flight quality before any analysis
#' qc <- flight_qc(terra::rast("flight.tif"), expected_bands = 5)
#' flight_qc_report(qc, "qc_report.html")
#'
#' # 2. Detect calibration panels and convert to reflectance
#' panels <- detect_panels(terra::rast("flight.tif"), n_panels = 2)
#' refl   <- calibrate_reflectance(terra::rast("flight.tif"), panels,
#'             known_reflectance = c(panel_1 = 0.50, panel_2 = 0.20))
#'
#' # 3. Auto-segment trial plots
#' plots <- segment_plots(refl, wavelengths = wl)
#' sf::st_write(plots, "trial_plots.gpkg")
#'
#' # 4. Build a time-series stack and compute VIs
#' fs <- flight_stack(raster_list, dates, wavelengths = wl)
#' vi <- compute_vi(fs, indices = c("NDVI", "NDRE", "LCI"))
#'
#' # 5. Integrate over the season
#' auc <- integrate_season(vi, index = "NDVI")
#'
#' # 6. Find which window best predicts yield
#' best <- optimal_window(vi, ground_truth = yield_df,
#'                         trait = "yield_t_ha", plots = plots)
#'
#' # 7. Fuse with field data and model
#' plot_data <- extract_plot(vi, plots)
#' fused     <- fuse_ground(plot_data, yield_df, traits = "yield_t_ha")
#' }
#'
#' @keywords internal
"_PACKAGE"
