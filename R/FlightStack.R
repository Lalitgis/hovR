#' FlightStack: A multi-date drone hyperspectral raster collection
#'
#' @description
#' The \code{FlightStack} class is the central data structure of \pkg{hovR}.
#' It holds an ordered collection of single-date hyperspectral rasters
#' (as \code{terra::SpatRaster} objects), their acquisition dates, wavelength
#' vectors, phenological stage labels, and optional metadata per flight.
#'
#' No existing R package (hsdar, hyperSpec, stars) provides this abstraction.
#' A FlightStack is the prerequisite for every temporal operation in the package:
#' \code{integrate_season()}, \code{vi_derivative()}, \code{optimal_window()},
#' and \code{phenostage_tag()}.
#'
#' @slot rasters  Named list of \code{terra::SpatRaster} objects, one per flight
#' @slot dates    \code{Date} vector, length == length(rasters)
#' @slot wavelengths Numeric vector of band centre wavelengths (nm)
#' @slot bbch     Optional character vector of BBCH stage labels per flight
#' @slot meta     Optional data frame with per-flight metadata columns
#'
#' @name FlightStack
#' @aliases FlightStack-class
NULL

# ---- Constructor ------------------------------------------------------------

#' Create a FlightStack from a list of rasters and dates
#'
#' @param rasters A named or unnamed list of \code{terra::SpatRaster} objects.
#'   Each raster must have the same number of bands and cover the same spatial
#'   extent (or be coercible to the same grid via \code{terra::resample}).
#' @param dates A \code{Date} vector (or character vector coercible via
#'   \code{as.Date()}) of acquisition dates, one per raster. Must be the same
#'   length as \code{rasters}.
#' @param wavelengths A numeric vector of band centre wavelengths in nanometres.
#'   Length must equal \code{terra::nlyr(rasters[[1]])}.
#' @param bbch Optional character vector of BBCH phenological stage codes
#'   (e.g. \code{"21"}, \code{"55"}, \code{"89"}) per flight date. If
#'   \code{NULL}, stages can be added later with \code{phenostage_tag()}.
#' @param meta Optional data frame with one row per flight, containing
#'   any additional flight metadata (pilot, cloud cover, wind speed, etc.).
#' @param align Logical. If \code{TRUE} (default), rasters with slightly
#'   mismatched extents or resolutions are resampled to match the first raster.
#'   If \code{FALSE}, extents must be identical or an error is raised.
#'
#' @return A \code{FlightStack} object.
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' # Simulate two 5-band rasters at different dates
#' r1 <- rast(nrows=100, ncols=100, nlyr=5, vals=runif(100*100*5))
#' r2 <- rast(nrows=100, ncols=100, nlyr=5, vals=runif(100*100*5))
#' wl <- c(450, 550, 670, 720, 800)   # nm
#'
#' fs <- flight_stack(
#'   rasters     = list(r1, r2),
#'   dates       = c("2024-04-01", "2024-05-15"),
#'   wavelengths = wl,
#'   bbch        = c("21", "55")
#' )
#' print(fs)
#' }
#'
#' @export
flight_stack <- function(rasters,
                         dates,
                         wavelengths,
                         bbch  = NULL,
                         meta  = NULL,
                         align = TRUE) {

  # ---- Input validation ---------------------------------------------------
  if (!is.list(rasters) || length(rasters) == 0)
    cli::cli_abort("{.arg rasters} must be a non-empty list of SpatRaster objects.")

  dates <- tryCatch(as.Date(dates),
    error = function(e) cli::cli_abort("Cannot coerce {.arg dates} to Date: {e$message}"))

  if (length(dates) != length(rasters))
    cli::cli_abort(
      "{.arg dates} length ({length(dates)}) must equal {.arg rasters} length ({length(rasters)}).")

  n_bands <- as.integer(terra::nlyr(rasters[[1]]))

  if (length(wavelengths) != n_bands)
    cli::cli_abort(
      "{.arg wavelengths} length ({length(wavelengths)}) must equal band count ({n_bands}).")

  if (!is.null(bbch) && length(bbch) != length(rasters))
    cli::cli_abort(
      "{.arg bbch} length ({length(bbch)}) must equal {.arg rasters} length ({length(rasters)}).")

  if (!is.null(meta) && nrow(meta) != length(rasters))
    cli::cli_abort(
      "{.arg meta} must have one row per flight ({length(rasters)} rows required).")

  # ---- Check band counts consistent across all rasters --------------------
  band_counts <- vapply(rasters, function(r) as.integer(terra::nlyr(r)), integer(1))
  if (!all(band_counts == n_bands))
    cli::cli_abort(
      "All rasters must have the same number of bands. Found: {unique(band_counts)}.")

  # ---- Align extents if needed --------------------------------------------
  ref <- rasters[[1]]
  if (align && length(rasters) > 1) {
    rasters <- lapply(seq_along(rasters), function(i) {
      r <- rasters[[i]]
      needs_resample <- tryCatch({
        terra::compareGeom(ref, r, stopOnError = TRUE)
        FALSE
      }, error = function(e) TRUE)
      if (needs_resample) {
        cli::cli_inform("Flight {i} ({dates[i]}): resampling to match reference extent.")
        terra::resample(r, ref, method = "bilinear")
      } else {
        r
      }
    })
  }

  # ---- Sort chronologically -----------------------------------------------
  ord     <- order(dates)
  rasters <- rasters[ord]
  dates   <- dates[ord]
  if (!is.null(bbch)) bbch <- bbch[ord]
  if (!is.null(meta)) meta <- meta[ord, , drop = FALSE]

  # ---- Name rasters by date -----------------------------------------------
  names(rasters) <- as.character(dates)

  # ---- Set wavelength names on each raster --------------------------------
  wl_names <- paste0("B", round(wavelengths), "nm")
  rasters  <- lapply(rasters, function(r) { names(r) <- wl_names; r })

  # ---- Build object -------------------------------------------------------
  structure(
    list(
      rasters     = rasters,
      dates       = dates,
      wavelengths = wavelengths,
      bbch        = bbch,
      meta        = meta
    ),
    class = "FlightStack"
  )
}

# ---- S3 methods -------------------------------------------------------------

#' @export
print.FlightStack <- function(x, ...) {
  cli::cli_h1("FlightStack")
  cli::cli_bullets(c(
    "*" = "Flights    : {length(x$rasters)}",
    "*" = "Date range : {min(x$dates)} to {max(x$dates)}",
    "*" = "Bands      : {length(x$wavelengths)}  ({round(min(x$wavelengths))} - {round(max(x$wavelengths))} nm)",
    "*" = "BBCH stages: {if (is.null(x$bbch)) 'not set' else paste(x$bbch, collapse=', ')}",
    "*" = "Extent     : {paste(round(as.vector(terra::ext(x$rasters[[1]])), 2), collapse=' ')}",
    "*" = "CRS        : {terra::crs(x$rasters[[1]], describe=TRUE)$name}"
  ))
  invisible(x)
}

#' @export
summary.FlightStack <- function(object, ...) {
  df <- data.frame(
    flight     = seq_along(object$rasters),
    date       = object$dates,
    bbch       = if (is.null(object$bbch)) NA_character_ else object$bbch,
    n_bands    = vapply(object$rasters, function(r) as.integer(terra::nlyr(r)), integer(1)),
    nrow       = vapply(object$rasters, function(r) as.integer(terra::nrow(r)), integer(1)),
    ncol       = vapply(object$rasters, function(r) as.integer(terra::ncol(r)), integer(1))
  )
  print(df, row.names = FALSE)
  invisible(df)
}

#' @export
length.FlightStack <- function(x) length(x$rasters)

#' @export
`[.FlightStack` <- function(x, i) {
  flight_stack(
    rasters     = x$rasters[i],
    dates       = x$dates[i],
    wavelengths = x$wavelengths,
    bbch        = if (!is.null(x$bbch)) x$bbch[i] else NULL,
    meta        = if (!is.null(x$meta)) x$meta[i, , drop = FALSE] else NULL,
    align       = FALSE
  )
}
