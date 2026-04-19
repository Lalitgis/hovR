#' Vegetation index computation for FlightStack objects
#'
#' @description
#' Computes one or more vegetation indices across all flights in a
#' \code{FlightStack}. Uses wavelength-agnostic band matching so the same
#' formula works across sensors with different bandwidths (e.g. a 5-nm
#' Headwall vs a 10-nm Rikola). Returns a named list of \code{terra::SpatRaster}
#' objects, one per flight date.
#'
#' @section Supported built-in indices:
#' \describe{
#'   \item{NDVI}{Normalised Difference Vegetation Index — (NIR - Red) / (NIR + Red)}
#'   \item{NDRE}{Normalised Difference Red-Edge — (NIR - RedEdge) / (NIR + RedEdge)}
#'   \item{VARI}{Visible Atmospherically Resistant Index — (Green - Red) / (Green + Red - Blue)}
#'   \item{MCARI}{Modified Chlorophyll Absorption in Reflectance Index}
#'   \item{mNDVI705}{Modified NDVI at 705/750 nm — sensitive to chlorophyll in dense canopies}
#'   \item{LCI}{Leaf Chlorophyll Index — (NIR - RedEdge) / (NIR + Red)}
#'   \item{GNDVI}{Green NDVI — (NIR - Green) / (NIR + Green)}
#'   \item{CRI1}{Carotenoid Reflectance Index 1 — (1/Blue) - (1/Green)}
#'   \item{WBI}{Water Band Index — 900 nm / 970 nm}
#'   \item{PSRI}{Plant Senescence Reflectance Index — (Red - Blue) / RedEdge}
#' }

# ---- Index definitions (wavelength-based) ----------------------------------

.VI_DEFS <- list(

  NDVI = list(
    bands  = c(NIR = 800, Red = 670),
    tol    = 20,
    fun    = function(b) (b$NIR - b$Red) / (b$NIR + b$Red),
    desc   = "Normalised Difference Vegetation Index"
  ),

  NDRE = list(
    bands  = c(NIR = 800, RedEdge = 720),
    tol    = 20,
    fun    = function(b) (b$NIR - b$RedEdge) / (b$NIR + b$RedEdge),
    desc   = "Normalised Difference Red-Edge Index"
  ),

  VARI = list(
    bands  = c(Green = 550, Red = 670, Blue = 480),
    tol    = 25,
    fun    = function(b) (b$Green - b$Red) / (b$Green + b$Red - b$Blue),
    desc   = "Visible Atmospherically Resistant Index"
  ),

  MCARI = list(
    bands  = c(RE = 700, Red = 670, Green = 550),
    tol    = 20,
    fun    = function(b) ((b$RE - b$Red) - 0.2 * (b$RE - b$Green)) * (b$RE / b$Red),
    desc   = "Modified Chlorophyll Absorption Reflectance Index"
  ),

  mNDVI705 = list(
    bands  = c(B750 = 750, B705 = 705, B445 = 445),
    tol    = 15,
    fun    = function(b) (b$B750 - b$B705) / (b$B750 + b$B705 - 2 * b$B445),
    desc   = "Modified NDVI at 705/750 nm (chlorophyll in dense canopies)"
  ),

  LCI = list(
    bands  = c(NIR = 850, RedEdge = 710, Red = 670),
    tol    = 25,
    fun    = function(b) (b$NIR - b$RedEdge) / (b$NIR + b$Red),
    desc   = "Leaf Chlorophyll Index"
  ),

  GNDVI = list(
    bands  = c(NIR = 800, Green = 550),
    tol    = 20,
    fun    = function(b) (b$NIR - b$Green) / (b$NIR + b$Green),
    desc   = "Green Normalised Difference Vegetation Index"
  ),

  CRI1 = list(
    bands  = c(Blue = 510, Green = 550),
    tol    = 20,
    fun    = function(b) (1 / b$Blue) - (1 / b$Green),
    desc   = "Carotenoid Reflectance Index 1"
  ),

  WBI = list(
    bands  = c(B900 = 900, B970 = 970),
    tol    = 15,
    fun    = function(b) b$B900 / b$B970,
    desc   = "Water Band Index"
  ),

  PSRI = list(
    bands  = c(Red = 678, Blue = 500, RedEdge = 750),
    tol    = 25,
    fun    = function(b) (b$Red - b$Blue) / b$RedEdge,
    desc   = "Plant Senescence Reflectance Index"
  )
)

# ---- Band matcher -----------------------------------------------------------

#' @keywords internal
.match_band <- function(wavelengths, target_nm, tol_nm) {
  dists <- abs(wavelengths - target_nm)
  best  <- which.min(dists)
  if (dists[best] > tol_nm)
    cli::cli_abort(
      "No band within {tol_nm} nm of {target_nm} nm. \\
       Closest is {round(wavelengths[best], 1)} nm ({round(dists[best], 1)} nm away).")
  best
}

#' @keywords internal
.extract_bands <- function(raster, band_spec, wavelengths, tol) {
  lapply(stats::setNames(nm = names(band_spec)), function(nm) {
    idx <- .match_band(wavelengths, band_spec[[nm]], tol)
    raster[[idx]]
  })
}

# ---- Main function ----------------------------------------------------------

#' Compute vegetation indices across all flights
#'
#' @param stack A \code{FlightStack} object.
#' @param indices Character vector of index names. Use \code{list_vi()} to see
#'   all available built-in names. Default is \code{"NDVI"}.
#' @param custom_index Optional named list defining a custom index. Must contain:
#'   \describe{
#'     \item{\code{bands}}{Named numeric vector: band name → centre wavelength (nm)}
#'     \item{\code{tol}}{Numeric: tolerance in nm for band matching}
#'     \item{\code{fun}}{Function taking a named list of single-band rasters,
#'       returning a single-band SpatRaster}
#'     \item{\code{desc}}{Character: short description}
#'   }
#' @param na_val Numeric value to set for out-of-range pixels (default: \code{NA}).
#' @param clamp Logical. If \code{TRUE} (default for normalised indices),
#'   clamp output to \code{[-1, 1]}.
#'
#' @return A named list (by index name) of named lists (by flight date) of
#'   single-band \code{terra::SpatRaster} objects.
#'
#' @examples
#' \dontrun{
#' vi <- compute_vi(my_stack, indices = c("NDVI", "NDRE", "LCI"))
#' # Access NDVI raster for flight on 2024-05-15:
#' vi$NDVI[["2024-05-15"]]
#' }
#'
#' @seealso \code{\link{list_vi}}, \code{\link{integrate_season}},
#'   \code{\link{vi_derivative}}
#' @export
compute_vi <- function(stack,
                       indices      = "NDVI",
                       custom_index = NULL,
                       na_val       = NA_real_,
                       clamp        = TRUE) {

  if (!inherits(stack, "FlightStack"))
    cli::cli_abort("{.arg stack} must be a FlightStack object.")

  # Merge built-in and custom definitions
  defs <- .VI_DEFS
  if (!is.null(custom_index)) {
    .validate_custom_vi(custom_index)
    defs[["CUSTOM"]] <- custom_index
    if (!"CUSTOM" %in% indices) indices <- c(indices, "CUSTOM")
  }

  unknown <- setdiff(indices, names(defs))
  if (length(unknown) > 0)
    cli::cli_abort("Unknown index/indices: {unknown}. See {.fn list_vi} for options.")

  wl <- stack$wavelengths

  cli::cli_progress_bar("Computing VIs", total = length(indices) * length(stack$rasters))

  result <- lapply(stats::setNames(nm = indices), function(idx_name) {
    def <- defs[[idx_name]]
    per_date <- lapply(stats::setNames(nm = names(stack$rasters)), function(dt) {
      r    <- stack$rasters[[dt]]
      bnds <- .extract_bands(r, def$bands, wl, def$tol)
      vi   <- def$fun(bnds)
      names(vi) <- idx_name
      if (clamp && grepl("^N", idx_name))   # normalised indices → [-1,1]
        vi <- terra::clamp(vi, lower = -1, upper = 1)
      terra::NAflag(vi) <- na_val
      cli::cli_progress_update()
      vi
    })
    per_date
  })

  cli::cli_progress_done()
  structure(result, class = "VIStack", dates = stack$dates, indices = indices)
}

# ---- list_vi ----------------------------------------------------------------

#' List all available built-in vegetation indices
#'
#' @return A data frame with columns \code{name}, \code{description}, and
#'   \code{wavelengths} (comma-separated nm values).
#' @export
list_vi <- function() {
  df <- do.call(rbind, lapply(names(.VI_DEFS), function(nm) {
    def <- .VI_DEFS[[nm]]
    data.frame(
      name        = nm,
      description = def$desc,
      wavelengths = paste(sort(def$bands), collapse = ", "),
      tolerance   = def$tol,
      stringsAsFactors = FALSE
    )
  }))
  df
}

# ---- Custom index validator -------------------------------------------------

#' @keywords internal
.validate_custom_vi <- function(ci) {
  required <- c("bands", "fun", "tol", "desc")
  missing  <- setdiff(required, names(ci))
  if (length(missing) > 0)
    cli::cli_abort("custom_index missing required fields: {missing}.")
  if (!is.numeric(ci$bands) || is.null(names(ci$bands)))
    cli::cli_abort("custom_index$bands must be a named numeric vector.")
  if (!is.function(ci$fun))
    cli::cli_abort("custom_index$fun must be a function.")
}

# ---- VIStack print ----------------------------------------------------------

#' @export
print.VIStack <- function(x, ...) {
  cli::cli_h1("VIStack")
  cli::cli_bullets(c(
    "*" = "Indices : {attr(x, 'indices')}",
    "*" = "Flights : {length(attr(x, 'dates'))}",
    "*" = "Dates   : {paste(attr(x, 'dates'), collapse=', ')}"
  ))
  invisible(x)
}
