#' @importFrom methods new setClass setMethod
NULL

#' Container for cell boundary coordinates and their EFA fit
#'
#' Produced by \code{\link{compute_efa}}, which stores both the input
#' boundary matrices and every component of the fitted EFA model in a
#' single object with clearly-labeled slots.
#'
#' @slot boundaries Named list of \code{(n x 2)} coordinate matrices in the
#'   original tissue coordinate system (e.g. Xenium stage coordinates in
#'   microns), one entry per cell. Produced by
#'   \code{\link{extract_cell_boundaries}}.
#' @slot coefficients \code{(n_cells x 4 * n_harmonics)} numeric matrix of
#'   EFA coefficients (a, b, c, d per harmonic). Normalized when
#'   \code{compute_efa(normalize = TRUE)} (the default).
#' @slot raw_efa Named list of raw \code{Momocs::efourier()} output objects,
#'   one per cell. Each entry carries the boundary centroid offsets \code{a0}
#'   and \code{c0} (mean x and mean y of the original boundary points) used
#'   by \code{\link{reconstruct_cell}} to place reconstructed polygons back
#'   at the cell's original position in tissue space.
#' @slot n_harmonics Integer; number of harmonics fitted.
#' @slot valid Logical vector (length = number of cells); \code{TRUE} for
#'   cells with a successful EFA fit.
#' @export
setClass("CellEFA",
  representation(
    boundaries   = "list",
    coefficients = "matrix",
    raw_efa      = "list",
    n_harmonics  = "integer",
    valid        = "logical"
  )
)

setMethod("show", "CellEFA", function(object) {
  cat("CellEFA object\n")
  cat("  Cells     :", length(object@boundaries), "\n")
  cat("  Valid fits:", sum(object@valid), "/", length(object@valid), "\n")
  cat("  Harmonics :", object@n_harmonics, "\n")
  cat("  Coef dim  :", nrow(object@coefficients), "x",
      ncol(object@coefficients), "\n")
})
