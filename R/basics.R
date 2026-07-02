#' sfe_efa_bridge.R -- code created in dialogue with claude.ai
#' Opus 4.6 medium, June 3 2026
#'
#' Bridge between SpatialFeatureExperiment cell boundary polygons
#' and Elliptic Fourier Analysis (Momocs).
#'
#' Extracts cell segmentation polygons from an SFE object,
#' computes EFA coefficients per cell, stores them as a
#' reducedDim or colData, and provides reconstruction/visualization.
#'
#' Dependencies:
#'   BiocManager::install("SpatialFeatureExperiment")
#'   install.packages("Momocs")
#'   # optional: install.packages("pliman")

#library(SpatialFeatureExperiment)
#library(sf)
#library(Momocs)



# ============================================================
# 1. Extract ordered boundary coordinates from SFE
# ============================================================

#' Extract cell boundary coordinate matrices from an SFE object
#' @import SpatialFeatureExperiment
#' @import sf
#' @import Momocs
#'
#' @param sfe A SpatialFeatureExperiment object
#' @param geom_name Name of the colGeometry, typically "cellSeg"
#' @param n_points If not NULL, resample each boundary to this many
#'   equally-spaced points (recommended for EFA comparability)
#' @return A named list of (n x 2) matrices, one per cell
#' @export
extract_cell_boundaries <- function(sfe,
                                    geom_name = "cellSeg",
                                    n_points = 200) {
  seg_sf <- colGeometry(sfe, geom_name)
  cell_ids <- colnames(sfe)

  boundaries <- lapply(seq_len(nrow(seg_sf)), function(i) {
    geom <- st_geometry(seg_sf)[[i]]

    # Handle EMPTY geometries
    if (sf::st_is_empty(geom)) return(NULL)

    # For MULTIPOLYGON: take the largest polygon's exterior ring
    if (inherits(geom, "MULTIPOLYGON")) {
      polys <- lapply(geom, function(p) sf::st_polygon(list(p[[1]])))
      areas <- vapply(polys, function(p) as.numeric(sf::st_area(p)),
                      numeric(1))
      geom <- geom[[which.max(areas)]]
    }

    # Exterior ring is the first element of a POLYGON
    coords <- tryCatch(as.matrix(geom[[1]]), error = function(e) NULL)
    if (is.null(coords) || nrow(coords) < 4) return(NULL)

    # Remove the closing point (last == first for a closed ring)
    if (all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- coords[-nrow(coords), , drop = FALSE]
    }
    if (nrow(coords) < 3) return(NULL)
    coords
  })
  names(boundaries) <- cell_ids

  # Optional: resample to uniform point count
  if (!is.null(n_points)) {
    boundaries <- lapply(boundaries, function(coords) {
      if (is.null(coords)) return(NULL)
      resample_contour(coords, n = n_points)
    })
  }

  n_null <- sum(vapply(boundaries, is.null, logical(1)))
  if (n_null > 0) {
    message("  ", n_null, " cells have empty or degenerate geometries.")
  }

  boundaries
}

#' Resample a contour to n equally-spaced points by arc length
#'
#' @param coords (m x 2) matrix of contour coordinates
#' @param n target number of points
#' @return (n x 2) matrix
#' @export
resample_contour <- function(coords, n = 200) {
  # Close the contour temporarily for arc-length computation
  closed <- rbind(coords, coords[1, ])
  dx <- diff(closed[, 1])
  dy <- diff(closed[, 2])
  ds <- sqrt(dx^2 + dy^2)
  s <- c(0, cumsum(ds))
  total_length <- s[length(s)]

  # Target arc-length positions (exclude the closing point)
  s_target <- seq(0, total_length, length.out = n + 1)[-(n + 1)]

  # Interpolate x and y at target arc-length positions
  x_new <- stats::approx(s, closed[, 1], xout = s_target, rule = 2)$y
  y_new <- stats::approx(s, closed[, 2], xout = s_target, rule = 2)$y
  cbind(x_new, y_new)
}


# ============================================================
# 2. Compute EFA coefficients using Momocs
# ============================================================

#' Compute EFA for all cell boundaries
#'
#' @param boundaries Named list of (n x 2) coordinate matrices
#'   (output of extract_cell_boundaries)
#' @param n_harmonics Number of Fourier harmonics
#' @param normalize Logical; normalize coefficients for size,
#'   rotation, and starting-point invariance
#' @return A \code{\link{CellEFA}} object bundling the input boundaries and
#'   all fitted EFA components
#' @export
compute_efa <- function(boundaries,
                        n_harmonics = 12,
                        normalize = TRUE) {

  na_row <- rep(NA_real_, 4 * n_harmonics)

  # --- Validate each boundary before attempting EFA ---
  is_valid <- vapply(boundaries, function(coords) {
    if (is.null(coords) || !is.matrix(coords)) return(FALSE)
    if (nrow(coords) < 2 * n_harmonics + 1) return(FALSE)
    # Check for degenerate (zero-area) shapes: all collinear
    rng_x <- diff(range(coords[, 1]))
    rng_y <- diff(range(coords[, 2]))
    if (rng_x == 0 && rng_y == 0) return(FALSE)
    # Check for zero arc-length segments dominating
    dx <- diff(coords[, 1])
    dy <- diff(coords[, 2])
    ds <- sqrt(dx^2 + dy^2)
    if (sum(ds) == 0) return(FALSE)
    TRUE
  }, logical(1))

  n_invalid <- sum(!is_valid)
  if (n_invalid > 0) {
    message("  Skipping ", n_invalid, " degenerate cell boundaries ",
            "(of ", length(boundaries), " total).")
  }

  # --- Compute EFA with tryCatch per cell ---
  # Store both raw (for reconstruction) and normalized (for PCA)
  raw_efa <- vector("list", length(boundaries))
  norm_efa <- vector("list", length(boundaries))
  names(raw_efa) <- names(boundaries)
  names(norm_efa) <- names(boundaries)

  for (i in seq_along(boundaries)) {
    if (!is_valid[i]) next
    raw_efa[[i]] <- tryCatch({
      Momocs::efourier(boundaries[[i]],
                       nb.h = n_harmonics, smooth.it = 0)
    }, error = function(e) {
      NULL
    })
    if (!is.null(raw_efa[[i]]) && normalize) {
      norm_efa[[i]] <- tryCatch(
        Momocs::efourier_norm(raw_efa[[i]]),
        error = function(e) NULL
      )
    }
  }

  succeeded <- !vapply(raw_efa, is.null, logical(1))
  if (sum(!succeeded) > n_invalid) {
    message("  Additional ", sum(!succeeded) - n_invalid,
            " cells failed during efourier computation.")
  }
  message("  Successfully computed EFA for ", sum(succeeded),
          " / ", length(boundaries), " cells.")

  # --- Flatten into coefficient matrix ---
  # Use normalized coefficients for the matrix (PCA, clustering)
  # Fall back to raw if normalization was not requested
  efa_for_matrix <- if (normalize) norm_efa else raw_efa

  # efourier() returns $an,$bn,$cn,$dn
  # efourier_norm() returns $A,$B,$C,$D
  first_good <- which(succeeded)[1]
  ef1 <- efa_for_matrix[[first_good]]
  ef_names <- names(ef1)

  # Build accessor: try common Momocs conventions
  get_coefs <- function(ef) {
    if (is.null(ef)) return(na_row)
    # Try lowercase-with-n first (efourier default)
    a <- ef$an; b <- ef$bn; c <- ef$cn; d <- ef$dn
    # If NULL, try uppercase (efourier_norm convention)
    if (is.null(a)) { a <- ef$A; b <- ef$B; c <- ef$C; d <- ef$D }
    # If still NULL, try lowercase without n
    if (is.null(a)) { a <- ef$a; b <- ef$b; c <- ef$c; d <- ef$d }
    # Last resort: positional (first 4 list elements that are vectors)
    if (is.null(a)) {
      vecs <- Filter(function(x) is.numeric(x) && length(x) == n_harmonics, ef)
      if (length(vecs) >= 4) {
        a <- vecs[[1]]; b <- vecs[[2]]; c <- vecs[[3]]; d <- vecs[[4]]
      }
    }
    if (is.null(a) || length(a) != n_harmonics) return(na_row)
    c(an = a, bn = b, cn = c, dn = d)
  }

  # Report what we found for debugging
  message("  EFA result fields: ", paste(ef_names, collapse = ", "))

  coef_mat <- t(vapply(efa_for_matrix, get_coefs, numeric(4 * n_harmonics)))

  colnames(coef_mat) <- paste0(
    rep(c("a", "b", "c", "d"), each = n_harmonics),
    rep(seq_len(n_harmonics), times = 4)
  )
  rownames(coef_mat) <- names(boundaries)

  methods::new("CellEFA",
    boundaries   = boundaries,
    coefficients = coef_mat,
    raw_efa      = raw_efa,
    n_harmonics  = as.integer(n_harmonics),
    valid        = succeeded
  )
}


# ============================================================
# 3. Store EFA results back into the SFE object
# ============================================================

#' Add EFA coefficients as a reducedDim in the SFE object
#'
#' @param sfe A SpatialFeatureExperiment object
#' @param x A \code{\link{CellEFA}} object (output of compute_efa())
#' @param name Name for the reducedDim slot (default "EFA")
#' @return The modified SFE object
#' @export
store_efa_in_sfe <- function(sfe, x, name = "EFA") {
  # Store the full coefficient matrix as a reducedDim
  SingleCellExperiment::reducedDim(sfe, name) <- x@coefficients

  # Optionally store summary shape metrics in colData
  n_h <- x@n_harmonics
  coefs <- x@coefficients

  # Identify cells that have valid (non-NA) coefficients
  valid <- !is.na(coefs[, 1])

  # Harmonic power: Power_k = (a_k^2 + b_k^2 + c_k^2 + d_k^2) / 2
  power_mat <- matrix(NA_real_, nrow = nrow(coefs), ncol = n_h)
  for (k in seq_len(n_h)) {
    ak <- coefs[, paste0("a", k)]
    bk <- coefs[, paste0("b", k)]
    ck <- coefs[, paste0("c", k)]
    dk <- coefs[, paste0("d", k)]
    power_mat[, k] <- (ak^2 + bk^2 + ck^2 + dk^2) / 2
  }
  total_power <- rowSums(power_mat, na.rm = TRUE)
  total_power[!valid] <- NA_real_

  # Shape complexity: number of harmonics to reach 99% power
  complexity <- rep(NA_integer_, nrow(coefs))
  ellipticity <- rep(NA_real_, nrow(coefs))

  if (any(valid)) {
    cum_power <- t(apply(power_mat[valid, , drop = FALSE], 1, cumsum))
    complexity[valid] <- apply(cum_power, 1, function(cp) {
      tp <- cp[length(cp)]
      if (tp == 0) return(1L)
      min(which(cp >= 0.99 * tp))
    })
    ellipticity[valid] <- power_mat[valid, 1] / total_power[valid]
    ellipticity[valid & is.nan(ellipticity)] <- 1
  }

  sfe$efa_complexity <- complexity
  sfe$efa_ellipticity <- ellipticity
  message("  ", sum(!valid), " cells have NA shape descriptors ",
          "(degenerate boundaries).")

  sfe
}


# ============================================================
# 4. Reconstruction and visualization
# ============================================================

#' Remap EFA field names to the convention efourier_i expects
#' @param ef A list returned by efourier or efourier_norm
#' @return The same list with fields named an, bn, cn, dn, a0, c0
#' @export
remap_efa_fields <- function(ef) {
  # If already has $an, nothing to do

  if (!is.null(ef$an)) return(ef)
  # Normalized output uses A, B, C, D, ao, co
  if (!is.null(ef$A)) {
    ef$an <- ef$A; ef$bn <- ef$B; ef$cn <- ef$C; ef$dn <- ef$D
    ef$a0 <- ef$ao; ef$c0 <- ef$co
  }
  ef
}

#' Reconstruct a cell boundary from its EFA coefficients
#'
#' @param x A \code{\link{CellEFA}} object (output of compute_efa())
#' @param cell_id Character; which cell to reconstruct
#' @param n_harmonics_use How many harmonics to use in reconstruction
#'   (NULL = all available)
#' @param n_points Number of points in the reconstructed contour
#' @return (n_points x 2) matrix of reconstructed coordinates in the original
#'   tissue coordinate system
#' @export
reconstruct_cell <- function(x, cell_id,
                             n_harmonics_use = NULL,
                             n_points = 300) {
  ef <- x@raw_efa[[cell_id]]
  ef <- remap_efa_fields(ef)
  if (is.null(n_harmonics_use)) {
    n_harmonics_use <- x@n_harmonics
  }
  Momocs::efourier_i(ef, nb.h = n_harmonics_use, nb.pts = n_points)
}


# ============================================================
# 4b. Per-cell goodness-of-fit
# ============================================================

#' Compute nearest-neighbor distances from each row of A to B
#' @param A,B (n x 2) and (m x 2) coordinate matrices
#' @return numeric vector of length nrow(A)
nn_dists <- function(A, B) {
  # For each point in A, find min Euclidean distance to any point in B
  vapply(seq_len(nrow(A)), function(i) {
    dx <- B[, 1] - A[i, 1]
    dy <- B[, 2] - A[i, 2]
    sqrt(min(dx^2 + dy^2))
  }, numeric(1))
}

#' Signed area of a polygon (positive if CCW)
#' @param coords (n x 2) matrix, not closed
polygon_area <- function(coords) {
  n <- nrow(coords)
  x <- coords[, 1]; y <- coords[, 2]
  # Shoelace formula
  j <- c(2:n, 1)
  abs(sum(x * y[j] - x[j] * y)) / 2
}

#' Perimeter of a polygon
#' @param coords (n x 2) matrix, not closed
polygon_perimeter <- function(coords) {
  closed <- rbind(coords, coords[1, ])
  sum(sqrt(diff(closed[, 1])^2 + diff(closed[, 2])^2))
}

#' Compute per-cell EFA goodness-of-fit metrics
#'
#' For each cell, reconstructs the boundary from the stored EFA
#' coefficients and compares to the original boundary.
#'
#' @param x A \code{\link{CellEFA}} object (output of compute_efa())
#' @param n_harmonics_use Number of harmonics for reconstruction
#'   (NULL = all available)
#' @return A data.frame with one row per cell and columns:
#'   - mean_dev: mean distance from original to nearest reconstructed point
#'   - max_dev: max such distance (one-sided Hausdorff)
#'   - symmetric_hausdorff: max of both one-sided Hausdorff distances
#'   - rel_mean_dev: mean_dev / perimeter (dimensionless)
#'   - area_ratio: area(reconstruction) / area(original)
#'   - power_captured: fraction of total harmonic power in first
#'       n_harmonics_use harmonics (1.0 if using all)
#'
efa_goodness_of_fit <- function(x,
                                 n_harmonics_use = NULL) {
  cell_ids <- names(x@boundaries)
  n_cells <- length(cell_ids)
  if (is.null(n_harmonics_use)) {
    n_harmonics_use <- x@n_harmonics
  }

  mean_dev <- max_dev <- sym_hausdorff <- rep(NA_real_, n_cells)
  rel_mean_dev <- area_ratio <- power_frac <- rep(NA_real_, n_cells)
  names(mean_dev) <- names(max_dev) <- names(sym_hausdorff) <- cell_ids
  names(rel_mean_dev) <- names(area_ratio) <- names(power_frac) <- cell_ids

  for (i in seq_len(n_cells)) {
    cid <- cell_ids[i]
    orig <- x@boundaries[[cid]]
    ef <- x@raw_efa[[cid]]
    if (is.null(orig) || is.null(ef)) next

    # Reconstruct
    recon <- tryCatch(
      reconstruct_cell(x, cid,
                       n_harmonics_use = n_harmonics_use,
                       n_points = nrow(orig)),
      error = function(e) NULL
    )
    if (is.null(recon)) next

    # Nearest-neighbor distances: original -> reconstruction
    d_orig_to_recon <- nn_dists(orig, recon)
    d_recon_to_orig <- nn_dists(recon, orig)

    mean_dev[i] <- mean(d_orig_to_recon)
    max_dev[i] <- max(d_orig_to_recon)
    sym_hausdorff[i] <- max(max(d_orig_to_recon), max(d_recon_to_orig))

    perim <- polygon_perimeter(orig)
    rel_mean_dev[i] <- if (perim > 0) mean_dev[i] / perim else NA_real_

    area_orig <- polygon_area(orig)
    area_recon <- polygon_area(recon)
    area_ratio[i] <- if (area_orig > 0) area_recon / area_orig else NA_real_

    # Power fraction: how much power is in the first n_harmonics_use
    # vs. all n_harmonics available
    a <- ef$an; b <- ef$bn; cc <- ef$cn; d <- ef$dn
    if (!is.null(a)) {
      full_power <- sum((a^2 + b^2 + cc^2 + d^2) / 2)
      used_power <- sum((a[1:n_harmonics_use]^2 +
                         b[1:n_harmonics_use]^2 +
                         cc[1:n_harmonics_use]^2 +
                         d[1:n_harmonics_use]^2) / 2)
      power_frac[i] <- if (full_power > 0) used_power / full_power else 1
    }
  }

  data.frame(
    cell_id = cell_ids,
    mean_dev = mean_dev,
    max_dev = max_dev,
    symmetric_hausdorff = sym_hausdorff,
    rel_mean_dev = rel_mean_dev,
    area_ratio = area_ratio,
    power_captured = power_frac,
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )
}

#' Store GoF metrics in colData of an SFE object
#'
#' @param sfe SpatialFeatureExperiment object
#' @param gof_df Output of efa_goodness_of_fit()
#' @param prefix Column name prefix (default "efa_gof_")
#' @return Modified SFE object
store_gof_in_sfe <- function(sfe, gof_df, prefix = "efa_gof_") {
  metrics <- c("mean_dev", "max_dev", "symmetric_hausdorff",
               "rel_mean_dev", "area_ratio", "power_captured")
  for (m in metrics) {
    colData(sfe)[[paste0(prefix, m)]] <- gof_df[colnames(sfe), m]
  }
  sfe
}

#' Plot original vs reconstructed cell boundaries
#'
#' @param x A \code{\link{CellEFA}} object (output of compute_efa())
#' @param cell_ids Character vector of cell IDs to plot
#' @param n_harmonics_seq Harmonic counts to show progressive reconstruction
#' @param ncol Number of columns in the plot layout
#' @export
plot_efa_reconstruction <- function(x,
                                    cell_ids = NULL,
                                    n_harmonics_seq = c(1, 3, 6, 12),
                                    ncol = length(n_harmonics_seq)) {
  if (is.null(cell_ids)) {
    cell_ids <- names(x@boundaries)[seq_len(min(4L, length(x@boundaries)))]
  }

  n_cells <- length(cell_ids)
  n_steps <- length(n_harmonics_seq)
  par(mfrow = c(n_cells, n_steps + 1), mar = c(1, 1, 2, 1))

  for (cid in cell_ids) {
    orig <- x@boundaries[[cid]]
    # Plot original
    plot(orig, type = "l", asp = 1, axes = FALSE,
         main = paste0(cid, "\noriginal"), cex.main = 0.8)
    polygon(orig, border = "black", col = adjustcolor("grey80", 0.3))

    # Progressive reconstructions
    for (nh in n_harmonics_seq) {
      recon <- reconstruct_cell(x, cid, n_harmonics_use = nh)
      plot(orig, type = "n", asp = 1, axes = FALSE,
           main = paste0(nh, " harmonics"), cex.main = 0.8)
      polygon(orig, border = "grey70", col = NA, lty = 2)
      lines(rbind(recon, recon[1L, ]), col = "firebrick", lwd = 2)
    }
  }
}


# ============================================================
# 4c. Reconstruct boundaries as sf geometries in original space
# ============================================================

#' Convert EFA reconstructions to sf polygons in original coordinate space
#'
#' The reconstruction uses the DC terms (a0/c0) stored in the raw EFA output,
#' so the returned polygons are in the same coordinate system as the original
#' cell segmentation geometries (e.g. Xenium stage coordinates in microns).
#'
#' @param x A \code{\link{CellEFA}} object (output of compute_efa())
#' @param cell_ids Character vector of cell IDs to reconstruct (NULL = all valid)
#' @param n_harmonics_use Number of harmonics (NULL = all available)
#' @param n_points Number of points in the reconstructed contour
#' @return An sf data.frame with columns cell_id and geometry (POLYGON)
#' @examples
#' make_ellipse <- function(cx, cy, a, b, n = 200L) {
#'   theta <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)]
#'   cbind(x = cx + a * cos(theta), y = cy + b * sin(theta))
#' }
#' boundaries <- list(
#'   cellA = make_ellipse(100, 200, 15, 8),
#'   cellB = make_ellipse(150, 220, 10, 10),
#'   cellC = make_ellipse(120, 260, 12,  6)
#' )
#' x <- compute_efa(boundaries, n_harmonics = 6)
#' efa_to_sf(x)
#' @export
efa_to_sf <- function(x, cell_ids = NULL,
                      n_harmonics_use = NULL,
                      n_points = 300) {
  valid_ids <- names(x@raw_efa)[x@valid]
  if (is.null(cell_ids)) {
    cell_ids <- valid_ids
  } else {
    bad <- setdiff(cell_ids, valid_ids)
    if (length(bad) > 0) {
      message("  ", length(bad), " requested cell(s) have no valid EFA; skipping.")
      cell_ids <- intersect(cell_ids, valid_ids)
    }
  }

  polys <- lapply(cell_ids, function(cid) {
    coords <- tryCatch(
      reconstruct_cell(x, cid, n_harmonics_use, n_points),
      error = function(e) NULL
    )
    if (is.null(coords)) return(sf::st_polygon(list()))
    sf::st_polygon(list(rbind(coords, coords[1L, ])))
  })

  sf::st_sf(
    cell_id = cell_ids,
    geometry = sf::st_sfc(polys)
  )
}

#' Plot EFA-reconstructed and original boundaries in original coordinate space
#'
#' Overlays original cell boundaries (grey) and their EFA reconstructions
#' (red) in a single base-R plot using tissue/stage coordinates so spatial
#' relationships between cells are preserved.
#'
#' @param x A \code{\link{CellEFA}} object (output of compute_efa())
#' @param cell_ids Cell IDs to include (NULL = all valid, up to max_cells)
#' @param max_cells Maximum number of cells to plot when cell_ids is NULL
#' @param n_harmonics_use Number of harmonics (NULL = all available)
#' @param n_points Points in the reconstructed contour
#' @param col_orig Color for original boundaries
#' @param col_recon Color for reconstructed boundaries
#' @param lwd_recon Line width for reconstructed boundaries
#' @param ... Additional arguments passed to plot()
#' @return Invisibly returns an sf object of the reconstructed boundaries
#' @examples
#' make_ellipse <- function(cx, cy, a, b, n = 200L) {
#'   theta <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)]
#'   cbind(x = cx + a * cos(theta), y = cy + b * sin(theta))
#' }
#' boundaries <- list(
#'   cellA = make_ellipse(100, 200, 15, 8),
#'   cellB = make_ellipse(150, 220, 10, 10),
#'   cellC = make_ellipse(120, 260, 12,  6)
#' )
#' x <- compute_efa(boundaries, n_harmonics = 6)
#' plot_efa_in_space(x)
#' @export
plot_efa_in_space <- function(x,
                              cell_ids = NULL,
                              max_cells = 200L,
                              n_harmonics_use = NULL,
                              n_points = 300L,
                              col_orig = "grey50",
                              col_recon = "firebrick",
                              lwd_recon = 1.5,
                              ...) {
  valid_ids <- names(x@raw_efa)[x@valid]
  if (is.null(cell_ids)) {
    cell_ids <- valid_ids[seq_len(min(max_cells, length(valid_ids)))]
  } else {
    cell_ids <- intersect(cell_ids, valid_ids)
  }

  # Set axis limits from the union of original boundary extents
  orig_list <- Filter(Negate(is.null), x@boundaries[cell_ids])
  if (length(orig_list) == 0L) stop("No valid boundaries to plot.")
  all_orig <- do.call(rbind, orig_list)
  xlim <- range(all_orig[, 1L])
  ylim <- range(all_orig[, 2L])

  n_h <- if (is.null(n_harmonics_use)) x@n_harmonics else n_harmonics_use
  plot(NULL, xlim = xlim, ylim = ylim, asp = 1,
       xlab = "x", ylab = "y",
       main = sprintf("EFA reconstructions (%d harmonics, %d cells)",
                      n_h, length(cell_ids)),
       ...)

  for (cid in cell_ids) {
    orig <- x@boundaries[[cid]]
    if (!is.null(orig)) {
      polygon(orig[c(seq_len(nrow(orig)), 1L), ], border = col_orig,
              col = NA, lwd = 1)
    }
    recon <- tryCatch(
      reconstruct_cell(x, cid, n_harmonics_use, n_points),
      error = function(e) NULL
    )
    if (!is.null(recon)) {
      lines(recon[c(seq_len(nrow(recon)), 1L), ], col = col_recon,
            lwd = lwd_recon)
    }
  }

  legend("topright",
         legend = c("original", sprintf("EFA (%d harmonics)", n_h)),
         col    = c(col_orig, col_recon),
         lty    = 1L,
         lwd    = c(1, lwd_recon),
         bty    = "n")

  invisible(efa_to_sf(x, cell_ids, n_harmonics_use, n_points))
}


#' Interactive spatial plot of EFA reconstructions via ggplot2 and plotly
#'
#' Renders original and EFA-reconstructed cell boundaries as an interactive
#' plotly figure. Pan and zoom into any subregion to inspect fit quality for
#' individual cells; hover to see the cell ID. Requires the ggplot2 and plotly
#' packages (listed under Suggests).
#'
#' @param x A \code{\link{CellEFA}} object (output of compute_efa())
#' @param cell_ids Cell IDs to include (NULL = all valid, up to max_cells)
#' @param max_cells Maximum number of cells when cell_ids is NULL
#' @param n_harmonics_use Number of harmonics (NULL = all available)
#' @param n_points Points in the reconstructed contour
#' @param col_orig Color for original boundaries
#' @param col_recon Color for EFA-reconstructed boundaries
#' @return A plotly figure object
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("plotly",  quietly = TRUE)) {
#'   make_ellipse <- function(cx, cy, a, b, n = 200L) {
#'     theta <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)]
#'     cbind(x = cx + a * cos(theta), y = cy + b * sin(theta))
#'   }
#'   boundaries <- list(
#'     cellA = make_ellipse(100, 200, 15,  8),
#'     cellB = make_ellipse(150, 220, 10, 10),
#'     cellC = make_ellipse(120, 260, 12,  6)
#'   )
#'   x <- compute_efa(boundaries, n_harmonics = 3)
#'   plot_efa_interactive(x)
#' }
#' @export
plot_efa_interactive <- function(x,
                                 cell_ids      = NULL,
                                 max_cells     = 500L,
                                 n_harmonics_use = NULL,
                                 n_points      = 300L,
                                 col_orig      = "steelblue",
                                 col_recon     = "firebrick") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required: install.packages('ggplot2')")
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("plotly is required: install.packages('plotly')")

  valid_ids <- names(x@raw_efa)[x@valid]
  if (is.null(cell_ids)) {
    cell_ids <- valid_ids[seq_len(min(max_cells, length(valid_ids)))]
  } else {
    cell_ids <- intersect(cell_ids, valid_ids)
  }

  make_rows <- function(cid, mat, type) {
    if (is.null(mat)) return(NULL)
    data.frame(x = mat[, 1L], y = mat[, 2L],
               cell_id = cid, type = type,
               stringsAsFactors = FALSE)
  }

  orig_rows  <- lapply(cell_ids, function(cid)
    make_rows(cid, x@boundaries[[cid]], "original"))
  recon_rows <- lapply(cell_ids, function(cid) {
    coords <- tryCatch(
      reconstruct_cell(x, cid, n_harmonics_use, n_points),
      error = function(e) NULL
    )
    make_rows(cid, coords, "EFA")
  })

  plot_df       <- do.call(rbind, c(orig_rows, recon_rows))
  plot_df$type  <- factor(plot_df$type, levels = c("original", "EFA"))
  plot_df$group <- paste(plot_df$cell_id, plot_df$type, sep = ".")

  n_h <- if (is.null(n_harmonics_use)) efa_result$n_harmonics else n_harmonics_use

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x     = x,
      y     = y,
      group = group,
      color = type,
      text  = cell_id
    )
  ) +
    ggplot2::geom_polygon(fill = NA, linewidth = 0.4) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_manual(
      values = c(original = col_orig, EFA = col_recon),
      name   = NULL
    ) +
    ggplot2::labs(
      title = sprintf("EFA reconstructions (%d harmonics, %d cells)",
                      n_h, length(cell_ids)),
      x = "x", y = "y"
    ) +
    ggplot2::theme_minimal()

  plotly::ggplotly(p, tooltip = "text")
}


# ============================================================
# 5. PCA on EFA shape space
# ============================================================

#' PCA on EFA coefficients and optionally store in SFE
#'
#' @param sfe SFE object with EFA reducedDim already stored
#' @param efa_dim_name Name of the EFA reducedDim
#' @param n_pcs Number of PCs to retain
#' @param store_name Name for the PCA reducedDim (NULL = don't store)
#' @return prcomp result, invisibly
#' @export
efa_pca <- function(sfe,
                    efa_dim_name = "EFA",
                    n_pcs = 10,
                    store_name = "EFA_PCA") {
  coef_mat <- SingleCellExperiment::reducedDim(sfe, efa_dim_name)

  # If normalized, drop constant coefficients
  # (a1=1, b1=0, c1=0 after normalization)
  if (all(coef_mat[, "a1"] == 1) && all(coef_mat[, "b1"] == 0)) {
    coef_mat <- coef_mat[, -match(c("a1", "b1", "c1"), colnames(coef_mat))]
  }

  pca_res <- stats::prcomp(coef_mat, center = TRUE, scale. = TRUE)

  if (!is.null(store_name)) {
    n_pcs <- min(n_pcs, ncol(pca_res$x))
    SingleCellExperiment::reducedDim(sfe, store_name) <- pca_res$x[, 1:n_pcs]
  }

  # Return both the modified SFE and the prcomp result
  list(sfe = sfe, pca = pca_res)
}


# ============================================================
# 6. Convenience wrapper: full pipeline
# ============================================================

#' Run the full EFA pipeline on an SFE object
#'
#' @param sfe SpatialFeatureExperiment object with cell segmentation
#' @param geom_name Geometry name (default "cellSeg")
#' @param n_points Points to resample each boundary to
#' @param n_harmonics Number of Fourier harmonics
#' @param normalize Normalize coefficients
#' @param do_pca Run PCA on the shape space
#' @param n_pcs Number of PCs to retain
#' @param compute_gof Logical; compute goodness-of-fit metrics (default TRUE)
#' @return Modified SFE object with EFA and optionally EFA_PCA
#'   in reducedDims, and efa_complexity / efa_ellipticity in colData
#' @export
run_efa_pipeline <- function(sfe,
                             geom_name = "cellSeg",
                             n_points = 200,
                             n_harmonics = 12,
                             normalize = TRUE,
                             do_pca = TRUE,
                             n_pcs = 10,
                             compute_gof = TRUE) {
  message("Extracting cell boundaries...")
  boundaries <- extract_cell_boundaries(sfe, geom_name, n_points)

  message("Computing EFA (", n_harmonics, " harmonics)...")
  x <- compute_efa(boundaries, n_harmonics, normalize)

  message("Storing in SFE...")
  sfe <- store_efa_in_sfe(sfe, x)

  if (compute_gof) {
    message("Computing per-cell goodness-of-fit...")
    gof_df <- efa_goodness_of_fit(x)
    sfe <- store_gof_in_sfe(sfe, gof_df)
    message("  Median rel. mean deviation: ",
            round(median(gof_df$rel_mean_dev, na.rm = TRUE), 5))
    message("  Median area ratio: ",
            round(median(gof_df$area_ratio, na.rm = TRUE), 4))
  }

  if (do_pca) {
    message("Running PCA on shape space...")
    pca_out <- efa_pca(sfe, n_pcs = n_pcs)
    sfe <- pca_out$sfe
  }

  message("Done. New reducedDims: 'EFA'",
          if (do_pca) ", 'EFA_PCA'", ".")
  message("New colData columns: 'efa_complexity', 'efa_ellipticity'",
          if (compute_gof) ", 'efa_gof_*'", ".")

  sfe
}

