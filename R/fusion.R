# fusion.R -- Morphological fusion and Lack-of-Fit utilities
# Created on July 1, 2026

#' Morphologically fuse disconnected components in a binary mask
#'
#' Dilates a binary mask to bridge small gaps between components, maps the
#' dilated components back to the original mask to group connected pieces,
#' and takes their union to produce clean, fused geometries.
#'
#' @param mask A binary matrix or 2D logical matrix.
#' @param brush_size Integer; size of the disk brush used for dilation (must be odd).
#' @return A matrix of the same dimensions as \code{mask} with fused and labeled components.
#' @import EBImage
#' @export
#' @examples
#' # Create a small binary matrix with a gap
#' m <- matrix(0, nrow = 20, ncol = 20)
#' m[3:7, 5] <- 1
#' m[9:13, 5] <- 1  # Gapped vertical line
#' fused <- fuse_components(m, brush_size = 3)
#' print(fused[3:13, 5])
fuse_components <- function(mask, brush_size = 9) {
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("EBImage is required for morphological operations.")
  }
  
  # Convert mask to logical matrix
  mask_logical <- as.matrix(mask) > 0
  
  # Label original components
  labeled_orig <- EBImage::bwlabel(EBImage::Image(mask_logical))
  orig_mat <- labeled_orig@.Data
  
  # Dilate to bridge gaps
  dilated <- EBImage::dilate(EBImage::Image(mask_logical), EBImage::makeBrush(brush_size, shape = 'disc'))
  labeled_dilated <- EBImage::bwlabel(dilated)
  dil_mat <- labeled_dilated@.Data
  
  # Create a container for the fused mask
  fused_mat <- matrix(0L, nrow = nrow(mask), ncol = ncol(mask))
  
  # Unique non-zero dilated IDs
  dil_ids <- unique(as.vector(dil_mat))
  dil_ids <- dil_ids[dil_ids != 0]
  
  for (dil_id in dil_ids) {
    # Pixels belonging to this dilated component
    dil_pixels <- which(dil_mat == dil_id, arr.ind = TRUE)
    
    # Original component IDs covered by this dilated component
    orig_ids <- unique(orig_mat[dil_pixels])
    orig_ids <- orig_ids[orig_ids != 0]
    
    if (length(orig_ids) > 0) {
      # Take union of original components
      fused_mask_vec <- orig_mat %in% orig_ids
      fused_mask <- matrix(fused_mask_vec, nrow = nrow(orig_mat), ncol = ncol(orig_mat))
      
      # Fill internal cavities/holes
      fused_filled <- EBImage::fillHull(EBImage::Image(fused_mask))
      
      # Merge into the final matrix
      fused_mat[fused_filled@.Data > 0] <- as.integer(dil_id)
    }
  }
  fused_mat
}

#' Compute Jaccard Lack-of-Fit between two shapes
#'
#' Computes the Lack-of-Fit (1 - IoU / Jaccard Index) between two shape geometries.
#' The shapes can be provided as coordinate matrices or sf geometry objects.
#'
#' @param shape1 A coordinate matrix (n x 2) or an sf polygon object.
#' @param shape2 A coordinate matrix (m x 2) or an sf polygon object.
#' @return A numeric value representing the Lack-of-Fit (between 0 and 1).
#' @import sf
#' @export
#' @examples
#' library(sf)
#' # Create two overlapping squares
#' s1 <- st_polygon(list(matrix(c(0,0, 4,0, 4,4, 0,4, 0,0), ncol=2, byrow=TRUE)))
#' s2 <- st_polygon(list(matrix(c(2,0, 6,0, 6,4, 2,4, 2,0), ncol=2, byrow=TRUE)))
#' lof <- compute_jaccard_lof(s1, s2)
#' print(lof)  # Intersection is 8, Union is 24, IoU is 1/3, LoF is 2/3
compute_jaccard_lof <- function(shape1, shape2) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("sf package is required for Jaccard calculation.")
  }
  
  to_sf_poly <- function(shape) {
    if (inherits(shape, "sfg") || inherits(shape, "sfc") || inherits(shape, "sf")) {
      return(sf::st_geometry(shape)[[1]])
    }
    if (is.matrix(shape) || is.data.frame(shape)) {
      coords <- as.matrix(shape)
      if (nrow(coords) > 2 && !all(coords[1, ] == coords[nrow(coords), ])) {
        coords <- rbind(coords, coords[1, ])
      }
      return(sf::st_polygon(list(coords)))
    }
    stop("Unsupported shape format.")
  }
  
  poly1 <- to_sf_poly(shape1)
  poly2 <- to_sf_poly(shape2)
  
  iou <- tryCatch({
    inter_area <- sf::st_area(sf::st_intersection(poly1, poly2))
    union_area <- sf::st_area(sf::st_union(poly1, poly2))
    if (length(inter_area) == 0) inter_area <- 0
    if (length(union_area) == 0 || union_area == 0) union_area <- 1
    as.numeric(inter_area / union_area)
  }, error = function(e) {
    NA_real_
  })
  
  1 - iou
}

#' Find optimal translation shift between two point sets
#'
#' Computes the optimal integer translation shift (drow, dcol) that maps a set
#' of source points to a set of target points by minimizing the average distance to the
#' nearest target point.
#'
#' @param src_pts A matrix or data.frame of source coordinates (n x 2).
#' @param target_pts A matrix or data.frame of target coordinates (m x 2).
#' @param search_range A list of two integer vectors specifying the search range for rows and cols,
#'   e.g., \code{list(row = -20:20, col = -20:20)}.
#' @return A list with components:
#'   - shift: A numeric vector c(drow, dcol) of the optimal shift.
#'   - mean_dist: The minimum mean distance achieved.
#' @export
#' @examples
#' src <- matrix(c(10, 10, 20, 20), ncol = 2, byrow = TRUE)
#' target <- matrix(c(12, 11, 22, 21), ncol = 2, byrow = TRUE)
#' res <- find_optimal_translation(src, target, list(row = 0:5, col = 0:5))
#' print(res$shift)
find_optimal_translation <- function(src_pts, target_pts, search_range = list(row = -30:30, col = -30:30)) {
  src <- as.matrix(src_pts)
  tgt <- as.matrix(target_pts)
  
  best_shift <- c(0, 0)
  min_mean_dist <- Inf
  
  for (dr in search_range$row) {
    for (dc in search_range$col) {
      shifted_src <- src
      shifted_src[, 1] <- shifted_src[, 1] + dr
      shifted_src[, 2] <- shifted_src[, 2] + dc
      
      # Compute mean distance to nearest target point
      dists <- vapply(seq_len(nrow(shifted_src)), function(i) {
        dx <- tgt[, 1] - shifted_src[i, 1]
        dy <- tgt[, 2] - shifted_src[i, 2]
        sqrt(min(dx^2 + dy^2))
      }, numeric(1))
      
      mean_d <- mean(dists)
      if (mean_d < min_mean_dist) {
        min_mean_dist <- mean_d
        best_shift <- c(row = dr, col = dc)
      }
    }
  }
  
  list(shift = best_shift, mean_dist = min_mean_dist)
}
