# This file has Claude-generated code for cell shape analysis

# ============================================================
# Utility: get listw from colGraph (handles both nb and listw)
# ============================================================

#' Safely extract a listw object from an SFE colGraph
#'
#' SFE's colGraph() typically returns a listw directly, but this
#' helper handles both cases to avoid double-conversion errors.
#'
#' @param sfe SpatialFeatureExperiment object
#' @param graph_name Name of the colGraph
#' @return A spdep listw object
#' @export
get_listw <- function(sfe, graph_name = "knn") {
  g <- colGraph(sfe, graph_name)
  if (inherits(g, "listw")) return(g)
  if (inherits(g, "nb")) return(spdep::nb2listw(g, style = "W"))
  stop("Unexpected colGraph class: ", paste(class(g), collapse = ", "))
}

#' Extract the nb (neighbors) list from a colGraph
#'
#' @param sfe SpatialFeatureExperiment object
#' @param graph_name Name of the colGraph
#' @return A spdep nb object
#' @export
get_nb <- function(sfe, graph_name = "knn") {
  g <- colGraph(sfe, graph_name)
  if (inherits(g, "listw")) return(g$neighbours)
  if (inherits(g, "nb")) return(g)
  stop("Unexpected colGraph class: ", paste(class(g), collapse = ", "))
}


#' Build a spatial neighbor graph and compute Moran's I for shape features
#'
#' Tests whether cell morphology is spatially organized: do nearby
#' cells tend to have similar shapes? Uses Voyager's spatial
#' statistics infrastructure.
#'
#' @param sfe SFE object with EFA_PCA in reducedDims and shape
#'   metrics in colData (output of run_efa_pipeline)
#' @param graph_name Name for the colGraph (default "knn")
#' @param k Number of nearest neighbors
#' @param n_pcs Number of shape PCs to test
#' @param features Additional colData columns to test
#'   (default: efa_ellipticity, efa_complexity)
#' @return Modified SFE object with colGraph and Moran's I results
#' @export
spatial_shape_analysis <- function(sfe,
                                   graph_name = "knn",
                                   k = 6,
                                   n_pcs = 5,
                                   features = c("efa_ellipticity",
                                                "efa_complexity")) {
  if (!requireNamespace("Voyager", quietly = TRUE)) {
    stop("Voyager is required: BiocManager::install('Voyager')")
  }

  # --- Build spatial neighbor graph if not present ---
  if (!(graph_name %in% colGraphNames(sfe))) {
    message("Building k=", k, " nearest neighbor graph...")
    colGraph(sfe, graph_name) <- SpatialFeatureExperiment::findSpatialNeighbors(
      sfe, method = "knearneigh", k = k
    )
  }

  # --- Store shape PCs as colData for Moran's I ---
  efa_pca_mat <- SingleCellExperiment::reducedDim(sfe, "EFA_PCA")
  n_pcs <- min(n_pcs, ncol(efa_pca_mat))
  pc_names <- paste0("shape_PC", seq_len(n_pcs))
  for (j in seq_len(n_pcs)) {
    sfe[[pc_names[j]]] <- efa_pca_mat[, j]
  }

  # --- Global Moran's I on shape features ---
  all_features <- c(pc_names, intersect(features, names(colData(sfe))))
  message("Computing Moran's I for: ",
          paste(all_features, collapse = ", "))
  sfe <- Voyager::colDataMoransI(sfe, all_features,
                                  colGraphName = graph_name)

  # Extract and report
  morans <- vapply(all_features, function(f) {
    attr(sfe[[f]], "moran_I") %||% NA_real_
  }, numeric(1))

  # Moran's I may be stored in colFeatureData instead
  # Try to retrieve from colFeatureData if attr didn't work
  cfd <- tryCatch(Voyager::colFeatureData(sfe), error = function(e) NULL)
  if (!is.null(cfd) && "moran_sample" %in% names(cfd)) {
    for (f in all_features) {
      if (f %in% rownames(cfd)) {
        morans[f] <- cfd[f, "moran_sample"]
      }
    }
  }

  message("\nGlobal Moran's I for shape features:")
  for (f in all_features) {
    val <- if (is.na(morans[f])) "see colFeatureData(sfe)" else
      round(morans[f], 4)
    message("  ", f, ": ", val)
  }

  sfe
}

#' Compute shape heterogeneity in spatial neighborhoods
#'
#' For each cell, computes the variance of a shape metric
#' among its spatial neighbors. High local heterogeneity
#' suggests tissue boundaries or mixed microenvironments.
#'
#' @param sfe SFE object with colGraph
#' @param feature Shape feature to assess (default "efa_ellipticity")
#' @param graph_name Which colGraph to use
#' @return Modified SFE with neighborhood heterogeneity in colData
#' @export
neighborhood_shape_heterogeneity <- function (sfe, feature = "efa_ellipticity", graph_name = "knn") 
{
    nb <- get_nb(sfe, graph_name)
    x <- sfe[[feature]]
    if (is.null(x)) 
        stop("Feature '", feature, "' not found in colData.")
    local_sd <- vapply(seq_along(nb), function(i) {
        neighbors <- nb[[i]]
        if (length(neighbors) == 0 || (length(neighbors) == 1 && 
            neighbors[1] == 0L)) {
            return(NA_real_)
        }
        sd(x[neighbors], na.rm = TRUE)
    }, numeric(1))
    col_name <- paste0("nbhd_sd_", feature)
    sfe[[col_name]] <- local_sd
    message("Stored neighborhood heterogeneity as '", col_name, 
        "'.")
    message("  Median: ", round(median(local_sd, na.rm = TRUE), 
        4), "  Max: ", round(max(local_sd, na.rm = TRUE), 4))
    sfe
}


#' Compute local indicators of spatial association (LISA) for shape
#'
#' Identifies spatial clusters of morphologically similar cells.
#' Cells are classified as High-High, Low-Low (clusters),
#' High-Low, Low-High (outliers), or Not significant.
#' 
#' @param sfe SFE object (after spatial_shape_analysis)
#' @param feature Which shape feature to analyze (default "shape_PC1")
#' @param graph_name Which colGraph to use
#' @param nsim Number of Monte Carlo simulations for p-values
#' @param p_threshold Significance threshold
#' @return Modified SFE with LISA classifications in colData
#' @export   
shape_lisa = function (sfe, feature = "shape_PC1", graph_name = "knn", nsim = 499, 
    p_threshold = 0.05) 
{
    if (!requireNamespace("spdep", quietly = TRUE)) {
        stop("spdep is required: install.packages('spdep')")
    }

# get spatial weights
    listw <- get_listw(sfe, graph_name)
    x <- sfe[[feature]]
    if (is.null(x)) 
        stop("Feature '", feature, "' not found in colData.")
# handle NAs
    valid <- !is.na(x)
    if (any(!valid)) {
        message("  Excluding ", sum(!valid), " cells with NA values.")
    }
    message("Computing local Moran's I for '", feature, "' (", 
        nsim, " permutations)...")
    lmoran <- spdep::localmoran_perm(x, listw, nsim = nsim, alternative = "two.sided", 
        na.action = na.exclude)

# classify into LISA categories
    z_x <- base::scale(x)[, 1]
    lag_z <- spdep::lag.listw(listw, z_x)
    p_col <- grep("Pr", colnames(lmoran), value = TRUE)[1]
    p_val <- lmoran[, p_col]
    sig <- p_val < p_threshold
    lisa_class <- rep("Not significant", length(x))
    lisa_class[sig & z_x > 0 & lag_z > 0] <- "High-High"
    lisa_class[sig & z_x < 0 & lag_z < 0] <- "Low-Low"
    lisa_class[sig & z_x > 0 & lag_z < 0] <- "High-Low"
    lisa_class[sig & z_x < 0 & lag_z > 0] <- "Low-High"
    col_name <- paste0("lisa_", feature)
    sfe[[col_name]] <- factor(lisa_class, levels = c("High-High", 
        "Low-Low", "High-Low", "Low-High", "Not significant"))
    li_col_name <- paste0("local_moran_", feature)
    sfe[[li_col_name]] <- lmoran[, "Ii"]
    tab <- table(lisa_class)
    message("LISA clusters for '", feature, "':")
    for (cl in names(tab)) {
        message("  ", cl, ": ", tab[cl])
    }
    sfe
}
