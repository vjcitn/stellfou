
# ============================================================
# Example usage (not run)
# ============================================================
if (TRUE) {
  library(SpatialFeatureExperiment)
  library(SFEData)

  # --- Load a Xenium or MERFISH dataset with cell segmentations ---
  # sfe <- readXenium("/path/to/xenium_output")
  # or use a demo dataset:
  # sfe <- XeniumOutput("v2", data_subset = "small")
  sfepath = XeniumOutput("v2")
  sfe = readXenium(sfepath)

  # --- Run the full pipeline ---
  sfe <- run_efa_pipeline(sfe, n_harmonics = 12)

  # --- Inspect shape descriptors and goodness-of-fit ---
  head(colData(sfe)[, c("efa_complexity", "efa_ellipticity")])

  # GoF summary across all cells
  gof_cols <- grep("^efa_gof_", names(colData(sfe)), value = TRUE)
  summary(as.data.frame(colData(sfe)[, gof_cols]))

  # Distribution of relative mean deviation
  hist(sfe$efa_gof_rel_mean_dev, breaks = 50,
       main = "EFA reconstruction error\n(mean deviation / perimeter)",
       xlab = "Relative mean deviation", col = "grey80")

  # Identify cells with poor fits (may need more harmonics)
  poor_fit <- which(sfe$efa_gof_rel_mean_dev > 0.01)
  message(length(poor_fit), " cells have > 1% relative mean deviation")

  # --- Visualize progressive reconstruction for a few cells ---
  boundaries <- extract_cell_boundaries(sfe)
  efa_result <- compute_efa(boundaries, n_harmonics = 12)
  plot_efa_reconstruction(boundaries, efa_result,
                          cell_ids = names(boundaries)[1:3],
                          n_harmonics_seq = c(1, 3, 6, 12))

  # --- PCA on shape space is already in reducedDim ---
  # Can be used with Voyager for spatial autocorrelation of shape:
  library(Voyager)
  colGraph(sfe, "knn") <- findSpatialNeighbors(sfe, method = "knearneigh", k = 6)
  # Moran's I on shape PC1:
  sfe$shape_pc1 <- reducedDim(sfe, "EFA_PCA")[, 1]
  sfe <- colDataMoransI(sfe, "shape_pc1")

  # --- Correlate shape with gene expression ---
  # e.g., test whether shape complexity associates with a gene:
  # Use log1p(counts) if logcounts assay is not yet available
  count_mat <- assay(sfe, "counts")
  # Pick a gene present in the data
  test_gene <- intersect(c("VIM", "KRT8", "EPCAM"), rownames(sfe))[1]
  if (!is.na(test_gene)) {
    plot(sfe$efa_ellipticity,
         log1p(count_mat[test_gene, ]),
         xlab = "Ellipticity (EFA)",
         ylab = paste0(test_gene, " log1p(counts)"),
         pch = 16, cex = 0.3, col = adjustcolor("steelblue", 0.4))
  }
}


