<div id="main" class="col-md-9" role="main">

# Run the full EFA pipeline on an SFE object

<div class="ref-description section level2">

Run the full EFA pipeline on an SFE object

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
run_efa_pipeline(
  sfe,
  geom_name = "cellSeg",
  n_points = 200,
  n_harmonics = 12,
  normalize = TRUE,
  do_pca = TRUE,
  n_pcs = 10,
  compute_gof = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   sfe:

    SpatialFeatureExperiment object with cell segmentation

-   geom\_name:

    Geometry name (default "cellSeg")

-   n\_points:

    Points to resample each boundary to

-   n\_harmonics:

    Number of Fourier harmonics

-   normalize:

    Normalize coefficients

-   do\_pca:

    Run PCA on the shape space

-   n\_pcs:

    Number of PCs to retain

</div>

<div class="section level2">

## Value

Modified SFE object with EFA and optionally EFA\_PCA in reducedDims, and
efa\_complexity / efa\_ellipticity in colData

</div>

</div>
