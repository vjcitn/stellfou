<div id="main" class="col-md-9" role="main">

# sfe\_efa\_bridge.R – code created in dialogue with claude.ai Opus 4.6 medium, June 3 2026

<div class="ref-description section level2">

Bridge between SpatialFeatureExperiment cell boundary polygons and
Elliptic Fourier Analysis (Momocs).

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
extract_cell_boundaries(sfe, geom_name = "cellSeg", n_points = 200)
```

</div>

</div>

<div class="section level2">

## Arguments

-   sfe:

    A SpatialFeatureExperiment object

-   geom\_name:

    Name of the colGeometry, typically "cellSeg"

-   n\_points:

    If not NULL, resample each boundary to this many equally-spaced
    points (recommended for EFA comparability)

</div>

<div class="section level2">

## Value

A named list of (n x 2) matrices, one per cell

</div>

<div class="section level2">

## Details

Extracts cell segmentation polygons from an SFE object, computes EFA
coefficients per cell, stores them as a reducedDim or colData, and
provides reconstruction/visualization.

Dependencies: BiocManager::install("SpatialFeatureExperiment")
install.packages("Momocs") \# optional: install.packages("pliman")
Extract cell boundary coordinate matrices from an SFE object

</div>

</div>
