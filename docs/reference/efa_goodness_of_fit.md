<div id="main" class="col-md-9" role="main">

# Compute per-cell EFA goodness-of-fit metrics

<div class="ref-description section level2">

For each cell, reconstructs the boundary from the stored EFA
coefficients and compares to the original boundary.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
efa_goodness_of_fit(boundaries, efa_result, n_harmonics_use = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   boundaries:

    Named list of (n x 2) coordinate matrices

-   efa\_result:

    Output of compute\_efa()

-   n\_harmonics\_use:

    Number of harmonics for reconstruction (NULL = all available)

</div>

<div class="section level2">

## Value

A data.frame with one row per cell and columns: - mean\_dev: mean
distance from original to nearest reconstructed point - max\_dev: max
such distance (one-sided Hausdorff) - symmetric\_hausdorff: max of both
one-sided Hausdorff distances - rel\_mean\_dev: mean\_dev / perimeter
(dimensionless) - area\_ratio: area(reconstruction) / area(original) -
power\_captured: fraction of total harmonic power in first
n\_harmonics\_use harmonics (1.0 if using all)

</div>

</div>
