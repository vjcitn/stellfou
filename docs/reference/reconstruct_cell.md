<div id="main" class="col-md-9" role="main">

# Reconstruct a cell boundary from its EFA coefficients

<div class="ref-description section level2">

Reconstruct a cell boundary from its EFA coefficients

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
reconstruct_cell(efa_result, cell_id, n_harmonics_use = NULL, n_points = 300)
```

</div>

</div>

<div class="section level2">

## Arguments

-   efa\_result:

    Output of compute\_efa()

-   cell\_id:

    Character; which cell to reconstruct

-   n\_harmonics\_use:

    How many harmonics to use in reconstruction (NULL = all available)

-   n\_points:

    Number of points in the reconstructed contour

</div>

<div class="section level2">

## Value

(n\_points x 2) matrix of reconstructed coordinates

</div>

</div>
