<div id="main" class="col-md-9" role="main">

# Compute EFA for all cell boundaries

<div class="ref-description section level2">

Compute EFA for all cell boundaries

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
compute_efa(boundaries, n_harmonics = 12, normalize = TRUE)
```

</div>

</div>

<div class="section level2">

## Arguments

-   boundaries:

    Named list of (n x 2) coordinate matrices (output of
    extract\_cell\_boundaries)

-   n\_harmonics:

    Number of Fourier harmonics

-   normalize:

    Logical; normalize coefficients for size, rotation, and
    starting-point invariance

</div>

<div class="section level2">

## Value

A list with components: - coefficients: (n\_cells x 4\*n\_harmonics)
matrix - raw\_efa: list of raw efourier output per cell - n\_harmonics:
the number of harmonics used

</div>

</div>
