<div id="main" class="col-md-9" role="main">

# Plot original vs reconstructed cell boundaries

<div class="ref-description section level2">

Plot original vs reconstructed cell boundaries

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
plot_efa_reconstruction(
  boundaries,
  efa_result,
  cell_ids = NULL,
  n_harmonics_seq = c(1, 3, 6, 12),
  ncol = length(n_harmonics_seq)
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   boundaries:

    Named list of coordinate matrices

-   efa\_result:

    Output of compute\_efa()

-   cell\_ids:

    Character vector of cell IDs to plot

-   n\_harmonics\_seq:

    Harmonic counts to show progressive reconstruction

-   ncol:

    Number of columns in the plot layout

</div>

</div>
