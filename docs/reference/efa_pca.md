<div id="main" class="col-md-9" role="main">

# PCA on EFA coefficients and optionally store in SFE

<div class="ref-description section level2">

PCA on EFA coefficients and optionally store in SFE

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
efa_pca(sfe, efa_dim_name = "EFA", n_pcs = 10, store_name = "EFA_PCA")
```

</div>

</div>

<div class="section level2">

## Arguments

-   sfe:

    SFE object with EFA reducedDim already stored

-   efa\_dim\_name:

    Name of the EFA reducedDim

-   n\_pcs:

    Number of PCs to retain

-   store\_name:

    Name for the PCA reducedDim (NULL = don't store)

</div>

<div class="section level2">

## Value

prcomp result, invisibly

</div>

</div>
