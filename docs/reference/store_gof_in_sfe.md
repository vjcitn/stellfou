<div id="main" class="col-md-9" role="main">

# Store GoF metrics in colData of an SFE object

<div class="ref-description section level2">

Store GoF metrics in colData of an SFE object

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
store_gof_in_sfe(sfe, gof_df, prefix = "efa_gof_")
```

</div>

</div>

<div class="section level2">

## Arguments

-   sfe:

    SpatialFeatureExperiment object

-   gof\_df:

    Output of efa\_goodness\_of\_fit()

-   prefix:

    Column name prefix (default "efa\_gof\_")

</div>

<div class="section level2">

## Value

Modified SFE object

</div>

</div>
