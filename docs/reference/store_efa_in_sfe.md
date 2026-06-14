<div id="main" class="col-md-9" role="main">

# Add EFA coefficients as a reducedDim in the SFE object

<div class="ref-description section level2">

Add EFA coefficients as a reducedDim in the SFE object

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
store_efa_in_sfe(sfe, efa_result, name = "EFA")
```

</div>

</div>

<div class="section level2">

## Arguments

-   sfe:

    A SpatialFeatureExperiment object

-   efa\_result:

    Output of compute\_efa()

-   name:

    Name for the reducedDim slot (default "EFA")

</div>

<div class="section level2">

## Value

The modified SFE object

</div>

</div>
