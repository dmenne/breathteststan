# S3 method to exctract the residual standard deviation

Functions for S3 method defined in breathtestcore for `stan_fit` and
`stan_group fit`.

## Usage

``` r
# S3 method for class 'breathteststanfit'
sigma(object, ...)
```

## Arguments

- object:

  A Stan-based fit

- ...:

  Not used

## Value

A numeric value giving the sigma (= average residual standard deviation)
term from the Stan fit.
