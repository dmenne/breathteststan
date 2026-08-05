# S3 method to extract the residual standard deviation

Functions for S3 method defined in breathtestcore for `CmdStanFit`.

## Usage

``` r
# S3 method for class 'breathteststanfit'
sigma(object, ...)
```

## Arguments

- object:

  A Stan-based fit returning an S6 method

- ...:

  Not used

## Value

A numeric value giving the sigma (= average residual standard deviation)
term from the Stan fit.
