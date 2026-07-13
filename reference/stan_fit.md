# Bayesian Stan fit to 13C Breath Data

Fits exponential beta curves to 13C breath test series data using
Bayesian Stan methods. See
<https://menne-biomed.de/blog/breath-test-stan/> for a comparison
between single curve, mixed-model population and Bayesian methods.

## Usage

``` r
stan_fit(
  data,
  dose = 100,
  sample_minutes = 15,
  student_t_df = 10,
  chains = 2,
  iter = 1000,
  model = "breath_test_1",
  seed = 4711
)
```

## Arguments

- data:

  Data frame or tibble as created by
  [`cleanup_data`](https://rdrr.io/pkg/breathtestcore/man/cleanup_data.html),
  with mandatory columns `patient_id, group, minute` and `pdr`. It is
  recommended to run all data through
  [`cleanup_data`](https://rdrr.io/pkg/breathtestcore/man/cleanup_data.html)
  which will insert dummy columns for `patient_id` and `minute` if the
  data are distinct, and report an error if not. Since the Bayesian
  method is stabilized by priors, it is possible to fit single curves.

- dose:

  Dose of acetate or octanoate. Currently, only one common dose for all
  records is supported.

- sample_minutes:

  If mean sampling interval is \< sampleMinutes, data are subsampled
  using a spline algorithm

- student_t_df:

  When student_t_df \< 10, the student distribution is used to model the
  residuals. Recommended values to model typical outliers are from 3
  to 6. When student_t_df \>= 10, the normal distribution is used.

- chains:

  Number of chains for Stan

- iter:

  Number of iterations for each Stan chain

- model:

  Name of model; use `names(stanmodels)` for other models.

- seed:

  Optional seed for rstan

## Value

A list of classes "breathteststanfit" and "breathtestfit" with elements

- `coef` Estimated parameters as data frame in a key-value format with
  columns `patient_id, group, parameter, method` and `value`. Has an
  attribute AIC.

- `data` The effectively analyzed data. If density of points is too
  high, e.g. with BreathId devices, data are subsampled before fitting.

- `stan_fit` The Stan fit for use with `shinystan::launch_shiny` or
  extraction of chains.

## See also

Base methods `coef, plot, print`; methods from package
`broom: tidy, augment`.

## Examples

``` r
library(breathtestcore)
suppressPackageStartupMessages(library(dplyr))
d = breathtestcore::simulate_breathtest_data(n_records = 3) # default 3 records
data = breathtestcore::cleanup_data(d$data)
# Use more than 80 iterations and 4 chains for serious fits
fit = stan_fit(data, chains = 1, iter = 80)
plot(fit) # calls plot.breathtestfit

# Extract coefficients and compare these with those
# used to generate the data
options(digits = 2)
cf = coef(fit)
cf %>%
  filter(grepl("m|k|beta", parameter )) %>%
  select(-method, -group) %>%
  tidyr::spread(parameter, value) %>%
  inner_join(d$record, by = "patient_id") %>%
  select(patient_id, m_in = m.y, m_out = m.x,
         beta_in = beta.y, beta_out = beta.x,
         k_in = k.y, k_out = k.x)
#> # A tibble: 3 × 7
#>   patient_id  m_in m_out beta_in beta_out    k_in  k_out
#>   <chr>      <dbl> <dbl>   <dbl>    <dbl>   <dbl>  <dbl>
#> 1 rec_01        26  23.9    1.27     1.35 0.00999 0.0114
#> 2 rec_02        43  42.2    1.90     1.93 0.0116  0.0119
#> 3 rec_03        16  16.3    1.90     1.84 0.0129  0.0124
# For a detailed analysis of the fit, use the shinystan library
# \donttest{
library(shinystan)
#> Loading required package: shiny
#> 
#> This is shinystan version 2.7.0
# launch_shinystan(fit$stan_fit)
# }
# The following plots are somewhat degenerate because
# of the few iterations in stan_fit
suppressPackageStartupMessages(library(rstan))
stan_plot(fit$stan_fit, pars = c("beta[1]","beta[2]","beta[3]"))
#> ci_level: 0.8 (80% intervals)
#> outer_level: 0.95 (95% intervals)

stan_plot(fit$stan_fit, pars = c("k[1]","k[2]","k[3]"))
#> ci_level: 0.8 (80% intervals)
#> outer_level: 0.95 (95% intervals)

stan_plot(fit$stan_fit, pars = c("m[1]","m[2]","m[3]"))
#> ci_level: 0.8 (80% intervals)
#> outer_level: 0.95 (95% intervals)

```
