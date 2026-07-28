#' @title breathteststan
#' @description Stan fits to 13C breath test curves. This
#'   is an add-on package to package breathtestcore (https://github.com/dmenne/breathtestcore)
#'
#' @name breathteststan-package
#' @aliases breathteststan
# usethis namespace: start
#' @import Rcpp
#' @import dplyr
#' @import methods
#' @import rstantools
#' @import cmdstanr
#' @useDynLib breathteststan, .registration = TRUE
#' @importFrom instantiate stan_package_model
#' @importFrom stats rnorm rlnorm
#' @importFrom utils capture.output
#' @importFrom stringr str_extract str_match
#' @importFrom purrr map_df
#' @importFrom stats na.omit quantile
#' @importFrom tidyr separate_wider_regex pivot_wider
# usethis namespace: end
NULL
