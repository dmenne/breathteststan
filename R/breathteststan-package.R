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
#' @importFrom rlang .data
#' @import cmdstanr
#' @importFrom stats rnorm rlnorm
#' @importFrom utils capture.output
#' @importFrom stringr str_extract str_match
#' @importFrom purrr map_df
#' @importFrom stats na.omit quantile
#' @importFrom here here
#' @importFrom posterior as_draws_df
#' @importFrom tidyr separate_wider_regex pivot_wider
#' @importFrom breathtestcore t50_maes_ghoos tlag_maes_ghoos t50_maes_ghoos_scintigraphy
#' @importFrom breathtestcore t50_bluck_coward tlag_bluck_coward subsample_data
#' @importFrom breathtestcore cleanup_data coef_by_group
# usethis namespace: end
"_PACKAGE"
