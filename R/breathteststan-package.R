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
#' @importFrom instantiate stan_package_model
#' @importFrom stats rnorm rlnorm
#' @importFrom utils capture.output
#' @importFrom stringr str_extract str_match
#' @importFrom purrr map_df
#' @importFrom stats na.omit quantile
#' @importFrom rappdirs app_dir
#' @importFrom posterior as_draws_df
#' @importFrom tidyr separate_wider_regex pivot_wider
# usethis namespace: end
NULL

#' @export
.onLoad = function(libname, pkgname) {
  # https://wlandau.github.io/instantiate/
  # https://blog.r-hub.io/2020/03/12/user-preferences/#not-so-temporary-files3:
  if (FALSE) {
    cache_dir = rappdirs::app_dir("breathteststan", "dmenne")$cache()
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    Sys.setenv("CMDSTAN" = cache_dir)
    Sys.setenv("CMDSTAN_INSTALL" = "fixed")
  }
}
