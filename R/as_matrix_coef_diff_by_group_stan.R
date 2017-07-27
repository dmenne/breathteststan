#' @title S3 as.matrix for result of coef_diff_by_group
#' @description Generates a matrix that can be used with plotting functions from
#' package \code{\link[bayesplot]{mcmc_hist}}.
#'
#' @param x Result of a call to \code{coef_diff_by_group(fit)}
#' @param ... parameter name as string, e.g. \code{"m", "k", "beta", "t50_bluck_coward"}.
#' When missing, \code{"t50_maes_ghoos"} is assumed.
#'
#' @return mcmc array with columns of differences for use with functions from packages
#' bayesplot or coda
#' @examples
#' # \donttest{
#' library(dplyr)
#' library(breathtestcore)
#' library(ggplot2)
#' data("usz_13c", package = "breathtestcore")
#' data = usz_13c %>%
#'   dplyr::filter( patient_id %in%  c("norm_001", "norm_002", "norm_003",
#'                         "norm_004", "pat_001", "pat_002","pat_003")) %>%
#'   cleanup_data()
#' fit = stan_group_fit(data, iter = 300, chains = 1)
#' cf = coef_diff_by_group(fit)
#' str(cf)
#' # Calling without parameters gives Maes/Ghoos t50
#' bayesplot::mcmc_hist(as.matrix(cf))
#' # Use a function from the bayesplot universe
#' dens = bayesplot::mcmc_dens(as.matrix(cf, parameter = "m"))
#' # use suppressMessages to avoid a message "another scale"
#' suppressMessages(
#'   dens + geom_vline(xintercept = 0) + scale_x_continuous(limits= c(-20,10)))
#' # }
#' @importFrom tidyr spread
#' @export
as.matrix.coef_diff_by_group_stan = function(x, ...){
  parameter = list(...)$parameter
  if (is.null(parameter))
    parameter = "t50_maes_ghoos"
  chain = attr(x, "chain")
  if (is.null(chain))
    stop("Attribute chain is missing")
  p_meth = unique(chain$key)
  diff_groups = group1 = group2 = key = NULL # mute CRAN
  if (! parameter %in% p_meth)
    stop("There is no parameter ", parameter, " in chain")
  ret = chain %>%
    filter(key == parameter) %>%
    mutate(
      diff_groups = paste(group2, group1, sep = " - ")
    ) %>%
    group_by(diff_groups) %>%
    mutate(
      n = row_number()
    ) %>%
    select(n, diff_groups, diff)  %>%
    tidyr::spread(diff_groups, diff) %>%
    ungroup() %>%
    select(-n) %>%
    as.matrix()
}

