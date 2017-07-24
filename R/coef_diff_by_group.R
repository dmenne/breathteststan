#' @title Tabulates breath test parameter differences of groups from Stan group fit
#' @description Given a Stan fit with grouping to 13C breath test curves,
#' computes point estimated and Bayesian credible intervals for all group pair
#' differences, for examples of the half emptying time \code{t50}.
#'
#' @param fit Object of class \code{breathtestfit}, for example from
#' \code{\link{nlme_fit}}, \code{\link{nls_fit}} or \code{\link[breathteststan]{stan_fit}}
#' @param mcp_group Not used, always all pairs are compared
#' @param reference_group Not used
#' @param ... Not used
#'
#' @return A \code{tibble} with columns
#' \describe{
#'   \item{parameter}{Parameter of fit, e.g. \code{beta, k, m, t50}}
#'   \item{method}{Method used to compute parameter. \code{exp_beta} refers to primary
#'   fit parameters \code{beta, k, m}.
#'   \item{groups}{Which pairwise difference, e.g \code{solid - liquid}}
#'   \item{estimate}{Point estimate (chain mean) of the difference}
#'   \item{cred.low, cred.high}{Lower and upper 95% credible interval of difference.}
#' }
#' @examples
#' library(dplyr)
#' data("usz_13c")
#' data = usz_13c %>%
#'   dplyr::filter( patient_id %in%
#'     c("norm_001", "norm_002", "norm_003", "norm_004", "pat_001", "pat_002","pat_003")) %>%
#'   cleanup_data()
#' fit = stan_group_fit(data)
#' coef_diff_by_group(fit)
#' \dontrun{
#' fit = nlme_fit(data)
#' coef_diff_by_group(fit)
#' }
#' @importFrom stats confint relevel
#' @export
coef_diff_by_group.breathteststangroupfit =
  function(fit, mcp_group = NULL, reference_group = NULL, ...) {
  if (!inherits(fit, "breathteststangroupfit")) {
    stop("Function coef_diff_by_group: parameter 'fit' must inherit from class breathteststangropfit")
  }
  if (! is.null(mcp_group) ){
    stop("Function coeff_diff_by_group: mcp_group must be NULL for Stan group fits")
  }
  if (! is.null(mcp_group) ){
      stop("Function coeff_diff_by_group: reference_group must be NULL for Stan group fits")
  }
  cm = comment(fit$data)
  # Keep CRAN quite
  # No differences if there is only one group
  nlevels = length(unique(fit$data$group))
  if (nlevels <=1){
    return(NULL)
  }
  sig = as.integer(options("digits"))
  ch = fit$coef_chain
  diffs = t(combn(unique(ch$group),2))

  #comment(cf) = cm
  NULL
}
