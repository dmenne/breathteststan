#' @title Tabulates breath test parameter differences of groups
#' @description Given a fit to 13C breath test curves, computes absolute values and
#' their conficence intervals of parameters from the fit,
#' for examples of the half emptying time \code{t50}.
#'
#' @param fit Object of class \code{breathtestfit}, for example from
#' \code{\link{nlme_fit}}, \code{\link{nls_fit}} or \code{\link[breathteststan]{stan_fit}}
#' @param mcp_group "Tukey" (default) for all pairwise comparisons, "Dunnett" for
#' comparisons relative to the reference group.
#' @param reference_group Used as the first group and as reference group for
#' \code{mcp_group == "Dunnett"}
#' @param ... Not used
#'
#' @return A \code{tibble} with columns
#' \describe{
#'   \item{parameter}{Parameter of fit, e.g. \code{beta, k, m, t50}}
#'   \item{method}{Method used to compute parameter. \code{exp_beta} refers to primary
#'   fit parameters \code{beta, k, m}.
#'   \item{groups}{Which pairwise difference, e.g \code{solid - liquid}}
#'   \item{estimate}{Estimate of the difference}
#'   \item{cred.low, cred.high}{Lower and upper 95% confidence interval of difference.
#'   A comparison is significantly different from zero when both estimates have the same
#'   sign.}
#' }
#' @examples
#' library(dplyr)
#' data("usz_13c")
#' data = usz_13c %>%
#'   dplyr::filter( patient_id %in%
#'     c("norm_001", "norm_002", "norm_003", "norm_004", "pat_001", "pat_002","pat_003")) %>%
#'   cleanup_data()
#' fit = nls_group_fit(data)
#' coef_diff_by_group(fit)
#' \dontrun{
#' fit = nlme_fit(data)
#' coef_by_group(fit)
#' }
#' @importFrom stats confint relevel
#' @import multcomp
#' @export
coef_diff_by_group = function(fit, mcp_group = "Tukey", reference_group = NULL, ...) {
  UseMethod("coef_diff_by_group", fit)
}

#' @export
coef_diff_by_group.breathtestfit =
  function(fit, mcp_group = "Tukey", reference_group = NULL, ...) {
  if (!inherits(fit, "breathtestfit")) {
    stop("Function coef_diff_by_group: parameter 'fit' must inherit from class breathtestfit")
  }
  cf = coef(fit)
  if (is.null(cf))
    return(NULL)
  if (! mcp_group %in% c("Dunnett", "Tukey")){
    stop("Function coeff_diff_by_group: mcp_group must be 'Dunnett' or 'Tukey', but is ",
         mcp_group)
  }
  cm = comment(cf)
  # Keep CRAN quite
  . = confint = estimate.x = estimate.y = lhs = method = parameter = rhs = statistic =
    std.error = conf.low = conf.high = p.value = NULL
  cf = cf %>%
    mutate( # lme requires factors
      group = as.factor(.$group)
    )
  # No differences if there is only one group
  if (nlevels(cf$group) <=1){
    return(NULL)
  }
  if (!is.null(reference_group) && !(reference_group %in% levels(cf$group))) {
    stop("Function coeff_diff_by_group: reference_group must be a level in coef(fit)$group")
  }
  if (!is.null(reference_group)){
    cf$group = relevel(cf$group, reference_group)
  }
  sig = as.integer(options("digits"))
  cf = cf %>%
    group_by(parameter, method) %>%
    do({
      fit_lme = nlme::lme(value~group, random = ~1|patient_id, data = .)
      glh = multcomp::glht(fit_lme, linfct =  multcomp::mcp(group = mcp_group))
      broom::tidy(confint(glh))[,-2] %>%
        left_join(broom::tidy(summary(glh)), by = "lhs", copy = TRUE) %>%
      mutate(
        estimate.x = signif(estimate.x, sig),
        conf.low = signif(conf.low, sig),
        conf.high = signif(conf.high, sig),
        p.value = signif(p.value, sig)
      )
    }) %>%
    ungroup() %>%
    dplyr::select(-rhs, -estimate.y, -std.error, -statistic,
      estimate = estimate.x, groups = lhs)
  comment(cf) = cm
  cf
}
