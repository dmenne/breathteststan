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
#'   fit parameters \code{beta, k, m}.}
#'   \item{groups}{Which pairwise difference, e.g \code{solid - liquid}}
#'   \item{estimate}{Point estimate (chain mean) of the difference}
#'   \item{cred.low, cred.high}{Lower and upper 95 percent credible interval of difference.}
#' }
#' The chains of pairwise differences  are returned as a attribute \code{chain}
#' for use in plotting. See example below how to use these to display difference histograms.
#' @examples
#' \donttest{
#' library(dplyr)
#' library(breathtestcore)
#' data("usz_13c", package = "breathtestcore")
#' data = usz_13c %>%
#'   dplyr::filter( patient_id %in%
#'     c("norm_001", "norm_002", "norm_003", "norm_004", "pat_001", "pat_002","pat_003")) %>%
#'   cleanup_data()
#' fit = stan_group_fit(data, iter = 300, chains = 1) # Use more iterations!
#' cf = coef_diff_by_group(fit)
#' cc = attr(cf, "chain") %>%
#'    filter(key == "t50_maes_ghoos", abs(value) < 200) %>%
#'    mutate(
#'      groups = paste(group2, group1, sep = " - ")
#'    )
#' str(cc)
#' if (require(ggplot2)) {
#'   ggplot(cc, aes(x = value)) + geom_histogram() + facet_wrap(~groups)
#' }
#' # For comparison
#' fit = nlme_fit(data)
#' coef_diff_by_group(fit)
#' }
#' @importFrom stats confint relevel
#' @importFrom utils combn
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
  # Keep CRAN quiet
  . = cred.high = cred.low = estimate = group = group_i = key = method = parameter =
    value = value1 = value2 = NULL
  cm = comment(fit)
  # Keep CRAN quite
  # No differences if there is only one group
  nlevels = length(unique(fit$data$group))
  if (nlevels <=1){
    return(NULL)
  }
  sig = as.integer(options("digits"))

  # Local extractor function
  ex = function(par, i = NA) {
    if (is.na(i)) {
      as.vector(rstan::extract(fit$stan_fit, permuted = TRUE, pars = par)[[par]])
    } else
    {
      rstan::extract(fit$stan_fit, permuted = TRUE, pars = par)[[par]][,i]
    }
  }

  group_dict =
    (unique(fit$data %>% select(group, group_i)) %>%
       arrange(group_i))

  group_chain = group_dict %>%
    rowwise() %>%
    do (
      {
        data_frame(
          group = .$group,
          m = ex("mu_m") + ex("m_group",.$group_i),
          k = ex("mu_k") + ex("k_group",.$group_i),
          beta = ex("mu_beta") + ex("beta_group",.$group_i)
        )
      }
    )   %>%
    ungroup() %>%
    mutate(
      t50_maes_ghoos = t50_maes_ghoos(.),
      tlag_maes_ghoos = tlag_maes_ghoos(.),
      t50_maes_ghoos_scintigraphy = t50_maes_ghoos_scintigraphy(.),
      t50_bluck_coward = t50_bluck_coward(.),
      tlag_bluck_coward = tlag_bluck_coward(.)
    ) %>%
    tidyr::gather(key, value, -group) %>%
    ungroup()

  diffs = t(combn(group_dict$group,2))
  cf = data.frame()
  cf_chain = data.frame()
  for (i in 1:nrow(diffs)){
    group1 = diffs[i,1]
    group2 = diffs[i,2]
    dd = do_coef_by_group(group_chain, group1, group2)
    cf = rbind(cf, dd$d_summary)
    cf_chain = rbind(cf_chain, dd$d_chain)
  }
  attr(cf, "chain") = cf_chain
  comment(cf) = cm
  cf
  }



## Local function
do_coef_by_group = function(group_chain, group1, group2){
  # Keep CRAN cool
  cred.high = cred.low = estimate = group = key = method = parameter =
    value = value1 = value2 = NULL
  d_chain =
    cbind(group_chain %>% filter(group == group1) %>% select(key, value1 = value),
          group_chain %>% filter(group == group2) %>% select(value2 = value)) %>%
    na.omit() %>%
    mutate(
      group1 = group1,
      group2 = group2,
      value = value2-value1
    )

  d_summary = d_chain %>%
    group_by(group1, group2, key) %>%
    summarize(
      estimate = mean(value, na.rm = TRUE),
      cred.low = quantile(value, 0.0275, na.rm= TRUE),
      cred.high = quantile(value, 0.975, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      parameter = stringr::str_match(key, "k|m|beta|t50|tlag")[,1],
      method = stringr::str_match(key, "maes_ghoos_scintigraphy|maes_ghoos|bluck_coward|exp_beta")[,1],
      method = ifelse(is.na(method), "exp_beta", method),
      groups = paste(group2, group1,sep = " - ")
    ) %>%
    select(parameter, method, groups,  estimate, cred.low, cred.high)

  list(d_summary = d_summary, d_chain = d_chain)
}
