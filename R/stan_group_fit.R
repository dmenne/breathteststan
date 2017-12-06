#' @title Bayesian Stan fit to 13C Breath Data in Multiple Groups
#' @description Fits exponential beta curves to 13C breath test series
#' data using Bayesian Stan methods, by assuming fixed between group effects.
#' This model is overly parsiomonious. Do not use it
#' unless you check the results carefully and understand why
#' fits can be very bad.
#'
#' @param data Data frame or tibble as created by \code{\link{cleanup_data}},
#' with mandatory columns \code{patient_id, group, minute} and \code{pdr}.
#' It is recommended to run all data through \code{\link{cleanup_data}} which
#' will insert dummy columns for \code{patient_id} and \code{minute} if the
#' data are distinct, and report an error if not. Since the Bayesian method
#' is stabilized by priors, it is possible to fit single curves.
#' @param dose Dose of acetate or octanoate. Currently, only one common dose
#' for all records is supported.
#' @param sample_minutes If mean sampling interval is < sampleMinutes, data are subsampled
#' using a spline algorithm
#' @param student_t_df When student_t_df < 10, the student distribution is used to
#' model the residuals. Recommended values to model typical outliers are from 3 to 6.
#' When student_t_df >= 10, the normal distribution is used.
#' @param chains Number of chains for Stan
#' @param iter Number of iterations for each Stan chain
#' @param model Name of model; use \code{names(stanmodels)} for other models.
#'
#'
#' @return A list of classes "breathteststangroupfit", "breathteststanfit" and "breathtestfit"
#' with elements
#' \itemize{
#'   \item {\code{coef} Estimated parameters as data frame in a key-value format with
#'    columns \code{patient_id, group, parameter, method} and \code{value}.
#'    Has an attribute AIC.}
#'    \item {\code{data}  The effectively analyzed data. If density of points
#'    is too high, e.g. with BreathId devices, data are subsampled before fitting.}
#'    \item {\code{stan_fit} The Stan fit for use with \code{shinystan::launch_shiny}
#'    or extraction of chains. }
#' }
#' @seealso Base methods \code{coef, plot, print}; methods from package
#'  \code{broom: tidy, augment}.
#' @examples
#' \donttest{
#' library(breathtestcore)
#' library(dplyr)
#' data("usz_13c", package = "breathtestcore")
#' data = usz_13c %>%
#'   dplyr::filter( patient_id %in%
#'        c("norm_001", "norm_002", "norm_003", "norm_004",
#'          "pat_001", "pat_002","pat_003")) %>%
#'   cleanup_data()
#' fit = stan_group_fit(data, chains = 1, iter = 100)
#' plot(fit) # calls plot.breathtestfit
#' coef(fit)
#' }
#'
#' @export
#'
stan_group_fit = function(data, dose = 100, sample_minutes = 15, student_t_df = 10,
                    chains = 2, iter = 1000, model = "breath_test_group_1") {

  if (length(unique(data$group)) < 2)
    stop("Use stan_fit if there is only one group")
  # Avoid notes on CRAN
  value = patient_id = group = minute = pdr = NULL
  stat = estimate = . = k = key =  m = q_975 = NULL
  cm  = comment(data)
  data = breathtestcore::subsample_data(data, sample_minutes) %>%
    mutate(
      pat_i =  as.integer(as.factor(patient_id)),
      group_i =  as.integer(as.factor(group))
    )
  n_pat = max(data$pat_i)
  n_group = max(data$group_i)
  data_list = list(
    n = nrow(data),
    n_pat = n_pat,
    n_group = n_group,
    student_t_df = 5,
    dose = 100,
    pat_i = data$pat_i,
    group_i = data$group_i,
    minute = data$minute,
    pdr = data$pdr)

  # Note: as.array is required to handle the case of n_pat = 1
  init = rep(list(list(
    m_pat_raw = as.array(rnorm(n_pat,0,2)),
    m_group = rnorm(n_group,0,.1),
    sigma_m_pat = abs(rnorm(1, 10,1)),
    mu_m = rnorm(1, 40, 2),

    k_pat_raw = as.array(rnorm(n_pat,0,.0001)),
    k_group = rnorm(n_group,0,.0001),
    sigma_k_pat = abs(rnorm(1, 0,.0001)),
    mu_k = rnorm(1, 40, 3),

    beta_pat_raw = as.array(rnorm(n_pat,0,.1)),
    beta_group = rnorm(n_group,0,.1),
    sigma_beta_pat = abs(rnorm(1, 0, .1)),
    mu_beta = rnorm(1, 2,0.1),

    sigma = abs(rnorm(1, 10, 1))

  )),chains)


  if (!exists("stanmodels"))
    stop("stanmodels not found")
  mod = stanmodels[[model]]
  if (is.null(mod))
    stop("stanmodels ", model, " not found")
  options(mc.cores = min(chains, max(parallel::detectCores()/2, 1)))
  capture.output({fit = suppressWarnings(
    rstan::sampling(mod, data = data_list, init = init,
                    control = list(adapt_delta = 0.9),
                    iter =  iter, chains = chains)
  )})


  # Local extractor function
  ex = function(par, i = NA) {
    if (is.na(i)) {
      as.vector(rstan::extract(fit, permuted = TRUE, pars = par)[[par]])
    } else
    {
      rstan::extract(fit, permuted = TRUE, pars = par)[[par]][,i]
    }
  }

  coef_chain = data %>%
    select(-minute, -pdr) %>%
    distinct() %>%
    rowwise() %>%
    do (
      {
        data_frame(
          patient_id = .$patient_id,
          group = .$group,
          m = ex("mu_m") + ex("m_group",.$group_i) + ex("m_pat", .$pat_i),
          k = ex("mu_k") + ex("k_group",.$group_i) + ex("k_pat", .$pat_i),
          beta = ex("mu_beta") + ex("beta_group",.$group_i) + ex("beta_pat", .$pat_i)
        )
      }
    ) %>%
    ungroup() %>%
    mutate(
      t50_maes_ghoos = breathtestcore::t50_maes_ghoos(.),
      tlag_maes_ghoos = breathtestcore::tlag_maes_ghoos(.),
      t50_maes_ghoos_scintigraphy = breathtestcore::t50_maes_ghoos_scintigraphy(.),
      t50_bluck_coward = breathtestcore::t50_bluck_coward(.),
      tlag_bluck_coward = breathtestcore::tlag_bluck_coward(.)
    ) %>%
    rename(m_exp_beta = m, k_exp_beta = k, beta_exp_beta = beta) %>%
    tidyr::gather(key, value, -patient_id, -group) %>%
    na.omit() %>%
    ungroup()

  cf = coef_chain %>%
    group_by(patient_id, group, key) %>%
    summarize(
      estimate = mean(value),
      q_0275 = quantile(value, 0.0275),
      q_25 = quantile(value, 0.25),
      q_75 = quantile(value, 0.75),
      q_975 = quantile(value, 0.975)
    ) %>%
    ungroup() %>%
    mutate(
      parameter = stringr::str_match(key, "k|m|beta|t50|tlag")[,1],
      method = stringr::str_match(key, "maes_ghoos_scintigraphy|maes_ghoos|bluck_coward|exp_beta")[,1]
    ) %>%
    select(-key) %>%
    tidyr::gather(stat, value, estimate:q_975)

  ret = list(coef = cf, data = data, stan_fit = fit, coef_chain = coef_chain)
  class(ret) = c("breathteststangroupfit", "breathteststanfit", "breathtestfit")
  comment(ret) = cm # Recover comment
  ret
}

