library(testit)

assert("Data that cannot be fitted with nls_list/nlme work with stan_fit", {
  # with this seed, cf[10] does not fit with nls_list
  library(breathtestcore) # Debug
  library(breathteststan)

  chains = 1
  student_t_df = 10
  dose = 100
  iter = 400
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(seed = 100)$data)
  comment(data) = "comment"
  fit = stan_fit(
    data,
    dose = dose,
    student_t_df = student_t_df,
    chains = chains,
    iter = iter
  )
  # https://yihui.org/en/2026/05/testthat-to-testit/
  (inherits(fit, "breathtestfit"))
  (inherits(fit, "breathteststanfit"))
  (inherits(fit$stan_fit, "CmdStanMCMC"))
  (inherits(fit$stan_fit, "CmdStanFit"))
  (identical(names(fit), c("coef", "data", "stan_fit", "coef_chain")))
  (all.equal(names(fit$data), names(data)))
  (sigma(fit) > 0.8)
  (identical(comment(fit), "comment"))

  cf = fit$coef
  (identical(unique(cf$group), "A"))
  (identical(unique(cf$parameter), c("beta", "k", "m", "t50", "tlag")))
  (identical(
    unique(cf$stat),
    c("estimate", "q_0275", "q_25", "q_75", "q_975")
  ))
  (nrow(cf) > 395)
  (identical(ncol(cf), 6L))

  cf = coef(fit) # This is the "mean" group only
  (identical(unique(cf$group), "A"))
  (identical(unique(cf$parameter), c("beta", "k", "m", "t50", "tlag")))
  (nrow(cf) > 79)
  (identical(ncol(cf), 5L))
})
