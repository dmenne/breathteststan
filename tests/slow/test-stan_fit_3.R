library(testit)

assert("A single record can be fitted", {
  library(breathtestcore)
  library(breathteststan)
  chains = 1
  student_t_df = 10
  dose = 100
  iter = 300
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(n_records = 1, seed = 100)$data)
  fit = stan_fit(
    data,
    dose = dose,
    student_t_df = student_t_df,
    chains = chains,
    iter = iter
  )
  (inherits(fit, "breathtestfit"))
  (inherits(fit, "breathteststanfit"))
  (inherits(fit$stan_fit, "CmdStanMCMC"))
  (inherits(fit$stan_fit, "CmdStanFit"))
})
