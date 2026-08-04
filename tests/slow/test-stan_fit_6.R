assert("Bad records are skipped when multiple in one file ", {
  library(breathtestcore)
  library(breathteststan)
  chains = 1
  student_t_df = 10
  dose = 100
  iter = 300
  # This file has two bad records
  filename = btcore_file("short_record.xml")
  xml_data = read_any_breathtest(filename)
  data = cleanup_data(xml_data)
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
