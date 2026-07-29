test_that("Data that cannot be fitted with nls_list/nlme work with stan_fit", {
  # with this seed, cf[10] does not fit with nls_list
  library(breathtestcore)

  chains = 1
  student_t_df = 10
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(seed = 100)$data)
  comment(data) = "comment"
  capture.output({
    fit = stan_fit(
      data,
      dose = dose,
      student_t_df = student_t_df,
      chains = chains,
      iter = iter
    )
  })
  expect_s3_class(fit, "breathtestfit")
  expect_s3_class(fit, "breathteststanfit")
  expect_r6_class(fit$stan_fit, "CmdStanMCMC")
  expect_r6_class(fit$stan_fit, "CmdStanFit")
  expect_identical(names(fit), c("coef", "data", "stan_fit", "coef_chain"))
  expect_equal(names(fit$data), names(data))
  expect_gte(sigma(fit), 0.8)
  expect_identical(comment(fit), "comment")

  cf = fit$coef
  expect_identical(unique(cf$group), "A")
  expect_identical(unique(cf$parameter), c("beta", "k", "m", "t50", "tlag"))
  expect_identical(
    unique(cf$stat),
    c("estimate", "q_0275", "q_25", "q_75", "q_975")
  )
  expect_gte(nrow(cf), 395)
  expect_equal(ncol(cf), 6)

  cf = coef(fit) # This is the "mean" group only
  expect_identical(unique(cf$group), "A")
  expect_identical(unique(cf$parameter), c("beta", "k", "m", "t50", "tlag"))
  expect_gte(nrow(cf), 79)
  expect_equal(ncol(cf), 5)
})
