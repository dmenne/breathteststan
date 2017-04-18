context("Bayesian fit with Stan")

test_that("stanmodels exist", {
  expect_is(breathteststan:::stanmodels,"list")
  expect_is(breathteststan:::stanmodels$breath_test_1, "stanmodel")
  expect_gte(length(breathteststan:::stanmodels), 1L)
})

test_that("Data that cannot be fitted with nls_list/nlme work with stan_fit", {
  # with this seed, cf[10] does not fit with nls_list
  library(breathtestcore)
  chains = 1
  student_df = 10
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(seed = 100)$data)
  fit = stan_fit(data, dose = dose, student_df = student_df,
                 chains = chains, iter = iter  )
  expect_is(fit, "breathtestfit")
  expect_is(fit, "breathteststanfit")
  expect_is(fit$stan_fit, "stanfit" )
})

test_that("A single record can be fitted", {
  library(breathtestcore)
  chains = 1
  student_df = 10
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(n_records = 1, seed = 100)$data)
  fit = stan_fit(data, dose = dose, student_df = student_df,
                 chains = chains, iter = iter  )
  expect_is(fit, "breathtestfit")
})

