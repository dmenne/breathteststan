context("Bayesian fit with Stan")

test_that("stanmodels exist", {
  expect_is(breathteststan:::stanmodels,"list")
  expect_is(breathteststan:::stanmodels$breath_test_1, "stanmodel")
  expect_gte(length(breathteststan:::stanmodels), 1L)
})

test_that("Data that cannot be fitted with nls_list/nlme work with stan_fit", {
  # with this seed, cf[10] does not fit with nls_list
  library(breathtestcore)

#  library(rstan)
#  library(dplyr)
#  library(rstan)
#  library(stringr)
#  library(testthat)
  chains = 1
  student_t_df = 10
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(seed = 100)$data)
  fit = stan_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter  )
  expect_is(fit, "breathtestfit")
  expect_is(fit, "breathteststanfit")
  expect_is(fit$stan_fit, "stanfit" )
  cf = fit$coef
  expect_identical(unique(cf$group), "A")
  expect_identical(unique(cf$parameter), c("beta", "k", "m", "t50", "tlag"))
  expect_identical(unique(cf$stat), c("estimate", "q_0275", "q_25", "q_75", "q_975"))
  expect_equal(nrow(cf), 400)
  expect_equal(ncol(cf), 6)

  cf = coef(fit) # This is the "mean" group only
  expect_identical(unique(cf$group), "A")
  expect_identical(unique(cf$parameter), c("beta", "k", "m", "t50", "tlag"))
  expect_equal(nrow(cf), 80)
  expect_equal(ncol(cf), 5)
})

test_that("A single record can be fitted", {
  skip_on_cran()
  library(breathtestcore)
  chains = 1
  student_t_df = 10
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(n_records = 1, seed = 100)$data)
  fit = stan_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter  )
  expect_is(fit, "breathtestfit")
})

test_that("Multiple chains return valid results similar to nlme", {
  skip_on_cran()
  library(breathtestcore)
  chains = 2
  student_t_df = 10
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(n_records = 1, seed = 100)$data)
  fit = stan_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter  )
  fit_nlme = nlme_fit(data, dose = dose)
  cf = coef(fit) %>%
    left_join(coef(fit_nlme), by = c("patient_id", "parameter", "method")) %>%
    filter(method == "exp_beta") %>%
    mutate(rel_diff = 2*abs(value.x - value.y)/(value.x + value.y)) %>%
    select(parameter, rel_diff) %>%
    summarize(
      rel_diff = mean(rel_diff)
    )
  expect_lt(cf$rel_diff, 0.03)
})

test_that("Non-gaussian residuals with student_t_df work with outliers", {
  skip_on_cran()

  library(breathtestcore)
  chains = 1
  student_t_df = 3
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(n_records = 1, seed = 100,
                                         student_t_df = student_t_df)$data)
  fit = stan_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter  )
  fit_nlme = nlme_fit(data, dose = dose)
  cf = coef(fit) %>%
    left_join(coef(fit_nlme), by = c("patient_id", "parameter", "method")) %>%
    filter(method == "exp_beta") %>%
    mutate(rel_diff = 2*abs(value.x - value.y)/(value.x + value.y)) %>%
    select(parameter, rel_diff) %>%
    summarize(
      rel_diff = mean(rel_diff)
    )
  expect_lt(cf$rel_diff, 0.05)
})

