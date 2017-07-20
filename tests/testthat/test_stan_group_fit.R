context("Bayesian fit of multiple groups")

test_that("Exception when there is only one group", {
  data = cleanup_data(simulate_breathtest_data(seed = 100)$data)
  expect_error(stan_group_fit(data), "only one group")

})

test_that("Multiple records per patients return multiple groups", {
  skip_on_cran()
  library(breathteststan)
  library(breathtestcore)
  library(rstan)
  library(dplyr)
  library(rstan)
  library(stringr)
  library(testthat)
  chains = 1
  student_t_df = 10
  dose = 100
  iter = 500
  sample_minutes = 15
  model = "breath_test_group_1"
  data("usz_13c")
  set.seed(4711)
  data = usz_13c %>%
    dplyr::filter( patient_id %in%
                     c("norm_001", "norm_002", "norm_003")) %>%
    cleanup_data()
  fit = stan_group_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter, model = model  )

  if (FALSE){
    fit_nlme = nlme_fit(data, dose = dose)
    cf = coef(fit)
    expect_equal(unique(cf$group), c("liquid_normal", "solid_normal"))

    cf = coef(fit) %>%
      left_join(coef(fit_nlme), by = c("patient_id", "parameter", "method", "group")) %>%
      filter(method == "exp_beta") %>%
      mutate(rel_diff = 2*abs(value.x - value.y)/(value.x + value.y))   %>%
      select(parameter, rel_diff) %>%
      summarize(
        rel_diff = mean(rel_diff)
      )
    expect_lt(cf$rel_diff, 0.03)
  }
})


