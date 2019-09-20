context("Bayesian fit of multiple groups")

test_that("Exception when there is only one group", {
  data = breathtestcore::simulate_breathtest_data(seed = 100)$data
  data = breathtestcore::cleanup_data(data)
  expect_error(stan_group_fit(data), "only one group")

})


