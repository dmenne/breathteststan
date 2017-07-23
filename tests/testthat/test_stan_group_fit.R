context("Bayesian fit of multiple groups")

test_that("Exception when there is only one group", {
  data = breathtestcore::cleanup_data(simulate_breathtest_data(seed = 100)$data)
  expect_error(stan_group_fit(data), "only one group")

})


test_that("Multiple records per patient return multiple groups (CRAN version)", {
  data("usz_13c", package = "breathtestcore")
  set.seed(4711)
  data = usz_13c %>%
    dplyr::filter( patient_id %in%
                     c("norm_001", "norm_002", "norm_003")) %>%
    breathtestcore::cleanup_data()
  comment(data) = "comment"
  fit = stan_group_fit(data, iter = 300)
  expect_identical(comment(fit), "comment")
  expect_is(fit, "breathteststangroupfit")
  expect_identical(names(fit), c("coef", "data", "stan_fit"))
})

test_that("Multiple records per patient return multiple groups (long version)", {
  skip_on_cran() # long

#  library(breathtestcore)
#  library(breathteststan)
  library(dplyr)

  chains = 2
  student_t_df = 3  # Rough student distribution
  dose = 100
  iter = 500
  model = "breath_test_group_1"
  data("usz_13c", package = "breathtestcore")
  set.seed(4711)
  data = usz_13c %>%
    dplyr::filter( patient_id %in%
                     c("norm_001", "norm_002", "norm_003",
                       "pat_001", "pat_003", "pat_016")) %>%
    breathtestcore::cleanup_data()
  # fit nlme for comparison
  fit_nlme = breathtestcore::nlme_fit(data)
  # fit stan_group
  fit = stan_group_fit(data, dose = dose, student_t_df = student_t_df,
                       chains = chains, iter = iter, model = model  )

  cf = coef(fit)
  expect_equal(unique(cf$group), c("liquid_normal", "solid_normal", "patient"))
  expect_gt(sigma(fit), 0.5)

  cf = coef(fit) %>%
    left_join(coef(fit_nlme), by = c("patient_id", "parameter", "method", "group"))  %>%
    filter(method == "exp_beta") %>%
    mutate(rel_diff = 2*abs(value.x - value.y)/(value.x + value.y))   %>%
    select(parameter, rel_diff) %>%
    summarize(
      rel_diff = mean(rel_diff)
    )
  expect_lt(cf$rel_diff, 0.06)
})


