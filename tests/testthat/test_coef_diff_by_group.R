context("Coeffient differences by group")

test_that("Credible intervals are returned as stan_group_fit coefficients",{
  skip_on_cran() # Slow
  library(dplyr)
  library(breathtestcore)
  data("usz_13c", package = "breathtestcore")
  data = usz_13c %>%
    dplyr::filter( patient_id %in%
     c("norm_001", "norm_002", "norm_003", "norm_004", "pat_001", "pat_002","pat_003")) %>%
    cleanup_data()
  comment(data) = "comment"
  fit = stan_group_fit(data, iter = 300)
  expect_identical(names(fit), c("coef", "data", "stan_fit", "coef_chain"))

})

