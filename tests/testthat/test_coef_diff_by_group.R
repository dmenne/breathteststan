context("Coeffient differences by group")

# TODO: this calls coef_diff_by_group.breathtestfit, which include a post-hoc classic test
# for contrasts. Should be replaced by a full Bayesian estimate using the posterior


test_that("Result with default parameters is tbl_df with required columns",{
  skip_on_cran() # Slow
  library(dplyr)
  library(breathtestcore)
  data("usz_13c")
  data = usz_13c %>%
    dplyr::filter( patient_id %in%
     c("norm_001", "norm_002", "norm_003", "norm_004", "pat_001", "pat_002","pat_003")) %>%
    cleanup_data()
  fit = stan_fit(data)
  cf = breathtestcore::coef_diff_by_group(fit)
  expect_is(cf, "tbl_df")
  expect_identical(ncol(cf), 7L)
  expect_identical(nrow(cf), 24L)
  expect_lt(min(cf$p.value), 5e-10)
  expect_equal(unique(cf$groups),
     c("patient - liquid_normal", "solid_normal - liquid_normal", "solid_normal - patient"))
})

