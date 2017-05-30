context("Coeffients by group")

# TODO: this calls coef_by_group.breathtestfit, which include a post-hoc classic test
# for contrasts. Should be replaced by a full Bayesian estimate using the posterior

test_that("Result with default parameters is tbl_df with required columns",{
  skip_on_cran()
  library(dplyr)
  library(breathtestcore)
  data("usz_13c")
  data = usz_13c %>%
    dplyr::filter( patient_id %in%  c("norm_001", "norm_002", "norm_003",
            "norm_004", "pat_001", "pat_002","pat_003")) %>%
    cleanup_data()
  fit = stan_fit(data)
  cf = coef_by_group(fit) # S3 method
  expect_is(cf, "tbl_df")
  expect_identical(ncol(cf), 7L)
  expect_equal(names(cf), c("parameter", "method", "group", "estimate", "conf.low",
                 "conf.high", "diff_group"))
  expect_identical(nrow(cf), 24L)
  expect_identical(unique(cf$diff_group), c("a", "b", "ab", "c"))
  expect_equal(unique(cf$group),
     c("liquid_normal", "patient", "solid_normal"))
})

