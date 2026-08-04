library(testit)
library(breathtestcore)
library(breathteststan)

assert("Result with default parameters is tbl_df with required columns", {
  # This calls coef_by_group.breathtestfit, which include a post-hoc classic test
  # for contrasts.
  data("usz_13c", package = "breathtestcore")
  data = usz_13c |>
    dplyr::filter(
      patient_id %in%
        c(
          "norm_001",
          "norm_002",
          "norm_004",
          "norm_007",
          "pat_001",
          "pat_002",
          "pat_003"
        )
    ) |>
    cleanup_data()
  comment(data) = "comment"
  fit = stan_fit(data, iter = 500, chains = 1)
  (inherits(fit, "breathtestfit"))
  (inherits(fit, "breathteststanfit"))
  (inherits(fit$stan_fit, "CmdStanMCMC"))
  (inherits(fit$stan_fit, "CmdStanFit"))

  #  class(fit) = class(fit)[-1] # Remove class breathteststanfit
  cf = coef_by_group(fit) # S3 method
  #https://yihui.org/en/2026/05/testthat-to-testit/
  (inherits(cf, "tbl_df"))
  (inherits(cf, "coef_by_group"))
  (identical(ncol(cf), 7L))
  (all.equal(
    names(cf),
    c(
      "parameter",
      "method",
      "group",
      "estimate",
      "conf.low",
      "conf.high",
      "diff_group"
    )
  ))
  (identical(comment(data), "comment"))
  (identical(nrow(cf), 24L))
  print(sort(unique(cf$diff_group)))
  (identical(sort(unique(cf$diff_group)), c("a", "b", "c")))
  (identical(
    unique(cf$group),
    c("liquid_normal", "solid_normal", "solid_patient")
  ))
})
