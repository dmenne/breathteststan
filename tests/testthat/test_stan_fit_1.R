test_that("stanmodels exist", {
  expect_type(breathteststan:::stanmodels,"list")
  expect_s4_class(breathteststan:::stanmodels$breath_test_1, "stanmodel")
  expect_gte(length(breathteststan:::stanmodels), 1L)
})

