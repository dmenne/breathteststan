library(testthat)
library(breathteststan)

test_check("breathteststan", filter = "coef_diff_by_group")
#test_check("breathteststan")
