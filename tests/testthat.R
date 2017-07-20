library(testthat)
library(breathteststan)

test_check("breathteststan", filter = "stan_group_fit")
