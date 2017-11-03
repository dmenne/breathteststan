library(testthat)
library(breathteststan)
options(warn = 2)
# Only one test per file to avoid hanging 32-bit compile
#test_check("breathteststan", filter = "stan_fit")
test_check("breathteststan")
