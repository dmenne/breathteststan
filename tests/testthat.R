library(testthat)
Sys.unsetenv("R_TESTS") # https://github.com/r-lib/testthat/issues/603
options(Ncpus = parallelly::availableCores(omit = 1))

test_check("breathteststan")
