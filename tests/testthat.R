library(testthat)
Sys.unsetenv("R_TESTS") # https://github.com/r-lib/testthat/issues/603
options(Ncpus = parallel::detectCores(logical = TRUE))

test_check("breathteststan")
