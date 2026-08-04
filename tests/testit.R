# tests/testit.R
library(testit)
test_pkg("breathteststan", dir = "cran")

# tests under `slow/`; only run when not on CRAN
if (identical(tolower(Sys.getenv("NOT_CRAN")), "true")) {
  test_pkg("breathteststan", dir = "slow")
}

# tests under `ci/`; only run when on ci (generic)
# Might also check for "GITHUB_ACTIONS"
#if (identical(tolower(Sys.getenv("CI")), "true")) {
#  test_pkg("breathteststan", dir = "ci")
#}
