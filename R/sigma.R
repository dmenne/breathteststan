#' @title S3 method to extract the residual standard deviation
#' @description Functions for S3 method defined in breathtestcore for
#' \code{CmdStanFit}.
#' @param object A Stan-based fit returning an S6 method
#' @param ... Not used
#' @return A numeric value giving the sigma (= average residual standard deviation) term
#' from the Stan fit.
#' @importFrom stats sigma
#' @exportS3Method
sigma.CmdStanFit = function(object, ...) {
  object$stan_fit$summary("sigma", "mean")$mean
}
