#' Computes ecological weirdness for a fitted model or it estimated  coefficients
#'
#' @param x Bag.
#' @param measure `[string(1)]{""coef_sign", "n_crosses", "response_auc""}` \cr Measure used
#' to quantify "weirdness" in the model or coefficients, based on the coeffcients and the response
#' plots for each type of covariate with zone of influence in a model.
#' It can be one or multiple of these options:
#' - `"coef_sign"`: counts the number of coefficients whose signal is opposite to the ecologically expected signal;
#' - `"n_crosses"`: counting the number of crosses in signal for the coefficients of the same covariate;
#' - `"response_auc"`: computing the area under the response plot curve which is in the unexpected
#' direction.
#' @param expected_sign `[numeric(1)=-1]` \cr Expected sign of the coefficient. Either -1 (negative),
#' +1 (positive), or 0 (no effect).
#' @param zero_coefficient_limit `[numeric(1)=1e8]` \cr Value above which an estimated coefficient is considered
#' non-zero. Default is 1e-8. Useful for comparing coefficients which are expected to be zero (i.e. to have no effect).
#'
#' @example weirdness_example.R
#'
weirdness <- function(x,
                      measure = c("coef_sign", "n_crosses", "response_auc"),
                      expected_sign = -1,
                      zero_coefficient_limit = 1e-8,
                      radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
                      type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
                               "mfilter")[1],
                      radius_max = 10000) {
  UseMethod("weirdness")
}

# weirdness.bag <- function(x,
#                       measure = c("coef_sign", "n_crosses", "response_auc"),
#                       expected_sign = -1,
#                       zero_coefficient_limit <- 1e-8,
#                       radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
#                       type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
#                                "mfilter")[1],
#                       radius_max = 10000) {
#
#
#   # coefs
#   coefs <- x$coef
#
#   weirdness_measures <- list(n_coefs = length(coefficients),
#                              coef_sign = NULL,
#                              n_crosses = NULL,
#                              response_auc = NULL)
#
#   # compute measure for the signal of coefficients
#   if("coef_sign" %in% measure) {
#     weirdness_measures$coef_sign <- ifelse(expected_sign == 0,
#                                            abs(coefficients) > zero_coefficient_limit,
#                                            sum(coefficients*expected_sign < 0))
#   }
#
#   # possibility: sum of absolute values of standardized coefficients that are against the expected sign
#
#   # number of crosses for the response plot
#
#
# }
#
# # x = vector of coefficients here
# weirdness <- function(x,
#                       measure = c("coef_sign", "n_crosses", "response_auc"),
#                       expected_sign = -1,
#                       zero_coefficient_limit <- 1e-8,
#                       radii = c(100, 250, 500, 1000, 2500, 5000, 10000),
#                       type = c("circle", "Gauss", "rectangle", "exp_decay", "bartlett", "threshold",
#                                "mfilter")[1],
#                       radius_max = 10000) {
#
#
#   weirdness_measures <- list(n_coefs = length(x),
#                              coef_sign = NULL,
#                              n_crosses = NULL,
#                              response_auc = NULL)
#
#   # compute measure for the signal of coefficients
#   if("coef_sign" %in% measure) {
#     weirdness_measures$coef_sign <- ifelse(expected_sign == 0,
#                                            abs(x) > zero_coefficient_limit,
#                                            sum(x*expected_sign < 0))
#   }
#
#   # possibility: sum of absolute values of standardized coefficients that are against the expected sign
#
#   # number of crosses for the response plot
#
#
# }
