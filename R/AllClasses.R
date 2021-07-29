#' MR elastic net class
#'
#' @description
#' The output class of the glmnet_enet_mr() function. This code was written by Dr. Verena Zuber, <v.zuber(at)imperial.ac.uk>
#' who has kindly given us permission to incorporate it into the PathWAS package.
#' The original code for this can be found here: https://github.com/verena-zuber/demo_AMD/blob/master/mvMR_glmnet.R
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Estimate The causal estimates from the multivariable MR-Lasso method.
#' @slot Lambda1 The value of the "best" lambda tuning parameter used to compute elasticnet.
#' @slot Lambda2 The "best" alpha value. With the alpha value being the elasticnet mixing parameter, with 0 ≤ α ≤ 1.
#'
setClass("MRenet",
         representation(Exposure = "character",
                        Outcome = "character",
                        Estimate = "numeric",
                        Lambda1 = "numeric",
                        Lambda2 = "numeric")
)
