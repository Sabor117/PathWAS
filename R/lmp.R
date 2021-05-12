#' Extract P-value from regression model
#'
#' @description lmp returns the P-value of an input model
#'
#' @details A lower level function for outputting the P-value for a model, as this is not always immediately
#' extractable from just the one of the variables of the model.
#'
#' @param modelobject model output from lm function
#'
lmp = function (modelobject) {

  if (class(modelobject) != "lm"){

    stop("Not an object of class 'lm'. Womp womp.")

  }

  f = summary(modelobject)$fstatistic
  p = pf(f[1], f[2], f[3],
         lower.tail = F)

  attributes(p) = NULL
  return(p)

}
