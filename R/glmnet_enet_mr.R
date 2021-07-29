#' Run elastic net-based mendelian randomisation from input MR object
#'
#' @description
#' Lower level function. Takes input "MR input" object (created by mr_mvinput()) and performs elastic net mendelian
#' randomisation with this input. Outputs a model of exposures and betas.
#'
#' @details
#' This function was written by Dr. Verena Zuber, <v.zuber(at)imperial.ac.uk> who has kindly given us permission
#' to incorporate it into the PathWAS package.
#' The original code for this can be found here: https://github.com/verena-zuber/demo_AMD/blob/master/mvMR_glmnet.R
#' It works in a similar method to the MendelianRandomization R package in working from the same input object
#' (a matrix of betas and SEs for every SNP used with each exposure).
#'
#' @param object MR Input object. Contains a matrix of exposure SNPs to betas, exposure SNPs to SEs and then a list of outcome SNP betas and SEs.
#'

glmnet_enet_mr = function(object,
                          cv = TRUE,
                          lambda = 0.1,
                          alpha = 0.1,
                          cv.param = "lambda.1se"
                          ){


  bX = object@betaX
  bY = object@betaY

  if (cv==TRUE){

    a = seq(0.1, 0.9, 0.1)

    if(cv.param == "lambda.1se"){

      search = foreach(i = a, .combine = rbind) %do% {

        cv = cv.glmnet(bX, bY,
                       family = "gaussian",
                       nfold = 10,
                       type.measure = "mse",
                       intercept = FALSE,
                       alpha = i,
                       standardize = FALSE
                       )
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)

      }

      cv = search[search$cvm == min(search$cvm), ]
      bestlambda = cv$lambda.1se
      bestalpha = cv$alpha

    }

    if(cv.param=="lambda.min"){

      search = foreach(i = a, .combine = rbind) %do% {
        cv = cv.glmnet(bX, bY,
                       family = "gaussian",
                       nfold = 10,
                       type.measure = "mse",
                       intercept = FALSE,
                       alpha = i,
                       standardize = FALSE
                       )
        data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
      }
      cv = search[search$cvm == min(search$cvm), ]
      bestlambda = cv$lambda.min
      bestalpha = cv$alpha

    }

    enet.out = glmnet(bX, bY,
                      family = "gaussian",
                      intercept = FALSE,
                      lambda = bestlambda,
                      alpha = bestalpha,
                      standardize = FALSE
                      )
    enet.coeff = coef(enet.out)[2:(ncol(bX)+1)]
  }

  else{

    bestlambda = lambda
    bestalpha = alpha
    enet.out =  glmnet(bX, bY,
                       family = "gaussian",
                       intercept = FALSE,
                       lambda = bestlambda,
                       alpha = bestalpha,
                       standardize = FALSE
                       )
    enet.coeff = coef(enet.out)[2:(ncol(bX)+1)]

    }

  return(new("MRenet",
             Exposure = object@exposure,
             Outcome = object@outcome,
             Estimate =  enet.coeff,
             Lambda1 = bestlambda,
             Lambda2 = bestalpha
             )
         )
}
