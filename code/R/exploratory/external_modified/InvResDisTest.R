#' Invariant residual distribution test.
#' This is a modified version of the R file "InvariantResidualDistributionTest.R" 
#' by Heinze-Deml et al. (2018).
InvResDisTest <- function(Y, E, X,
                          alpha = 0.05,
                          verbose = FALSE,
                          fitmodel = "GAM", #  formally: fitWithGam
                          test = leveneAndWilcoxResidualDistributions,
                          colNameNoSmooth = NULL,
                          mtry = sqrt(NCOL(X)),
                          ntree = 100,
                          nodesize = 5,
                          maxnodes = NULL,
                          returnModel = FALSE){

  
  # Y <- check_input_single(Y, return_vec = TRUE, str = "Y")
  # E <- check_input_single(E, check_factor = TRUE, return_vec = TRUE, str = "E")
  # X <- check_input_single(X, return_vec = FALSE)
  
  if(!is.factor(E)){
    stop("InvariantResidualDistributionTest can only be applied if E is a factor.")
  }

  if(NCOL(E) > 1){
    stop("InvariantResidualDistributionTest can only be applied if E is univariate.")
  }

  n <- NROW(X)
  p <- NCOL(X)

  if(fitmodel=="GAM"){
    res <- gamResidualDistributions(X, Y, colNameNoSmooth, returnModel)
  } else if(fitmodel=="RF"){
    res <- rfResidualDistributions(X, Y, mtry, ntree, nodesize, maxnodes, returnModel)
  } else if(fitmodel=="engression"){
    eng_obj <- engression(X, Y)$residuals
    Yhat <- predict(eng_obj, X)
    res <- Y - Yhat
  } else {
    stop("Invalid model, pleas specify one of the following: GAM, RF, engression.")
  }

  # test whether residual distribution is identical in all environments E
  result <- test(Y, res$predicted, E, verbose)


  if(returnModel){
    result$model <- res$model
  }

  result
}

