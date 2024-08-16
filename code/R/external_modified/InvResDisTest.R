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
                          noise_dim = 5,
                          hidden_dim = 100,
                          num_layer = 3,
                          num_epochs = 1000,
                          silent = FALSE, 
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
    res <- engressionResidualDistributions(X, Y, noise_dim, hidden_dim, num_layer, num_epochs, silent)
  } else if(fitmodel=="GSIengression"){
    X <- scale(X, center = TRUE, scale = FALSE)
    Y <- scale(Y, center = TRUE, scale = FALSE)
    res <- gaussianengressionRD(X, Y, noise_dim = 0, hidden_dim, num_layer, num_epochs, silent, GSI_rep = 1)
  } else {
    stop("Invalid model, please specify one of the following: GAM, RF, engression.")
  }

  # test whether residual distribution is identical in all environments E
  result <- test(Y, res$predicted, E, verbose)
  
  if(fitmodel=="GSIengression"){
    betavec <- setNames(rep(0, 6), paste("beta_Z", 1:6, sep = "_"))
    names(res$beta_hat) <- paste0("beta_",colnames(as.matrix(X)))
    betavec[names(res$beta_hat)] <- res$beta_hat
    result$beta_hat <- betavec
  
    sdvec <- setNames(rep(0, 6), paste("sd_Z", 1:6, sep = "_"))
    names(res$sd_hat) <- paste0("sd_",colnames(as.matrix(X)))
    sdvec[names(res$sd_hat)] <- res$sd_hat
    result$sd_hat <- sdvec
  }
  
  if(returnModel){
    result$model <- res$model
  }

  result
}

