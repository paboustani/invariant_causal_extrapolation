#' Invariant conditional quantile prediction.
#'
#' @description Tests the null hypothesis that Y and E are independent given X.
#'
#' @param Y An n-dimensional vector.
#' @param E An n-dimensional vector. If \code{test = fishersTestExceedance}, E needs
#' to be a factor.
#' @param X A matrix or dataframe with n rows and p columns.
#' @param alpha Significance level. Defaults to 0.05.
#' @param verbose If \code{TRUE}, intermediate output is provided. Defaults to \code{FALSE}.
#' @param test Unconditional independence test that tests whether exceedence is
#' independent of E. Defaults to \code{fishersTestExceedance}.
#' @param mtry Random forest parameter: Number of variables randomly sampled as
#' candidates at each split. Defaults to \code{sqrt(NCOL(X))}.
#' @param ntree Random forest parameter: Number of trees to grow. Defaults to 100.
#' @param nodesize Random forest parameter: Minimum size of terminal nodes.  Defaults to 5.
#' @param maxnodes Random forest parameter: Maximum number of terminal nodes trees in the forest can have.
#' Defaults to NULL.
#' @param quantiles Quantiles for which to test independence between exceedence and E.
#' Defaults to \code{c(0.1, 0.5, 0.9)}.
#' @param returnModel If \code{TRUE}, the fitted quantile regression forest model
#' will be returned. Defaults to \code{FALSE}.
#'
#' @return A list with the following entries:
#' \itemize{
#'  \item \code{pvalue} The p-value for the null hypothesis that Y and E are independent given X.
#'  \item \code{model} The fitted quantile regression forest model if \code{returnModel = TRUE}.
#'  }
#'
#' @examples
#' # Example 1
#' n <- 1000
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * (X)^2 + rnorm(n)
#' InvariantConditionalQuantilePrediction(Y, as.factor(E), X)
#'
#' # Example 2
#' E <- rbinom(n, size = 1, prob = 0.2)
#' X <- 4 + 2 * E + rnorm(n)
#' Y <- 3 * E + rnorm(n)
#' InvariantConditionalQuantilePrediction(Y, as.factor(E), X)
#'
InvariantConditionalQuantilePrediction <- function(Y, E, X,
                                                  fitmodel = "RF", # also: engression, SIE
                                                  alpha = 0.05,
                                                  verbose = FALSE,
                                                  test = fishersTestExceedance,
                                                  mtry = sqrt(NCOL(X)),
                                                  ntree = 100,
                                                  nodesize = 5,
                                                  maxnodes = NULL,
                                                  quantiles = c(0.1, 0.5, 0.9),
                                                  returnModel = FALSE, 
                                                  noise_dim = 5, # engression specifications
                                                  hidden_dim = 100,
                                                  num_layer = 3,
                                                  dropout = 0.05,
                                                  batch_norm = TRUE,
                                                  num_epochs = 50,
                                                  lr = 10^(-3),
                                                  beta = 1,
                                                  silent = TRUE,
                                                  standardize = TRUE,
                                                  var_names = paste("Z", 1:6, sep = "_") ){

  Y <- check_input_single(Y, return_vec = TRUE, str = "Y")
  E <- check_input_single(E, check_factor = TRUE, return_vec = TRUE, str = "E")
  X <- check_input_single(X, return_vec = FALSE)
  
  
  n <- NROW(X)
  p <- NCOL(X)

  # train model Y ~ X using all data
  mat <- as.matrix(X)
  colnames(mat) <- paste("V", 1:ncol(mat), sep = "")
  
  if(fitmodel=="RF"){
    # fit model
    rfResult <- quantregForest(x = mat,
                               y = Y,
                               mtry = mtry,
                               ntree = ntree,
                               nodesize = nodesize,
                               maxnodes = maxnodes)
    
    # predict
    predicted <- predict(rfResult, newdata = mat, what = quantiles)
  } else if(fitmodel=="engression") {
    # fit model
    engressionResult <- engression( X = mat,
                                    Y = Y,
                                    noise_dim = noise_dim,
                                    hidden_dim = hidden_dim,
                                    num_layer = num_layer,
                                    dropout = dropout,
                                    batch_norm = batch_norm,
                                    num_epochs = num_epochs,
                                    lr = lr,
                                    beta = beta,
                                    silent = silent,
                                    standardize = standardize)
    
    # predict
    predicted <- predict(engressionResult, mat, type="quantile", quantiles = quantiles)
  }  else if(fitmodel=="SIE") {
    # fit model
    res <- SIE(X, Y, noise_dim = noise_dim, hidden_dim, num_layer, num_epochs, silent, GSI_rep = 1, quantiles = quantiles)
    predicted <- res$predicted
  } else {
    stop("Invalid model, please specify one of the following: RF, engression.")
  }

  # test whether residual distribution is identical in all environments E
  result <- test(Y, predicted, E, verbose)
  
  if(fitmodel=="SIE"){
    betavec <- setNames(rep(0, length(var_names)), paste("beta", var_names, sep = "_"))
    names(res$beta_hat) <- paste0("beta_",colnames(as.matrix(X)))
    betavec[names(res$beta_hat)] <- res$beta_hat
    result$beta_hat <- betavec
    
    sdvec <- setNames(rep(0, length(var_names)), paste("sd", var_names, sep = "_"))
    names(res$sd_hat) <- paste0("sd_",colnames(as.matrix(X)))
    sdvec[names(res$sd_hat)] <- res$sd_hat
    result$sd_hat <- sdvec
  }

  if(returnModel){
    result$model <- list(rfResult = rfResult)
  }

  result
}

