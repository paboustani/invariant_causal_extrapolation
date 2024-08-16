# Gaussian Single-Index Engression
gaussianengressionRD <- function(
    X,
    Y,
    noise_dim = 0,
    hidden_dim = 100,
    num_layer = 3,
    dropout = 0.05,
    batch_norm = TRUE,
    num_epochs = 50,
    lr = 10^(-3),
    beta = 1,
    silent = TRUE,
    standardize = TRUE,
    GSI_rep = 1
){

  n <- nrow(X)
  p <- ncol(as.matrix(X))
  
  # initialize g and beta
  quantiles <- seq(0.01, 0.99, by = 0.01)
  
  GSIinitfit <- gaussianengression( X = X, 
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
  
  qantile_hat <- predict(GSIinitfit, 
                         matrix(0, nrow = 1, ncol = 1),  
                         type="quantile", 
                         quantiles= quantiles )
  
  apply_inverse <- estimate_inverse(quantiles, qantile_hat)
  
  yinv_init <- apply_inverse(Y)
  
  LMinitfit <- lm(yinv_init ~ X)
  index_init <- LMinitfit$fitted.values
  
  index_hat <- index_init
  
  # train gauissian single index model
  for (rep in 1:GSI_rep) {
    GSIfit  <- gaussianengression(
      index_hat,
      Y,
      noise_dim = noise_dim,
      hidden_dim = hidden_dim,
      num_layer = num_layer,
      dropout = dropout,
      batch_norm = batch_norm,
      num_epochs = num_epochs,
      lr = lr,
      beta = beta,
      silent = silent,
      standardize = standardize
    )
    
    yhat <- predict(GSIfit, index_hat)
    
    apply_inverse <- estimate_inverse(index_hat, yhat)
    
    yinv_hat <- apply_inverse(Y)
    
    LMfit <- lm(yinv_hat ~ X)
    index_hat <- LMfit$fitted.values
  }
  
  # use trained GSI model to predict Y
  predicted <- predict(GSIfit, index_hat)
  # return index coefficients
  beta_hat <- summary(LMfit)$coefficients[-1, 1]
  sd_hat <- summary(LMfit)$coefficients[-1, 2]
  
  # if(!returnModel){
    list(predicted = predicted, beta_hat = beta_hat, sd_hat = sd_hat)
  # }else{
  #   list(predicted = predicted, beta_hat = beta_hat, model = engressionRes$model)
  # }
}
















