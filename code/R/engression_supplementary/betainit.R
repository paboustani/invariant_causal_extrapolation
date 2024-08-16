init_beta <- function(X, 
                      Y, 
                      init_method = "gaussian_quantile",
                      hidden_dim = 100,
                      num_layer = 3,
                      dropout = 0.05,
                      batch_norm = TRUE,
                      num_epochs = 50,
                      lr = 10^(-3),
                      beta = 1,
                      silent = TRUE,
                      standardize = TRUE, 
                      H = 10){
  
  X <- scale(X)
  Y <- scale(Y)
  
  n <- nrow(X)
  p <- ncol(as.matrix(X))
  
  if (init_method=="gaussian_quantile") {
    G_quantiles <- seq(0.001, 0.999, by = 0.001)
    
    GSIinitfit <- gaussianengression( X = X, 
                                      Y = Y, 
                                      noise_dim = 0, 
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
                           matrix(0, nrow = 1, ncol = p),  
                           type="quantile", 
                           quantiles= G_quantiles )
    
    monotone_data <- monotone(G_quantiles,qantile_hat,drop = TRUE)
    G_quantiles <- monotone_data$X
    qantile_hat <- monotone_data$Y
    
    apply_inverse <- estimate_inverse(G_quantiles, qantile_hat)
    
    yinv_init <- apply_inverse(Y)
    
    LMinitfit <- lm(scale( yinv_init ) ~ X)
    beta_init <- summary(LMinitfit)$coefficients[-1, 1]
    beta_init <- round(beta_init/norm(beta_init, type = "2"),4)
    index_init <- X %*% beta_init
    # index_init <- LMinitfit$fitted.values
  
    } else if(init_method=="SIR"){
    sir_output <- SIR(X, Y, H = H)
    beta_init <- sir_output$beta_hat
    index_init <- sir_output$index_hat
    
  } else if(init_method=="linear"){
    LMinitfit <- lm(Y ~ X)
    beta_init <- summary(LMinitfit)$coefficients[-1, 1]
    beta_init <- round(beta_init/norm(beta_init, type = "2"),4)
    index_init <- X %*% beta_init
    # index_init <- LMinitfit$fitted.values
    
  } else {
    stop("Invalid initialization method, please specify one of the following: gaussian_quantile, SIR, linear.")
  }
  
  list(index_init = index_init, beta_init = beta_init)
}


