# Gaussian Single-Index Engression
SIE <- function(
    X,
    Y,
    noise_dim = 5,
    hidden_dim = 100,
    num_layer = 3,
    dropout = 0.05,
    batch_norm = TRUE,
    num_epochs = 50,
    lr = 10^(-3),
    beta = 1,
    silent = TRUE,
    standardize = TRUE,
    GSI_rep = 1, 
    quantiles = NULL, 
    return_index = FALSE, 
    predict_from_index = NULL,
    mean_prediction = FALSE, 
    test_data = NULL,
    verbose = FALSE,
    icp_test = NULL,
    E = NULL, 
    shapiro = FALSE,
    init_method = "linear", # "gaussian_quantile", 
    return_yinv = FALSE, 
    monotone = FALSE){
  
  if (!is.null(icp_test)) {
    implemented_tests <- c("leveneAndWilcoxResidualDistributions")
    if(!(icp_test %in% implemented_tests) ) {
      stop("Please specify valid ICP test. Options: leveneAndWilcoxResidualDistributions. Otherwise set to 'NULL'.")
    }
    if(is.null(E)) {
      stop("Environment is missing.")
    }
    
    if(!is.factor(E)){
      stop("InvariantResidualDistributionTest can only be applied if E is a factor.")
    }
    
    if(NCOL(E) > 1){
      stop("InvariantResidualDistributionTest can only be applied if E is univariate.")
    }
  }
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  X_mean <- apply(X,2,mean)
  X_sd <- apply(X,2,sd)
  
  Y_mean <- apply(Y,2,mean)
  Y_sd <- apply(Y,2,sd)
  
  # X <- scale(X, center = TRUE, scale = FALSE)
  # Y <- scale(Y, center = FALSE, scale = FALSE)
  
  X  <- sweep(sweep(X,2,X_mean,FUN="-"),2,X_sd,FUN="/")
  Y <- sweep(sweep(Y,2,Y_mean,FUN="-"),2,Y_sd,FUN="/")

  n <- nrow(X)
  p <- ncol(as.matrix(X))
  
  # initialize g and beta
  init_output <-  init_beta(X, 
                            Y, 
                            init_method = init_method,
                            hidden_dim = hidden_dim,
                            num_layer = num_layer,
                            dropout = dropout,
                            batch_norm = batch_norm,
                            num_epochs = num_epochs,
                            lr = lr,
                            beta = beta,
                            silent = silent,
                            standardize = standardize)
  beta_init <- init_output$beta_init
  index_hat <- as.matrix( init_output$index_init )
  
  # train gauissian single index model
  for (rep in 1:GSI_rep) {
    if (monotone) {
      GSIfit  <- monotone_engression(index_hat,
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
                            standardize = standardize )
    } else {
      GSIfit  <- engression(index_hat,
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
                            standardize = standardize )
    }
    
    index_range <- diff(range(index_hat))
    index_grid <- seq(min(index_hat) - 0.2*index_range, 
                max(index_hat) + 0.2*index_range, 
                by=0.01)
    
    if (mean_prediction){
      yhat <- predict(GSIfit, index_grid)
    } else {
      yhat <- predict(GSIfit, index_grid, type="quantile", quantiles = 0.5)
    }
    
    monotone_data <- monotone(index_grid,yhat,drop = TRUE)
    index_grid <- monotone_data$X
    yhat <- monotone_data$Y
    
    apply_inverse <- estimate_inverse(index_grid, yhat, method = "linear")
    
    yinv_hat <- apply_inverse(Y)
    
    LMfit <- lm(scale(yinv_hat) ~ X)
    index_hat <- as.matrix( LMfit$fitted.values )
  }
  
  shapiro_pval <- NULL
  if (shapiro) {
    shapiro_pval <- shapiro.test(scale(yinv_hat)-index_hat)$p.value
  }
  
  icp_result <- NULL
  if (!is.null(icp_test)) {
    # test whether residual distribution is identical in all environments E
    icp_result <- leveneAndWilcoxResidualDistributions(Y = scale(yinv_hat),  
                                                       predicted = index_hat, E, verbose)
  }
  
  # return index coefficients 
  beta_hat <- round(summary(LMfit)$coefficients[-1, 1],3)
  sd_hat <- round(summary(LMfit)$coefficients[-1, 2],3)
  
  predicted <- NULL
  if (!is.null(test_data)) {
    GSIfit  <- engression(index_hat,
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
                          standardize = standardize)
    test_data <- as.matrix(test_data)
    test_data  <- sweep(sweep(test_data,2,X_mean,FUN="-"),2,X_sd,FUN="/")
    
    if (!is.null(quantiles))  {
      predicted <- predict(GSIfit, as.matrix( test_data ) %*% as.vector( beta_hat ), 
                           type="quantile", quantiles = quantiles)
    } else if (mean_prediction){
      predicted <- predict(GSIfit, as.matrix( test_data ) %*% as.vector( beta_hat ))
    } else if (!mean_prediction) {
      predicted <- predict(GSIfit, as.matrix( test_data ) %*% as.vector( beta_hat ), 
                           type="quantile", quantiles = 0.5)
    }
    
    predicted <- sweep(sweep(as.matrix( predicted ),2, Y_sd,FUN="*"),2,Y_mean,FUN="+")
  } else {
    if (!is.null(quantiles))  {
      GSIfit  <- engression(index_hat,
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
                            standardize = standardize)
      predicted <- predict(GSIfit, as.matrix( X ) %*% as.vector( beta_hat ), 
                           type="quantile", quantiles = quantiles)
      
    } else if (mean_prediction){
      GSIfit  <- engression(index_hat,
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
                            standardize = standardize)
      predicted <- predict(GSIfit, as.matrix( X ) %*% as.vector( beta_hat ))
      
    } else if (!mean_prediction) {
      GSIfit  <- engression(index_hat,
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
                            standardize = standardize)
      predicted <- predict(GSIfit, as.matrix( X ) %*% as.vector( beta_hat ), 
                           type="quantile", quantiles = 0.5)
      
    }
    
    predicted <- sweep(sweep(as.matrix( predicted ),2, Y_sd,FUN="*"),2,Y_mean,FUN="+")
  }
  
  predicted_from_index <- NULL
  if (!is.null(predict_from_index)) {
    GSIfit  <- engression(index_hat,
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
                          standardize = standardize)
    predict_from_index <- as.matrix(predict_from_index)
    # predict_from_index  <- sweep(sweep(predict_from_index,2,X_mean,FUN="-"),2,X_sd,FUN="/")
    
    if (!is.null(quantiles))  {
      predicted_from_index <- predict(GSIfit, predict_from_index, 
                           type="quantile", quantiles = quantiles)
    } else if (mean_prediction){
      predicted_from_index <- predict(GSIfit, predict_from_index)
    } else if (!mean_prediction) {
      predicted_from_index <- predict(GSIfit, predict_from_index, 
                           type="quantile", quantiles = 0.5)
    }
    
    predicted_from_index <- sweep(sweep(as.matrix( predicted_from_index ),2, Y_sd,FUN="*"),2,Y_mean,FUN="+")
  }
  
  
  index <- NULL
  if (return_index) {
    index <- index_hat
  }
  
  yinv <- NULL
  if (return_yinv) {
    yinv <- scale(yinv_hat)
  }
  
  list(predicted = predicted, beta_hat = beta_hat, sd_hat = sd_hat, 
       icp_result = icp_result, shapiro_pval = shapiro_pval, index = index, 
       yinv = yinv, predicted_from_index=predicted_from_index)
}

















