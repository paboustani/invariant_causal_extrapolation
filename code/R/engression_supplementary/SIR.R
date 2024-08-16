SIR <- function(X, Y, H = 10) {
  X <- scale(as.matrix(X))
  Y <- scale(Y)
  
  # Add a small jitter to Y to ensure unique quantile breaks
  Y_jittered <- Y + rnorm(length(Y), mean = 0, sd = 1e-6)
  
  # Slice the response variable
  slices <- cut(Y_jittered, breaks = quantile(Y_jittered, probs = seq(0, 1, length.out = H + 1)), include.lowest = TRUE)
  
  # Calculate the means of X in each slice
  slice_means <- sapply(levels(slices), function(slice) colMeans(X[slices == slice, , drop = FALSE]))
  
  # Compute the covariance matrix of the slice means
  cov_slice_means <- cov(t(slice_means))
  
  # Eigen decomposition of the covariance matrix
  eig <- eigen(cov_slice_means)
  
  # Select the leading eigenvector
  beta_hat <- round(eig$vectors[,1], 4)
  
  # Predict the index from the input X and beta_init via linear prediction
  index_hat <- X %*% beta_hat
  
  return(list(beta_hat = beta_hat, index_hat = index_hat))
}
