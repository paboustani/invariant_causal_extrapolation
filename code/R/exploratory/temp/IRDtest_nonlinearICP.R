InvResDisTest <- function(Y, XS, E, alpha = 0.05, test) {
  # Pool the data from all environments and fit a model to predict Y with XS
  fitted_model <- gam(Y ~ s( XS ))
  
  pv <- 1
  t <- 0
  
  for (e in unique(E)) {
    # Use a two-sample test to assess whether the residuals of samples from environment e
    # have the same distribution as the residuals of samples from environments in the index set E′
    E_prime <- E[E != e]  # E'
    
    # Perform the test (subroutine)
    residuals_e <- fitted_model$residuals[E == e]
    residuals_e_prime <- fitted_model$residuals[E %in% E_prime]
    
    if (test == "ks") {
      # Use Kolmogorov-Smirnov test
      pv_e <- ks.test(residuals_e, residuals_e_prime)$p.value
    } else if (test == "leveneAndWilcox") {
      # Use Wilcoxon test for expectation
      pv_e <- wilcox.test(residuals_e, residuals_e_prime)$p.value
    } else {
      stop("Invalid test specified.")
    }
    
    t <- t + 1
    pv <- min(pv, pv_e)
    
    if (length(unique(E)) == 2) {
      break
    }
  }
  
  # Apply a Bonferroni correction for the number of performed tests t
  pv <- t * pv
  
  if (test == "leveneAndWilcox") {
    # Levene's test for homogeneity of variance
    pv_l <- levene.test(fitted_model$residuals, E)$p.value
    pv <- min(pv, pv_l)
  }
  
  return(pv)
}