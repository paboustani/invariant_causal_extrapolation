# # Function to create the inverse function
# estimate_inverse <- function(x_data, y_data) {
#   # Ensure data is sorted
#   sorted_indices <- order(y_data)
#   x_data_sorted <- x_data[sorted_indices]
#   y_data_sorted <- y_data[sorted_indices]
#   
#   # Generalized inverse function
#   apply_inverse <- function(y_values) {
#     # Initialize a vector to store the results
#     x_est_values <- numeric(length(y_values))
#     
#     # Compute the inverse for each y_value in y_values
#     for (i in seq_along(y_values)) {
#       y_value <- y_values[i]
#       if (y_value <= min(y_data_sorted)) {
#         x_est_values[i] <- min(x_data_sorted)
#       } else if (y_value >= max(y_data_sorted)) {
#         x_est_values[i] <- max(x_data_sorted)
#       } else {
#         x_est_values[i] <- x_data_sorted[which(y_data_sorted >= y_value)[1]]
#       }
#     }
#     
#     return(x_est_values)
#   }
#   
#   # Return the invert function
#   return(apply_inverse)
# }

monotone <- function(X, Y, drop = TRUE) {
  # Check if the relation is increasing or decreasing
  if (is.na(cor(X, Y))) {
    corr <- 0
  } else {
    corr <- cor(X, Y)
  }
  is_increasing <- corr >= 0
  
  # Create a data frame and sort it based on X
  data <- data.frame(X, Y)
  data_sorted <- data[order(data$X), ]
  
  if (drop) {
    if (is_increasing) {
      # Implementing the "drop" for increasing relationship
      Y_dropped <- data_sorted$Y
      to_keep <- rep(TRUE, length(Y_dropped))
      max_Y <- Y_dropped[1]
      for (i in 2:length(Y_dropped)) {
        if (Y_dropped[i] <= max_Y) {
          to_keep[i] <- FALSE
        } else {
          max_Y <- Y_dropped[i]
        }
      }
      data_sorted <- data_sorted[to_keep, ]
    } else {
      # Implementing the "drop" for decreasing relationship
      Y_dropped <- data_sorted$Y
      to_keep <- rep(TRUE, length(Y_dropped))
      min_Y <- Y_dropped[1]
      for (i in 2:length(Y_dropped)) {
        if (Y_dropped[i] >= min_Y) {
          to_keep[i] <- FALSE
        } else {
          min_Y <- Y_dropped[i]
        }
      }
      data_sorted <- data_sorted[to_keep, ]
    }
  } else {
    # Apply the appropriate step function based on the relation
    if (is_increasing) {
      # Implementing the "lifting" step function for increasing relationship
      Y_lifted <- data_sorted$Y
      for (i in 2:length(Y_lifted)) {
        if (Y_lifted[i] < Y_lifted[i - 1]) {
          Y_lifted[i] <- Y_lifted[i - 1]
        }
      }
      data_sorted$Y <- Y_lifted
    } else {
      # Implementing the "reducing" step function for decreasing relationship
      Y_reduced <- data_sorted$Y
      for (i in 2:length(Y_reduced)) {
        if (Y_reduced[i] > Y_reduced[i - 1]) {
          Y_reduced[i] <- Y_reduced[i - 1]
        }
      }
      data_sorted$Y <- Y_reduced
    }
  }
  
  data_sorted
}

# # Function to compute the generalized inverse
# estimate_inverse <- function(X, Y, method = "linear") {
#   # Check if the relation is increasing or decreasing
#   correlation <- cor(X, Y)
#   is_increasing <- correlation > 0
#   
#   # Ensure Y is strictly monotone in the correct direction
#   if (is_increasing) {
#     if (any(diff(Y) <= 0)) {
#       stop("Y is not strictly increasing.")
#     }
#   } else {
#     if (any(diff(Y) >= 0)) {
#       stop("Y is not strictly decreasing.")
#     }
#   }
#   
#   if (method == "linear") {
#     # Create a function for the generalized inverse using linear interpolation
#     apply_inverse <- approxfun(Y, X, method = "linear", rule = 2)
#   } else if (method == "stepwise") {
#     # Create a stepwise function for the generalized inverse
#     apply_inverse <- function(y) {
#       # Initialize the output vector
#       x_values <- numeric(length(y))
#       
#       # Loop through each y value to find the corresponding x
#       for (i in seq_along(y)) {
#         if (is_increasing) {
#           # Find the index of the largest Y that is less than or equal to the current y
#           index <- max(which(Y <= y[i]))
#         } else {
#           # Find the index of the smallest Y that is greater than or equal to the current y
#           index <- min(which(Y <= y[i]))
#         }
#         x_values[i] <- X[index]
#       }
#       
#       return(x_values)
#     }
#   } else if (method == "spline") {
#     # Fit a monotone spline based on the relationship
#     if (is_increasing) {
#       fit <- scam(X ~ s(Y, bs = "mpi", k = length(Y)))
#     } else {
#       fit <- scam(X ~ s(Y, bs = "mpd", k = length(Y)))
#     }
#     
#     # Create a function for the generalized inverse using the spline fit
#     apply_inverse <- function(y) {
#       predict(fit, newdata = data.frame(Y = y))
#     }
#   } else {
#     stop("Unknown method. Use 'linear', 'stepwise', or 'spline'.")
#   }
#   
#   return(apply_inverse)
# }


estimate_inverse <- function(X, Y, method = "linear") {
  # Check if the relation is increasing or decreasing
  if (is.na(cor(X, Y))) {
    correlation <- 0
  } else {
    correlation <- cor(X, Y)
  }
  is_increasing <- correlation >= 0
  
  # Ensure Y is strictly monotone in the correct direction
  if (is_increasing) {
    if (any(diff(Y) <= 0)) {
      stop("Y is not strictly increasing.")
    }
  } else {
    if (any(diff(Y) >= 0)) {
      stop("Y is not strictly decreasing.")
    }
  }
  
  if (method == "linear") {
    # Create a function for the generalized inverse using linear interpolation
    apply_inverse <- approxfun(Y, X, method = "linear", rule = 2)
  } else if (method == "stepwise") {
    # Create a stepwise function for the generalized inverse using stepfun
    apply_inverse <- stepfun(Y, c(X[1], X))
  } else if (method == "spline") {
    # Fit a monotone spline based on the relationship
    if (is_increasing) {
      fit <- scam(X ~ s(Y, bs = "mpi", k = length(Y)))
    } else {
      fit <- scam(X ~ s(Y, bs = "mpd", k = length(Y)))
    }
    
    # Create a function for the generalized inverse using the spline fit
    apply_inverse <- function(y) {
      predict(fit, newdata = data.frame(Y = y))
    }
  } else {
    stop("Unknown method. Use 'linear', 'stepwise', or 'spline'.")
  }
  
  return(apply_inverse)
}
