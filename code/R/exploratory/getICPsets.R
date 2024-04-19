get_accepted <- function(pvals, S_sets, alpha = 0.05) {
  # Initialize an empty list to store the sets
  selected_sets <- list()
  
  # Iterate through each element of pvals
  for (i in seq_along(pvals)) {
    # Check if the p-value is not smaller than 0.05
    if (pvals[i] >= alpha) {
      # Add the corresponding set/vector from S_sets to the list
      selected_sets[[length(selected_sets) + 1]] <- S_sets[[i]]
    }
  }
  
  # Return the list of selected sets
  return(selected_sets)
}

get_union <- function(selected_sets) {
  # Initialize an empty set to store the union
  union_set <- NULL
  
  # Iterate through each set in selected_sets and update the union
  for (set in selected_sets) {
    union_set <- union(union_set, set)
  }
  
  # Return the union set
  return(union_set)
}

get_intersection <- function(selected_sets) {
  # Check if the selected_sets list is not empty
  if (length(selected_sets) == 0) {
    return(NULL)  # If the list is empty, return NULL
  }
  
  # Initialize the intersection set with the first set in selected_sets
  intersection_set <- selected_sets[[1]]
  
  # Iterate through each set in selected_sets and update the intersection
  for (i in 2:length(selected_sets)) {
    intersection_set <- intersect(intersection_set, selected_sets[[i]])
  }
  
  # Return the intersection set
  return(intersection_set)
}

evaluate_sets <- function(accepted_sets, possible_sets) {
  # Initialize matrix with FALSE values
  result_matrix <- matrix(FALSE, nrow = length(possible_sets), ncol = length(accepted_sets))
  colnames(result_matrix) <- paste("AccSet",1:length(accepted_sets),sep = "_")
  max_set <- numeric(length(possible_sets))
  max_occ <- numeric(length(possible_sets))
  
  # Iterate over each possible set i and each accepted set j
  for (i in 1:length(possible_sets)) {
    for (j in 1:length(accepted_sets)) {
      # Set cell value to FALSE for each combination of i and j
      result_matrix[i, j] <- FALSE
      
      # Create all subsets of possible set i
      subsets_i <- lapply(1:length(possible_sets[[i]]), function(m) combn(possible_sets[[i]], m, simplify = FALSE))
      
      # Check if each subset of i is a subset of accepted set j
      for (s in unlist(subsets_i, recursive = FALSE)) {
        if (all(s %in% accepted_sets[[j]])) {
          if (max_set[i] < length(s)) {
            max_set[i] <- length(s)
            max_occ[i] <- 1
          } else if (length(s)==max_set[i]){
            max_occ[i] <- max_occ[i] + 1
          }
          
          # If subset s of i is also a subset of accepted set j, set cell value to TRUE
          result_matrix[i, j] <- TRUE
        }
      }
    }
  }
  
  # Create a new column "number_sets" which takes the row sum across all TRUE answers
  result_matrix <- cbind(result_matrix, number_sets = rowSums(result_matrix))
  
  # Create a new column "number_vars" which counts the length of the vector of the possible set in the corresponding row
  result_matrix <- cbind(result_matrix, number_vars = sapply(possible_sets, length), 
                         max_set = max_set, max_occ = max_occ)
  
  return(result_matrix)
}

get_defining <- function(accepted_sets, S_sets){
  eval_mat <- evaluate_sets(accepted_sets, S_sets)
  L <- length(accepted_sets)
  cond1 <- eval_mat[,"number_sets"]==L
  var_min <- min(eval_mat[cond1,"number_vars"][eval_mat[cond1,"number_vars"]!=1])
  cond2 <- eval_mat[,"number_vars"]==var_min
  occ_max <- max(eval_mat[(cond1 & cond2),"max_occ"])
  cond3 <- eval_mat[,"max_occ"]==occ_max
  set_max <- max(eval_mat[(cond1 & cond2 & cond3),"max_set"])
  cond4 <- eval_mat[,"max_set"]==set_max
  
  selected_rows <- eval_mat[(cond1 & cond2 & cond3 & cond4),]
  
  if ( is.null(dim(selected_rows)[1])==FALSE ){
    acc_sets <- list()
    for (i in 1:dim(selected_rows)[1]){
      set_row <- selected_rows[i,]
      rowidx <- which(apply(eval_mat, 1, function(row) all(row == as.vector(set_row))))
      for (j in seq_along(rowidx)) {
        acc_sets[[j]] <- S_sets[[rowidx[j]]]
      }
    }
  } else {
    set_row <- selected_rows
    rowidx <- which(apply(eval_mat, 1, function(row) all(row == as.vector(set_row))))
    acc_sets <- S_sets[[rowidx]]
  }
  return(acc_sets)
}