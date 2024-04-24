# Define the DAG

# number of nodes
d <- 6

# Define the DAG
dag <- list(
  Z_1 = character(0),
  Z_2 = "Z_1",
  Z_3 = c("Z_1", "Z_2"),
  Z_4 = "Z_3",
  Z_5 = character(0),
  Z_6 = c("Z_3", "Z_4", "Z_5")
)

# coeff_signs list
dag_signs <- list(
  Z_1 = "NA",
  Z_2 = "+",
  Z_3 = c("+", "-"),
  Z_4 = "-",
  Z_5 = "NA",
  Z_6 = c("+", "-", "+")
)

# Obtain all possible combinations of regressors for a target of choice
obtain_S <- function(target, d){
  # Define the set of nodes
  nodes <- paste("Z", 1:d, sep = "_")
  
  # Exclude the target variable from the set of nodes
  remainings <- nodes[nodes != target]
  
  # Initialize an empty list to store combinations
  all_combos <- list()
  
  # Generate all possible combinations of the remaining nodes
  for (i in 1:length(remainings)) {
    combos <- combn(remainings, i)
    for (j in 1:ncol(combos)) {
      all_combos[[length(all_combos) + 1]] <- as.vector(combos[, j])
    }
  }
  
  return(all_combos)
}

# Functions to print parents, children, ancestors and descendants of a node in a DAG

# Function to return parents of a node in a DAG
get_parents <- function(node, dag) {
  parents <- dag[[node]]  # Get parents of the node
  return(parents)
}

# Function to find ancestors of a node in a DAG
find_ancestors <- function(node, dag) {
  ancestors <- character(0)  # Initialize an empty vector to store ancestors
  
  # Helper function to recursively find ancestors
  find_ancestors_recursive <- function(node, dag, ancestors) {
    # Find parents of the current node
    parents <- dag[[node]]
    
    # If the node has no parents, return the current ancestors
    if (is.null(parents) || length(parents) == 0) {
      return(ancestors)
    }
    
    # Iterate through parents and add them to ancestors
    for (parent in parents) {
      if (!(parent %in% ancestors)) {
        ancestors <- c(ancestors, parent)
        # Recursively find ancestors of the parent
        ancestors <- find_ancestors_recursive(parent, dag, ancestors)
      }
    }
    return(ancestors)
  }
  
  # Call the recursive function to find ancestors
  ancestors <- find_ancestors_recursive(node, dag, ancestors)
  
  return(ancestors)
}

# Function to find descendants of a node in a DAG
find_descendants <- function(node, dag) {
  descendants <- character(0)  # Initialize an empty vector to store descendants
  
  # Helper function to recursively find descendants
  find_descendants_recursive <- function(node, dag, descendants) {
    # Find children of the current node
    children <- names(dag)[sapply(dag, function(x) node %in% x)]
    
    # If the node has no children, return the current descendants
    if (length(children) == 0) {
      return(descendants)
    }
    
    # Iterate through children and add them to descendants
    for (child in children) {
      if (!(child %in% descendants)) {
        descendants <- c(descendants, child)
        # Recursively find descendants of the child
        descendants <- find_descendants_recursive(child, dag, descendants)
      }
    }
    return(descendants)
  }
  
  # Call the recursive function to find descendants
  descendants <- find_descendants_recursive(node, dag, descendants)
  
  return(descendants)
}

# Function to find children of a node in a DAG
find_children <- function(node, dag) {
  children <- character(0)  # Initialize an empty vector to store children
  
  # Find children of the node
  children <- names(dag)[sapply(dag, function(x) node %in% x)]
  
  return(children)
}

# Function to print the coefficients of the parents of a node in a DAG
get_coeff <- function(node, coeffs){
  return(coeffs[[node]])
}

# Multiplicative or additive effects

# Function to simulate data from the DAG
g_k <- function(parents, coeff_signs, id, multiplic) {
  # Convert coefficient signs to numeric values (+1 for "+" and -1 for "-")
  coeff_values <- ifelse(coeff_signs == "+", 1, -1)
  
  # apply coefficients, function f_id and multiplication
  if (length(coeff_values) == 1){
    parents <- parents * coeff_values
    data_value <- f(parents, id)
  } else {
    parents <- parents %*% diag(coeff_values)
    parents <- matrix(data = f(parents, id), nrow = dim(parents)[1], ncol = dim(parents)[2])
    
    if (multiplic == TRUE) {
      # Calculate product of parental variables
      data_value <- apply(parents, 1, prod)
    } else {
      # Calculate sum of parental variables
      data_value <- rowSums(parents)
    }
  }
  
  return(data_value)
}

# Shift- or do-Interventions
interviene <- function(target, shift, epsilon){
  if(shift==TRUE){
    # shift-interventions
    target <- target + epsilon
  } else {
    # do-interventions
    target <- epsilon
  }
  return(target)
}

# Tail behavior of the noise
gen_eta <- function(n, eta_df){
  noise <- rt(n = n, df = eta_df)
  return(noise)
}

# Strength of interventions

gen_epsilon <- function(n, epsilon_df, meanshift, strength){
  epsilon <- rt(n = n, df = epsilon_df)
  epsilon <- strength*(epsilon + meanshift)
  return(epsilon)
}

# Non-linearities
f <- function(x, id) {
  if (id == 1) {
    return( x )
    
  } else if (id == 2) {
    return( pmax(0, x) )
    
  } else if (id == 3) {
    return( sign(x) * sqrt(abs(x)) )
    
  } else if (id == 4) {
    return( sin(2 * pi * x) )
    
  } else if (id == "softplus"){
    return( log(1+exp(x)) )
    
  } else if (id == "square"){
    return( x^2/2 )
    
  } else if (id == "cubic"){
    return( x^3/3 )
    
  } else if (id == "log"){
    return( if (x<=2) { (x-2)/3 + log(3) } else { log(x) } )
  
  } else {
    stop("Specify valid function class.")
      
  }
}

# Location of interventions
int_loc <- function(target,dag,interv){
  
  location <- character(0)
  
  if (interv=="all"){
    # all ancestors and descendants
    location <- names(dag)[names(dag) != target]
    
  } else if (interv=="rand"){
    # one of the ancestors at random
    ancestors <- c(find_ancestors(node = target,dag = dag))
    
    if (length(ancestors) == 0){
      location <- character(0)
    } else{
      location <- sample(ancestors,1)
    }
    
    # on of the descendants at random
    descendants <- c(find_descendants(node = target,dag = dag))
    
    if (length(descendants) == 0){
      location <- append(location, character(0))
    } else{
      location <- append(location, sample(descendants,1))
    }
    
  } else {
    # one of the parents at random
    parents <- c(get_parents(node = target,dag = dag))
    
    if (length(parents) == 0){
      location <- character(0)
    } else{
      location <- sample(parents,1)
    }
    
    # on of the children at random
    children <- c(find_children(node = target,dag = dag))
    
    if (length(children) == 0){
      location <- append(location, character(0))
    } else{
      location <- append(location, sample(children,1))
    }
  }
  
  return(location)
}

# simulate data (with and without intevention)
sim_data <- function(n, d, dag, coeff_signs, target, eta_df, f_id, multiplic, 
                     E, shift, interv, epsilon_df, meanshift, strength, verbose = FALSE) {
  
  if (verbose) cat("Entironment is set to E =", E, "\n")
  doc <- c(n,E)
  
  X <- matrix(NA, nrow = n, ncol = d+1)
  X[,d+1] <- rep(E,n)
  colnames(X) <- c(paste("Z", 1:d, sep = "_"),"E")
  if (verbose) print(X)
  
  if (E == 1) {
    doc <- append(doc, rep(FALSE, d))
    if (verbose) cat("Simulating without Interventions. \n")
    
    for (k in 1:d) {
      node <- paste("Z", k, sep = "_")
      if (verbose) cat("\n Starting simulating node =", node, "\n")
      parents <- get_parents(node, dag)
      if (verbose) cat("Parents =", parents, "\n")
      
      if(length(parents) == 0) {
        X[,k] <- gen_eta(n, eta_df)
        if (verbose) cat("No parents, generated eta. \n")
        if (verbose) print(X)
        
      } else {
        signs <- get_coeff(node = node, coeffs = coeff_signs)
        if (verbose) cat("Parental signs =", signs, "\n")
        
        X[,k] <- g_k(parents = X[,parents], coeff_signs = signs,
                     id = f_id, multiplic = multiplic) + gen_eta(n, eta_df)
        if (verbose) cat("There are parents and added eta. \n")
        if (verbose) print(X)
        
      }
    }
  } else {
    if (verbose) cat("Simulating with Interventions. \n")
    
    location <- int_loc(target = target, dag = dag, interv = interv)
    if (verbose) cat("Interviening on location(s) =", location, "\n")
    
    for (k in 1:d) {
      doc <- append(doc, FALSE)
      node <- paste("Z", k, sep = "_")
      if (verbose) cat("\n Simulating node =", node, "\n")
      parents <- get_parents(node, dag)
      if (verbose) cat("parents =", parents, "\n")
      
      if (length(parents) == 0) {
        X[,k] <- gen_eta(n, eta_df)
        if (verbose) cat("No parents; generated eta. \n")
        if (verbose) print(X)
        
        if (node %in% location){
          doc[2+k] <- TRUE
          epsilon <- gen_epsilon(n, epsilon_df, meanshift, strength)
          if (verbose) print(epsilon)
          X[,k] <- interviene(target = X[,k], shift = shift, epsilon = epsilon)
          if (verbose) cat("No parents and node is in interventions; added interviened via epsilon. \n")
          if (verbose) print(X)
          
        }
      } else {
        signs <- get_coeff(node = node, coeffs = coeff_signs)
        if (verbose) cat("parental signs =", signs, "\n")
        
        X[,k] <- g_k(parents = X[,parents], 
                     coeff_signs = signs,
                     id = f_id, multiplic = multiplic) + gen_eta(n, eta_df)
        if (verbose) cat("There are parents; generated transform plus eta. \n")
        if (verbose) print(X)
        
        if (node %in% location){
          doc[2+k] <- TRUE
          epsilon <- gen_epsilon(n, epsilon_df, meanshift, strength)
          if (verbose) print(epsilon)
          X[,k] <- interviene(target = X[,k], shift = shift, epsilon = epsilon)
          if (verbose) cat("Node is in interventions; interviened via epsilon. \n")
          if (verbose) print(X)
          
        }
      }
    }
  }
  return(list(sample = X, documentation = doc))
}

# simulate data
gen_sample <- function(n, d, dag, signs, target = "random", verbose = FALSE){
  
  # Target variable
  if (!target %in% paste("Z", 1:d, sep = "_")) {
    if (target != "random") {
      stop("Specify valid target.")
    }
    target <- paste("Z", sample(1:d, 1), sep = "_")
  } 
  
  # Tail behavior of the noise
  eta_df <- sample(c(2, 3, 5, 10, 20, 50, 100), 1)
  
  # Multiplicative or additive effects
  multiplic <- sample(c("product", "sum"), 1)
  
  # Shift- or do-Interventions
  shift <- sample(c(TRUE, FALSE),1)
  
  # Strength of interventions
  epsilon_df <- sample(c(2, 3, 5, 10, 20, 50, 100), 1)
  meanshift <- sample(c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10), 1)
  strength <- sample(c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10), 1)
  
  # Non-linearities
  id <- sample(c(1, 2, 3, 4), 1)
  
  # Location of interventions
  n_vec <- sample(1:n, 1)
  n_vec <- append(n_vec, sample(1:(n-n_vec), 1))
  n_vec <- append(n_vec, n - sum(n_vec))
  
  output <- list(sample = numeric(0), documentation = character(0))
  
  for (E in 1:3){
    interv <- sample(c("all", "rand", "close"), 1)
    
    # sample from dag
    out <- sim_data(n_vec[E], d, dag, coeff_signs = signs, target, eta_df, f_id = id, 
                    multiplic, E, shift, interv, epsilon_df, meanshift, strength, 
                    verbose = verbose)
    
    output$sample <- rbind(output$sample, out$sample)
    if (E==1) {
      output$documentation <- rbind(output$documentation,c(out$documentation,"None"))
    } else {
      output$documentation <- rbind(output$documentation,c(out$documentation,interv))
    }
  }
  
  # generate ID
  ID <- paste(sample(letters, 1), sample(10^5:10^6, 1), sep = "-")
  
  # output$sample <- cbind(output$sample, rep(ID,n))
  output$documentation <- cbind(rep(ID,3), output$documentation)
  
  # create documentation 
  stats <- c(ID, target, eta_df, multiplic, shift, epsilon_df, meanshift, strength, id, n)
  names(stats) <- c("ID", "target", "eta_df", "multiplic", "shift", "epsilon_df", 
                    "meanshift", "strength", "F_id", "N")
  colnames(output$documentation) <- c("ID", "N", "E", "Z1_int", "Z2_int", "Z3_int", 
                                      "Z4_int", "Z5_int", "Z6_int", "interv")
  
  return(list(sample = output$sample,statistics = stats, documentation = output$documentation))
}