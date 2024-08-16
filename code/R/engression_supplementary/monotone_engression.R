monotone_engression <- function(X,Y,  noise_dim=5, hidden_dim=100, num_layer=3, dropout=0.05, batch_norm=TRUE, num_epochs=1000,lr=10^(-3),beta=1, silent=FALSE, standardize=TRUE){
  
  if (is.data.frame(X)) {
    if (any(sapply(X, is.factor)))   warning("Data frame contains factor variables. Mapping to numeric values. Dummy variables would need to be created explicitly by the user.")
    X = dftomat(X)
  }
  
  if (is.vector(X) && !is.numeric(X)) X <- as.numeric(X)
  if (is.vector(X) && is.numeric(X)) X <- matrix(X, ncol = 1)
  if(is.vector(Y)) Y= matrix(Y, ncol=1)
  for (k in 1:ncol(Y)) Y[,k] = as.numeric(Y[,k])
  
  if(dropout<=0 & noise_dim==0){
    warning("dropout and noise_dim cannot both be equal to 0 as model needs to be stochastic. setting dropout to 0.5")
    dropout = 0.5
  }
  
  muX = apply(X,2,mean)
  sddX = apply(X,2,sd)
  if(any(sddX<=0)){
    warning("predictor variable(s) ", colnames(X)[which(sddX<=0)]," are constant on training data -- results might be unreliable")
    sddX = pmax(sddX, 10^(03))
  }
  muY = apply(Y,2,mean)
  sddY = apply(Y,2,sd)
  if(any(sddY<=0)){
    warning("target variable(s) ", colnames(Y)[which(sddY<=0)]," are constant on training data -- results might be unreliable")
  }
  
  if(standardize){
    X  = sweep(sweep(X,2,muX,FUN="-"),2,sddX,FUN="/")
    Y = sweep(sweep(Y,2,muY,FUN="-"),2,sddY,FUN="/")
  }
  eng = monotone_engressionfit(X,Y, noise_dim=noise_dim,hidden_dim=hidden_dim,num_layer=num_layer,dropout=dropout, batch_norm=batch_norm, num_epochs=num_epochs,lr=lr,beta=beta, silent=silent)
  engressor = list(engressor = eng$engressor, lossvec= eng$lossvec,  muX=muX, sddX=sddX,muY=muY, sddY=sddY, standardize=standardize, noise_dim=noise_dim,hidden_dim=hidden_dim,num_layer=num_layer,dropout=dropout, batch_norm=batch_norm, num_epochs=num_epochs,lr=lr)
  class(engressor) = "engression"
  return(engressor)
}

monotone_engressionfit <- function(X, Y, noise_dim=100, hidden_dim=100, num_layer=3, dropout=0.01, batch_norm=TRUE, num_epochs=200, lr=10^(-3), beta=1, silent=FALSE){
  in_dim = dim(X)[2]
  out_dim = dim(Y)[2]
  
  # Define custom Linear layer with non-negative weights
  MonotoneLinear <- nn_module(
    initialize = function(in_features, out_features) {
      self$in_features = in_features
      self$out_features = out_features
      self$weight = nn_parameter(torch_abs(torch_randn(out_features, in_features)))
      self$bias = nn_parameter(torch_zeros(out_features))
    },
    forward = function(input) {
      nnf_linear(input, self$weight, self$bias)
    }
  )
  
  if(num_layer <= 2){
    if(!batch_norm){
      model = nn_sequential(
        MonotoneLinear(in_dim + noise_dim, hidden_dim), nn_dropout(dropout), nn_relu(),
        MonotoneLinear(hidden_dim, out_dim)
      )
    } else {
      model = nn_sequential(
        MonotoneLinear(in_dim + noise_dim, hidden_dim), nn_relu(), nn_batch_norm1d(hidden_dim),
        MonotoneLinear(hidden_dim, out_dim)
      )
    }
  } else {
    if(!batch_norm){
      hid = nn_sequential(
        MonotoneLinear(hidden_dim, hidden_dim), nn_relu()
      )
      if(num_layer > 3) {
        for (lay in 3:num_layer) {
          hid = nn_sequential(hid, nn_sequential(
            MonotoneLinear(hidden_dim, hidden_dim), nn_relu()
          ))
        }
      }
      model = nn_sequential(
        nn_sequential(MonotoneLinear(in_dim + noise_dim, hidden_dim), nn_dropout(dropout), nn_relu()), hid,
        MonotoneLinear(hidden_dim, out_dim)
      )
    } else {
      hid = nn_sequential(
        MonotoneLinear(hidden_dim, hidden_dim), nn_relu(), nn_batch_norm1d(hidden_dim)
      )
      if(num_layer > 3) {
        for (lay in 3:num_layer) {
          hid = nn_sequential(hid, nn_sequential(
            MonotoneLinear(hidden_dim, hidden_dim), nn_relu(), nn_batch_norm1d(hidden_dim)
          ))
        }
      }
      model = nn_sequential(
        nn_sequential(MonotoneLinear(in_dim + noise_dim, hidden_dim), nn_relu(), nn_batch_norm1d(hidden_dim)), hid,
        MonotoneLinear(hidden_dim, out_dim)
      )
    }
  }
  model$train()
  
  optimizer = optim_adam(model$parameters, lr=lr)
  
  n = dim(X)[1]
  lossvec = matrix(nrow=num_epochs, ncol=3)
  colnames(lossvec) = c("energy-loss", "E(|Y-Yhat|)", "E(|Yhat-Yhat'|)")
  printat = pmax(1, floor((seq(1, num_epochs, length=11))))
  
  for (iter in 1:num_epochs) {
    optimizer$zero_grad()
    if(noise_dim > 0){
      xt = torch_tensor(cbind(X, matrix(rnorm(n * noise_dim), ncol=noise_dim)), dtype=torch_float(), requires_grad=TRUE)
      xpt = torch_tensor(cbind(X, matrix(rnorm(n * noise_dim), ncol=noise_dim)), dtype=torch_float(), requires_grad=TRUE)
      yt = torch_tensor(Y, dtype=torch_float(), requires_grad=TRUE)
    } else {
      xt = torch_tensor(X, dtype=torch_float(), requires_grad=TRUE)
      xpt = torch_tensor(X, dtype=torch_float(), requires_grad=TRUE)
      yt = torch_tensor(Y, dtype=torch_float(), requires_grad=TRUE)
    }
    la = energylossall(yt, model(xt), model(xpt))
    lossvec[iter, ] = signif(c(sapply(la, as.numeric)), 3)
    if(beta == 1) {
      loss = energyloss(yt, model(xt), model(xpt))
    } else {
      loss = energylossbeta(yt, model(xt), model(xpt), beta)
    }
    loss$backward()
    optimizer$step()
    if(!silent){
      cat("\r ", round(100 * iter / num_epochs), "% complete, epoch: ", iter)
      if(iter %in% printat){
        cat("\n")
        print(lossvec[iter, ])
      }
    }
  }
  if(batch_norm) model$train(mode=FALSE)
  
  if(noise_dim > 0){
    engressor = function(x) as.matrix(model(torch_tensor(cbind(x, matrix(rnorm(nrow(x) * noise_dim), ncol=noise_dim)), dtype=torch_float())), ncol=out_dim)
  } else {
    engressor = function(x) as.matrix(model(torch_tensor(x, dtype=torch_float())), ncol=out_dim)
  }
  return(list(engressor=engressor, lossvec=lossvec))
}
