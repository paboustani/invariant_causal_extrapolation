engressionResidualDistributions <- function(
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
    standardize = TRUE
){

  n <- nrow(X)
  p <- ncol(X)

  # train model Y ~ X using all data and predict
  engressionRes <- engression(
    X,
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
  
  predicted = predict(engressionRes, X)

  # if(!returnModel){
    list(predicted = predicted)
  # }else{
  #   list(predicted = predicted, model = engressionRes$model)
  # }
}

