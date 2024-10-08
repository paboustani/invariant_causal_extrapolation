#' Engression Fit Function
#'
#' This function fits an Engression model to the provided data. It allows for the tuning of 
#' several parameters related to model complexity and training. The function is not meant to 
#' be exported but can be used within the package or for internal testing purposes.
#'
#' @param X A matrix or data frame representing the predictors.
#' @param Y A matrix representing the target variable(s).
#' @param noise_dim The dimension of the noise introduced in the model (default: 100).
#' @param hidden_dim The size of the hidden layer in the model (default: 100).
#' @param num_layer The number of layers in the model (default: 3).
#' @param dropout The dropout rate to be used in the model in case no batch normalization is used (default: 0.01)
#' @param batch_norm A boolean indicating whether to use batch-normalization (default: TRUE).
#' @param num_epochs The number of epochs to be used in training (default: 200).
#' @param lr The learning rate to be used in training (default: 10^-3).
#' @param beta The beta scaling factor for energy loss (default: 1).
#' @param silent A boolean indicating whether to suppress output during model training (default: FALSE).
#'
#' @return A list containing the trained engression model and a vector of loss values.
#'
#' @keywords internal
#' 
gaussianengressionfit <- function(X,Y, noise_dim=0, hidden_dim=100, num_layer=3, dropout=0.01,batch_norm=TRUE, num_epochs=200,lr=10^(-3), beta=1,  silent=FALSE){
    in_dim = dim(X)[2]
    out_dim = dim(Y)[2]
    if(num_layer<=2){
        if(!batch_norm){
            model = nn_sequential( nn_linear(in_dim+noise_dim,hidden_dim),nn_dropout(dropout), nn_elu(), nn_linear(hidden_dim,out_dim))
        }else{
            model = nn_sequential( nn_linear(in_dim+noise_dim,hidden_dim), nn_elu(),nn_batch_norm1d(hidden_dim), nn_linear(hidden_dim,out_dim))

        }
    }else{
        if(!batch_norm){
            hid =  nn_sequential(nn_linear(hidden_dim, hidden_dim),nn_elu())
            if(num_layer>3) for (lay in 3:num_layer) hid = nn_sequential(hid,nn_sequential(nn_linear(hidden_dim, hidden_dim),nn_elu()) )
            model = nn_sequential( nn_sequential(nn_linear(in_dim+noise_dim,hidden_dim),nn_dropout(dropout), nn_elu()),hid, nn_linear(hidden_dim,out_dim))  
        }else{
            hid =  nn_sequential(nn_linear(hidden_dim, hidden_dim),nn_elu(),nn_batch_norm1d(hidden_dim))
            if(num_layer>3) for (lay in 3:num_layer) hid = nn_sequential(hid,nn_sequential(nn_linear(hidden_dim, hidden_dim),nn_elu(),nn_batch_norm1d(hidden_dim)) )
            model = nn_sequential( nn_sequential(nn_linear(in_dim+noise_dim,hidden_dim), nn_elu(),nn_batch_norm1d(hidden_dim)),hid, nn_linear(hidden_dim,out_dim)) 
        }
    }
    model$train()

    optimizer = optim_adam(model$parameters,lr=lr)

    n= dim(X)[1]
    lossvec = matrix(nrow=num_epochs, ncol=3)
    colnames(lossvec) = c("energy-loss","E(|Y-Yhat|)","E(|Yhat-Yhat'|)")
    printat = pmax(1,floor((seq(1,num_epochs, length=11))))

    for (iter in 1:num_epochs){
        optimizer$zero_grad()
        if(noise_dim>0){
            stop("Engression can only be implemented with gaussian noise if noise_dim is set to 0.")
        }else{
            xt = torch_tensor( X + torch_randn_like(X) , dtype=torch_float(),requires_grad=TRUE) # 
            xpt = torch_tensor( X + torch_randn_like(X) , dtype=torch_float(),requires_grad=TRUE) # 
            yt = torch_tensor(Y, dtype=torch_float(),requires_grad=TRUE)
        }
        la = energylossall(yt,model(xt),model(xpt))
        lossvec[iter, ] = signif(c(sapply(la, as.numeric)),3 )
        if(beta==1) loss = energyloss(yt,model(xt),model(xpt)) else loss= energylossbeta(yt,model(xt),model(xpt),beta)
        loss$backward()
        optimizer$step()
        if(!silent){
            cat("\r ", round(100*iter/num_epochs), "% complete, epoch: ", iter)
            if(iter %in% printat){cat("\n");  print(lossvec[iter,])}
        }
    } 
    if(batch_norm) model$train(mode=FALSE)

    if(noise_dim>0){
      engressor = function(x) as.matrix(model( torch_tensor( x + torch_randn_like(x), dtype=torch_float())),ncol=out_dim) #  
    }else{
      engressor = function(x) as.matrix(model( torch_tensor( x + torch_randn_like(x), dtype=torch_float())),ncol=out_dim) #  + torch_randn_like(x)
    }
    return(list(engressor=engressor, lossvec=lossvec))
}
