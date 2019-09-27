SGL <- function(data, index, type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.1, nlam = 20, gamma = 0.8, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL){
  
  if(!(type %in% c("linear","logit","cox"))){
    print("'type' must be one of 'linear', 'logit', or 'cox'!")
    return(NA)
  }
  
  ## centering (and potentially scaling)
  out <- center_scale(data$x, standardize)
  data$x <- out$x
  X.transform <- out$X.transform

  if(type == "linear"){
    intercept <- mean(data$y)
    data$y <- data$y - intercept
    Sol <- oneDim(data, index, thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, reset = reset, alpha = alpha)
    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, intercept = intercept, X.transform = X.transform)
  }

  if(type == "logit"){
    Sol <- oneDimLogit(data, index, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset)
    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, intercept = Sol$intercept, X.transform = X.transform)

  }

  if(type == "cox"){
    Sol <- oneDimCox(data, index, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset) 
    Sol = list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, X.transform = X.transform)
  }
  class(Sol) = "SGL"
    return(Sol)
}


cvSGL <- function(data, index = rep(1,ncol(data$x)), type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.05, nlam = 20, gamma = 0.8, nfold = 10, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL, foldid = NULL){

  if(!(type %in% c("linear","logit","cox"))){
    print("'type' must be one of 'linear', 'logit', or 'cox'!")
    return(NA)
  }
  
  ## centering (and potentially scaling)
  out <- center_scale(data$x, standardize)
  data$x <- out$x
  X.transform <- out$X.transform
  
  ## Forming cv folds
  if(is.null(foldid)){
    foldid = form_folds(nrow(data$x), nfold)
  }
  nfold <- length(unique(foldid))

  if(type == "linear"){
    intercept <- mean(data$y)
    data$y <- data$y - intercept
    Sol <- linCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, reset = reset, alpha = alpha, foldid = foldid)
    Sol$fit = list(beta = Sol$fit$beta, lambdas = Sol$fit$lambdas, intercept = intercept, step = step)
    Sol$prevals <- Sol$prevals + intercept
    }
  if(type == "logit"){
    Sol <- logitCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset, foldid = foldid)
  }
  if(type == "cox"){
    Sol <- coxCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset, foldid = foldid)

  }

  Sol = list(fit = Sol$fit, lldiff = Sol$lldiff, lambdas = Sol$lambdas, type = type, llSD = Sol$llSD, prevals = Sol$prevals, foldid = foldid)

class(Sol) = "cv.SGL"

    return(Sol)
}


## This function always centers X (and scales if standardize == TRUE)
## It also returns the centering/scaling numbers (for use in prediction with new data)
center_scale <- function(X, standardize){
  means <- apply(X,2,mean)
  X <- t(t(X) - means)
  
  X.transform <- list(X.means = means)
  
  if(standardize == TRUE){
    var <- apply(X,2,function(x)(sqrt(sum(x^2))))
    X <- t(t(X) / var)
    X.transform$X.scale <- var
  }
  else{
    X.transform$X.scale <- 1
  }
  
  return(list(x = X, X.transform = X.transform))
}

## This function forms the CV folds

form_folds <- function(n_obs, nfold){
  folds <- cut(seq(1,n_obs),breaks=nfold,labels=FALSE)
  folds_final <- sample(folds, replace = FALSE)
  return(folds_final)
}
