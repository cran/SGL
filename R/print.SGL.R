print.SGL = function(x, digits = max(3, getOption("digits") - 3), ...){
  num.nonzero <- apply(x$beta,2, function(z){sum(z != 0)})
  cat("\n regression type: ", x$type, "\n\n")
  print(cbind(lambdas = x$lambdas, num.nonzero = num.nonzero))
}
