\name{SGL-package}
\alias{SGL-package}
\docType{package}
\title{Fit a GLM (or cox model) with a combination of lasso and group lasso regularization
}
\description{Fit a regularized generalized linear model via penalized maximum likelihood.  The model is fit for a path of values of the penalty parameter. Fits linear, logistic and Cox models.
}
\details{
\tabular{ll}{
Package: \tab SGL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-3-12\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
Only 4 functions:
\code{SGL}
\code{cvSGL}
\code{predictSGL}
\code{plot.cvSGL}
}
\author{Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani

Maintainer: Noah Simon <nsimon@stanford.edu>
}
\references{Simon, N., Friedman, J., Hastie T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://www-stat.stanford.edu/~nsimon/SGL.pdf}
}
