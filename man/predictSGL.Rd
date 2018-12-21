\name{predictSGL}
\alias{predictSGL}
\title{Outputs Predicted Responses from an SGL Model for New Observations}
\description{
Outputs predicted response values for new user input observations at a specified \code{lambda} value}
\usage{
predictSGL(x, newX, lam)
}
\arguments{
  \item{x}{fitted \code{"SGL"} object}
  \item{newX}{covariate matrix for new observations whose responses we wish to predict}
  \item{lam}{the index of the lambda value for the model with which we desire to predict}
}
\details{Predicted outcomes are given}
\references{Simon, N., Friedman, J., Hastie T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://faculty.washington.edu/nrsimon/SGLpaper.pdf}
}
\author{Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani\cr
Maintainer: Noah Simon <nrsimon@uw.edu>
}
\seealso{\code{SGL} and \code{cvSGL}.}
\examples{
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
data = list(x = X, y = y)
Fit = SGL(data, index, type = "linear")
X.new = matrix(rnorm(n * p), ncol = p, nrow = n)
predictSGL(Fit, X.new, 5)
}
\keyword{models}
\keyword{regression}
