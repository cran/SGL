\name{plot.cv.SGL}
\alias{plot.cv.SGL}
\title{plots the cross-validated error curve produced by cv.SGL}
\description{
Plots the cross-validated error curve, and confidence bounds for each \code{lambda} in our regularization path.}
\usage{
\method{plot}{cv.SGL}(x, ...)
}
\arguments{
  \item{x}{fitted \code{"cv.SGL"} object}
  \item{...}{additional arguments to be passed to plot}
}
\details{
A cross validated deviance plot is produced. More regularized models are to the right (less regularized to the left)}
\references{Simon, N., Friedman, J., Hastie T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://faculty.washington.edu/nrsimon/SGLpaper.pdf}
}
\author{Noah Simon, Jerome Friedman, Trevor Hastie, and Rob Tibshirani\cr
Maintainer: Noah Simon <nrsimon@uw.edu>
}
\seealso{\code{SGL} and \code{cv.SGL}.}
\examples{
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
data = list(x = X, y = y)
cvFit = cvSGL(data, index, type = "linear")
plot(cvFit)
}
\keyword{models}
\keyword{regression}
