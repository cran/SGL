\name{print.SGL}
\alias{print.SGL}
\title{prints a summary of the SGL solution path}
\description{
Prints a short summary of the SGL solution path.}
\usage{
\method{print}{SGL}(x, digits, ...)
}
\arguments{
  \item{x}{fitted \code{"SGL"} object}
  \item{digits}{significant digits in printout}
  \item{...}{additional print arguments}
}
\details{
The time of regression run, followed by a 2-column matrix with rows \code{ lambdas} and \code{ num.nonzero}. \code{ lambdas} gives the lambda-value of each fit. \code{ num.nonzero} gives the the number of non-zero coefficients.}
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
fit = SGL(data, index, type = "linear")
print(fit)
}
\keyword{models}
\keyword{regression}
