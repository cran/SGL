\name{SGL}
\alias{SGL}

\title{Fit a GLM with a combination of lasso and group lasso regularization}
\description{Fit a regularized generalized linear model via penalized maximum likelihood.  The model is fit for a path of values of the penalty parameter. Fits linear, logistic and Cox models.}
\usage{
SGL(data, index, type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.1, nlam = 20, gamma = 0.8, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL)
}

\arguments{
  \item{data}{For \code{type="linear"} should be a list with $x$ an input matrix of dimension n-obs by p-vars, and $y$ a length $n$ response vector. For \code{type="logit"} should be a list with $x$, an input matrix, as before, and $y$ a length $n$ binary response vector. For \code{type="cox"} should be a list with x as before, \code{time}, an n-vector corresponding to failure/censor times, and \code{status}, an n-vector indicating failure (1) or censoring (0). 
}
  \item{index}{A p-vector indicating group membership of each covariate}
  \item{type}{model type: one of ("linear","logit", "cox")}
  \item{maxit}{Maximum number of iterations to convergence}
  \item{thresh}{Convergence threshold for change in beta}
  \item{min.frac}{The minimum value of the penalty parameter, as a fraction of the maximum value}
  \item{nlam}{Number of lambda to use in the regularization path}
  \item{gamma}{Fitting parameter used for tuning backtracking (between 0 and 1)}
  \item{standardize}{Logical flag for variable standardization prior to fitting the model.}
  \item{verbose}{Logical flag for whether or not step number will be output}
 \item{step}{Fitting parameter used for inital backtracking step size (between 0 and 1)}
 \item{reset}{Fitting parameter used for taking advantage of local strong convexity in nesterov momentum (number of iterations before momentum term is reset)}
\item{alpha}{The mixing parameter. \code{alpha} = 1 is the lasso penalty. \code{alpha} = 0 is the group lasso penalty.}
\item{lambdas}{A user specified sequence of lambda values for fitting. We recommend leaving this NULL and letting SGL self-select values}
}
\details{
The sequence of models along the regularization path is fit by accelerated generalized gradient descent.
}
\value{An object with S3 class \code{"SGL"}
\item{beta}{A p by \code{nlam} matrix, giving the penalized MLEs for the nlam different models, where the index corresponds to the penalty parameter \code{lambda}}
\item{lambdas}{The actual sequence of \code{lambda} values used (penalty parameter)}
\item{type}{Response type (linear/logic/cox)}
\item{intercept}{For some model types, an intercept is fit}
\item{X.transform}{A list used in \code{predict} which gives the empirical mean and variance of the x matrix used to build the model}
\item{lambdas}{A user specified sequence of lambda values for fitting. We recommend leaving this NULL and letting SGL self-select values}
}
\references{Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://www-stat.stanford.edu/~nsimon/SGL.pdf}
}
\author{Noah Simon, Jerry Friedman, Trevor Hastie, and Rob Tibshirani\cr
Maintainer: Noah Simon \email{nsimon@stanford.edu}
}

\seealso{\code{cv.SGL}}
\examples{
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
data = list(x = X, y = y)
fit = SGL(data, index, type = "linear")
}
\keyword{model}
\keyword{regression}
