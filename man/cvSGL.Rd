\name{cvSGL}
\alias{cvSGL}

\title{Fit and Cross-Validate a GLM with a Combination of Lasso and Group Lasso Regularization}
\description{Fits and cross-validates a regularized generalized linear model via penalized maximum likelihood.  The model is fit for a path of values of the penalty parameter, and a parameter value is chosen by cross-validation. Fits linear, logistic and Cox models.}
\usage{
cvSGL(data, index = rep(1, ncol(data$x)), type = "linear", maxit = 1000, thresh = 0.001,
min.frac = 0.05, nlam = 20, gamma = 0.8, nfold = 10, standardize = TRUE,
verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL,
foldid = NULL)
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
  \item{nfold}{Number of folds of the cross-validation loop}
  \item{standardize}{Logical flag for variable standardization (scaling) prior to fitting the model.}
  \item{verbose}{Logical flag for whether or not step number will be output}

 \item{step}{Fitting parameter used for inital backtracking step size (between 0 and 1)}
 \item{reset}{Fitting parameter used for taking advantage of local strong convexity in nesterov momentum (number of iterations before momentum term is reset)}

  \item{alpha}{
The mixing parameter. \code{alpha} = 1 is the lasso penalty.}

\item{lambdas}{A user inputted sequence of lambda values for fitting. We recommend leaving this NULL and letting SGL self-select values}

\item{foldid}{An optional user-pecified vector indicating the cross-validation fold in which each observation should be included. Values in this vector should range from 1 to nfold. If left unspecified, SGL will randomly assign observations to folds}
}

\details{
The function runs \code{SGL} \code{nfold}+1 times; the initial run is to find the \code{lambda} sequence, subsequent runs are used to compute the cross-validated error rate and its standard deviation.
}
\value{An object with S3 class \code{"cv.SGL"}
  \item{lldiff}{
An \code{nlam} vector of cross validated negative log likelihoods (squared error loss in the \code{linear} case, along the regularization path)
}
  \item{llSD}{
An \code{nlame} vector of approximate standard deviations of \code{lldiff}
}
  \item{lambdas}{The actual list of \code{lambda} values used in the regularization path.
}
  \item{type}{Response type (linear/logic/cox)}
  \item{fit}{A model fit object created by a call to \code{SGL} on the entire dataset
}
  \item{foldid}{A vector indicating the cross-validation folds that each observation is assigned to}
  \item{prevals}{A matrix of prevalidated predictions for each observation, for each lambda-value}
}
\references{Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2011)
  \emph{A Sparse-Group Lasso}, \cr
  \url{http://faculty.washington.edu/nrsimon/SGLpaper.pdf}
  }
\author{Noah Simon, Jerry Friedman, Trevor Hastie, and Rob Tibshirani\cr
Maintainer: Noah Simon \email{nrsimon@uw.edu}
}
\seealso{\code{SGL}
}
\examples{
set.seed(1)
n = 50; p = 100; size.groups = 10
index <- ceiling(1:p / size.groups)
X = matrix(rnorm(n * p), ncol = p, nrow = n)
beta = (-2:2)
y = X[,1:5] \%*\% beta + 0.1*rnorm(n)
data = list(x = X, y = y)
cvFit = cvSGL(data, index, type = "linear")
}
\keyword{model}
\keyword{regression}
