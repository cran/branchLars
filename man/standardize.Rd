\name{standardize}
\alias{standardize}
\title{Standardize the Data Set}
\description{
  Function \code{standardize} standardizes the predictors to have mean 0 and norm 1, and standardizes the 
  response vector to have mean 0. 
}
\usage{
standardize(X, y)
}
\arguments{
  \item{X}{matrix of predictors}
  \item{y}{response vector}
}
\value{
  The function \code{standardize} returns:
  \item{y}{standardized response vector}
  \item{X}{standardized matrix of predictors}
  \item{normX}{a vector of norm values of the original predictors}
}
\author{Li Hua Yue}
\seealso{\code{\link{branchLars}}}
\examples{
data(diabetes)
X <- as.matrix(diabetes[, 1:10])
y <- diabetes$Y
Xy <- standardize(X, y)              
Xs <- Xy$X
ys <- Xy$y
normX <- Xy$normX

costs <- list()
costs[[1]] <- c(0, 5, 10, 20, rep(30, 6))      
costs[[2]] <- list(1, 1, 2, 3, c(4,5), c(4,6), c(4,7),c(4,8), c(4,9), c(4,10))
names(costs[[2]]) <- colnames(X)
r10 <- branchLars(Xs, ys, lambda=10, costs)
unstandardize(r10, normX)
}
\keyword{manip}

