\name{unstandardize}
\alias{unstandardize}
\title{Unstandardize the Regression Coefficient}
\description{
  Function \code{unstandardize} changes the estimated regression coefficients from standardized scale 
  to original scale.
}
\usage{
unstandardize(object, normX)
}
\arguments{
  \item{object}{a branchLars object returned by function \code{branchLars}}
  \item{normX}{the norm value of the original predictors, which is returned by the function 
               \code{standardize} or  \code{standardizeG}}
}
\value{
  The function \code{unstandardize} returns:
  \item{coef}{the estimated coefficients for all variables in original scale}
}
\author{Li Hua Yue}
\seealso{\code{\link{branchLars}}, \code{\link{standardize}}, \code{\link{standardizeG}}}
\examples{
data(diabetes)
X <- as.matrix(diabetes[, 1:10])
y <- diabetes$Y
Xy <- standardize(X, y)              
Xs <- Xy$X
ys <- Xy$y
normx <- Xy$normX

costs <- list()
costs[[1]] <- c(0, 5, 10, 20, rep(30, 6))      
costs[[2]] <- list(1, 1, 2, 3, c(4,5), c(4,6), c(4,7),c(4,8), c(4,9), c(4,10))
names(costs[[2]]) <- colnames(X)

r10 <- branchLars(Xs, ys, lambda=10, costs)
unstandardize(r10, normx)

XyG <- standardizeG(Xs, normx, costs)              
XsG <- XyG$Xs[, -12]    # exclude the squared SEX term
normX <- XyG$normX[-12]
costsG <- XyG$cost
costsG[[2]] <- costsG[[2]][-12]

r60 <- branchLars(XsG, ys, lambda=60, costsG, gamma=1)
unstandardize(r60, normX)
}
\keyword{manip}

