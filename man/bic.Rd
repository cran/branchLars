\name{bic}
\alias{bic}
\title{BIC Calculation for a BranchLars Object}
\description{
  \code{bic} calculates the BIC value for a model selected by \code{branchLars} function. 
}
\usage{
bic(object)
}
\arguments{
  \item{object}{the branchLars object returned by the function \code{branchLars}}
}
\value{
  it returns the BIC value for a model selected by branchLars function.
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
costs <- list()
costs[[1]] <- c(0, 5, 10, 20, rep(30, 6))       # costs (including group cost)
costs[[2]] <- list(1, 1, 2, 3, c(4,5), c(4,6),  # cost pointers for each variable
              c(4,7),c(4,8), c(4,9), c(4,10))
names(costs[[2]]) <- colnames(Xs)

r10 <- branchLars(Xs, ys, lambda=10, costs)
bic(r10)
}
\keyword{utilities}

