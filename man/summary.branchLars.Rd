\name{summary.branchLars}
\alias{summary.branchLars}
\title{Summary of the Optimal Total Loss}
\description{
  \code{summary.branchLars} shows the optimal total loss and its three components: RSS, L1 penalty, and total cost 
  of collecting those predictors selected by branchLars.
}
\usage{
\method{summary}{branchLars}(object, \dots)
}
\arguments{
  \item{object}{the branchLars object returned by the function \code{branchLars}}
  \item{\dots}{additional arguments for summary generic}
}
\value{
  It returns the optimal total loss and its three components.
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
costs[[1]] <- c(0, 5, 10, 20, rep(30, 6))      
costs[[2]] <- list(1, 1, 2, 3, c(4,5), c(4,6), c(4,7),c(4,8), c(4,9), c(4,10))
names(costs[[2]]) <- colnames(X)
r10 <- branchLars(Xs, ys, lambda=10, costs)
summary(r10)
unstandardize(r10, Xy$normX)
}
\keyword{regression}

