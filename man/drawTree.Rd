\name{drawTree}
\alias{drawTree}
\title{Draw the Search Tree of Branching Lars Process}
\description{
  Draw the tree showing the searching process. The red colored path is the optimal path.
}
\usage{
drawTree(object)
}
\arguments{
  \item{object}{the branchLars object returned by the function \code{branchLars}}
}
\value{
  A plot of the search tree is returned.
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

r300 <- branchLars(Xs, ys, lambda=300, costs)
drawTree(r300)
r80 <- branchLars(Xs, ys, lambda=80, costs)
drawTree(r80)
}
\keyword{tree}

