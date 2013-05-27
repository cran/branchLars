\name{branchLars}
\alias{branchLars}
\title{Select Cost-efficient Variables Using Branching LARS Method}
\description{
  For a specific lasso regularization paramter \code{lambda}, the functions select cost-efficient variables 
  by Branching LARS search method.
}
\usage{
branchLars(X, y, lambda, cost, gamma = 1, sumC=0)
}
\arguments{
  \item{X}{the standardized matrix of predictors}
  \item{y}{the standardized response vector}
  \item{lambda}{lasso regularization parameter}
  \item{cost}{the cost of collecting the potential predictors.}
  \item{gamma}{The weight on cost. It reflects the reluctance of the user to select high cost predictors.
   The larger the value, the more reluctant the user will be. It can be any real value greater than or 
   equal to 0. Default is 1.}
  \item{sumC}{the sum of the costs of the previously selected variables using the same \code{lambda} and 
   a smaller \code{gamma}}
}
\value{
  A "branchLars" object is returned, which could be used in \code{print}, \code{summary}, \code{predict}, 
  \code{unstandardize}, \code{drawTree}, \code{bic}, \code{cp}, and \code{lambdaOpt} functions. 
}
\author{Li Hua Yue}
\seealso{\code{\link{buildCostFun}}, \code{\link{standardize}}, \code{\link{standardizeG}}, 
         \code{\link{unstandardize}}, \code{\link{print.branchLars}},
         \code{\link{summary.branchLars}}, \code{\link{predict.branchLars}}, 
         \code{\link{drawTree}}, \code{\link{bic}}, \code{\link{cp}},
         \code{\link{lambdaOpt}}}
\examples{
data(diabetes)
X <- as.matrix(diabetes[, 1:10])
y <- diabetes$Y
Xy <- standardize(X, y)              
Xs <- Xy$X
ys <- Xy$y
costs <- list()
costs[[1]] <- c(0, 5, 10, 20, rep(30, 6))      # costs (including group cost)
costs[[2]] <- list(1, 1, 2, 3, c(4,5), c(4,6), # cost pointers for each variable
              c(4,7),c(4,8), c(4,9), c(4,10))
names(costs[[2]]) <- colnames(X)

result1 <- branchLars(Xs, ys, lambda=20, costs, gamma=0.3)
unstandardize(result1, Xy$normX)
result1$evals
drawTree(result1)

result2 <- branchLars(Xs, ys, lambda=20, costs, gamma=1, sumC=result1$sumC)
unstandardize(result2, Xy$normX)
result2$evals
result3 <- branchLars(Xs, ys, lambda=20, costs, gamma=7, sumC=result2$sumC)
unstandardize(result3, Xy$normX)
result3$evals
}
\keyword{regression}

