\name{lambdaOpt}
\alias{lambdaOpt}
\title{Choose Optimal Tuning Parameter Lambda}
\description{
  Given a user defined \code{gamma} value, the function \code{lambdaOpt} find the optimal tuning 
  parameter \code{lambda} which will be used to select an optimal Branching LARS model.
}
\usage{
lambdaOpt(X, y, cost, gamma=1, lower=0, upper=Inf, method="BIC", 
          tol=.Machine$double.eps^0.25/3)
}
\arguments{
  \item{X}{matrix of predictors}
  \item{y}{response vector}
  \item{cost}{the cost of collecting the potential predictors.}
  \item{gamma}{The weight on cost. It reflects the reluctance of the user to select high cost predictors.
   The larger the value, the more reluctant the user will be. It can be any real value greater than or 
   equal to 0. Default is 1.}
  \item{lower}{the lower limit of the lambda range used to search for the optimum. Default is 0.} 
  \item{upper}{the upper limit of the lambda range used to search for the optimum. Default is infinity.}
  \item{method}{the model selection criterion. Can choose between BIC and Cp. Default is BIC. }
  \item{tol}{the desired length of the interval of uncertainty of the final result}
}
\details{
  The golden section search is implemented by this function. If the objective funtion of Branching LARS for a 
  particular \code{gamma} is not unimodal, the function \code{lambdaOpt} may return a local, but perhaps non-global,
  minimun BIC (or Cp) value and the corresponding local optimal \code{lambda} value. Choosing proper lower and upper limits 
  is essential in this case, and deviding the search range to several smaller parts may help.
}
\value{
  \code{lambdaOpt} returns:
  \item{Optimal_Lambda}{the optimal value of tuning parameter \code{lambda}}
  \item{BIC (or Cp)}{the corresponding minimun BIC (or Cp) value}
  \item{sumC}{the sum of the costs of the variables selected using the optimal \code{lambda}}
}
\author{Li Hua Yue}
\seealso{\code{\link{branchLars}}, \code{\link{bic}}, \code{\link{cp}}}
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

opt <- lambdaOpt(Xs, ys, costs, gamma=1)
modelOpt <- branchLars(Xs, ys, lambda=opt$Optimal_Lambda, costs, gamma=1, sumC=opt$sumC)
unstandardize(modelOpt, Xy$normX)
}
\keyword{optimize}

