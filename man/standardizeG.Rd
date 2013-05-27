\name{standardizeG}
\alias{standardizeG}
\title{Add the Standardized Squared and Two-way Interaction Terms}
\description{
  Function \code{standardizeG} standardizes the squared and two-way 
  interaction terms of the predictors, and add them to the design matrix. 
  It also constructs the \code{cost} variable used in \code{branchLars} 
  for the more complicated design matrix if user inputs a simpler \code{cost} variable 
  for the smaller design matrix used in the main effect model.  
}
\usage{
standardizeG(Xs, normX, cost="NULL")
}
\arguments{
  \item{Xs}{standardized design matrix used in main effect model}
  \item{normX}{a vector of norm values of the original predictors}
  \item{cost}{the cost of collecting the data. It is the \code{cost} used in main effect model,
              see \code{buildCostFun}. Users may ignore this input if they have their own cost 
              structure for a full model.}
}
\value{
  The function \code{standardizeG} returns:
  \item{Xs}{standardized matrix of predictors including squared and two-way interaction terms}
  \item{normX}{a vector of norm values of all predictors including squared and two-way 
              interaction terms}
  \item{cost}{the \code{cost} used in the input of \code{branchLars} for a full model}
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
normx <- Xy$normX

costs <- list()
costs[[1]] <- c(0, 5, 10, 20, rep(30, 6))      
costs[[2]] <- list(1, 1, 2, 3, c(4,5), c(4,6), c(4,7),c(4,8), c(4,9), c(4,10))
names(costs[[2]]) <- colnames(X)

XyG <- standardizeG(Xs, normx, costs)              
XsG <- XyG$Xs[, -12]    # exclude the squared SEX term
normX <- XyG$normX[-12]
costsG <- XyG$cost
costsG[[2]] <- costsG[[2]][-12]

r60 <- branchLars(XsG, ys, lambda=60, costsG, gamma=1)
unstandardize(r60, normX)
}
\keyword{manip}

