\name{buildCostFun}
\alias{buildCostFun}
\title{Build the Cost Calculating Function}
\description{
  This function is called by \code{branchLars} to build a cost calculating function \code{costFun(alpha)} 
  which takes \code{alpha} as input and outputs two results: the total cost of collecting data for the 
  predictors selected by the model; a vector containing the costs of each variables that have not been 
  searched in the Branching LARS search process.   
  The vector \code{alpha} has length \code{p}, where \code{p} is the number of predictors. It consists
  of -1's, 0's and 1's, where -1 means the corresponding predictor has not been searched, while 0 means 
  the predictor is not selected in the model, and 1 means that predictor is selected. The user could 
  write his(her) own function \code{buildCostFun(cost)} to incorporate more general kind of cost.
}
\usage{
buildCostFun(cost)
}
\arguments{
  \item{cost}{The cost of collecting the data. It is a two-element list. The first element is a vector 
  containing all the individual costs (including the group cost). The second element is also a list with
  length \code{p}, where \code{p} equals the number of predictors, and the \code{i}th element has the name
  of the \code{i}th predictor and contains the cost pointers pointing to the cost values in the first 
  element of \code{cost}. User could use more general kind of cost value provided that the user 
  write his(her) own function \code{buildCostFun}}.
}
\value{
  A function \code{costFun(alpha)} is returned.
}
\author{Li Hua Yue}
\seealso{\code{\link{branchLars}}}
\examples{
data(diabetes)
X <- as.matrix(diabetes[, 1:10])

# The X matrix consists of 10 varialbes. 
# Suppose S1 to S6 belong to one group, and we assign the costs to be
# AGE SEX BPI BP   S1 S2 S3 S4 S5 S6 
#                  (group cost: 20) 
#  0   0   5  10   30 30 30 30 30 30   

costs <- list()
costs[[1]] <- c(0, 5, 10, 20, rep(30, 6))      # costs (including group cost)
costs[[2]] <- list(1, 1, 2, 3, c(4,5), c(4,6), # cost pointers for each variable
              c(4,7), c(4,8), c(4,9), c(4,10))
names(costs[[2]]) <- colnames(X)

costFun <- buildCostFun(costs)
alpha <- c(0, rep(1, 5), rep(-1, 4))
names(alpha) <- colnames(X)
costFun(alpha)
}
\keyword{methods}
