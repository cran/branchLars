\name{branchLars-internal}
\alias{Loss12}
\alias{changeOrderX}
\alias{branchOpt}
\alias{buildChoose}
\title{Internal BranchLars Functions}
\description{
  Internal branchLars functions
}
\usage{
Loss12(lambda, coef, X, y)
changeOrderX(X, cost, Bin, orderLars)
branchOpt(X, y, lambda, M, alpha, step, preSolution, preLarsEntry, obj0,
          bestObj, costFun, gamma, preNode, sumC, Bin)
buildChoose(X, y, cost, gamma, method)
}
\details{
  These are not to be called by the user. The \code{Loss12} function calculates the lasso loss 
  (residual sum of squares and the loss due to L1 panelty). Function \code{changeOrderX} change 
  the order of the variables in X matrix based on their costs and LARS entries. Function \code{branchOpt} 
  is the recursive funtions inplementing Branching LARS algorithm for \code{branchLars}. 
  The \code{buildChoose} function builds a function using only \code{lambda} and \code{sumC} 
  as inputs and outputting the BIC(or Cp) and \code{sumC} value of a new branchLars object. This function 
  is called by the function \code{lambdaOpt} which searches for the optimal value of \code{lambda} 
  based on BIC(or Cp) criteria.
}
\author{Li Hua Yue}
\keyword{internal}
