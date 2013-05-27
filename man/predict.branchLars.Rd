\name{predict.branchLars}
\alias{predict.branchLars}
\title{Makeing Prediction Using New Dataset}
\description{
  \code{predict.branchLars} predicts the response value using a new dataset  
 based on an optimal model returned by branchLars.
}
\usage{
\method{predict}{branchLars}(object, newdata, \dots)
}
\arguments{
  \item{object}{the branchLars object returned by the function \code{branchLars}}
  \item{newdata}{a new dataset used to predict future values}
  \item{\dots}{additional arguments for predict generic}
}
\value{
  Given a new dataset, it returns the predict values based on an optimal model returned by \code{branchLars}.
}
\author{Li Hua Yue}
\seealso{\code{\link{branchLars}}}
\keyword{utilities}

