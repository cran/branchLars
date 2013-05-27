\name{diabetes}
\alias{diabetes}
\docType{data}
\title{Blood and Other Measurements in Diabetics}
\description{
  The diabetes data frame has 442 rows and 11 columns. These are the data used 
  in the Efron et al "Least Angle Regression" paper.
}
\usage{data(diabetes)}
\format{
  A data frame with 442 observations on the following 11 variables.
  \describe{
    \item{\code{AGE}}{age of the pations}
    \item{\code{SEX}}{a numeric vector}
    \item{\code{BMI}}{body mass index}
    \item{\code{BP}}{average blood pressure}
    \item{\code{S1}}{blood serum measurement: tc}
    \item{\code{S2}}{blood serum measurement: ldl}
    \item{\code{S3}}{blood serum measurement: hdl}
    \item{\code{S4}}{blood serum measurement: tch}
    \item{\code{S5}}{blood serum measurement: ltg}
    \item{\code{S6}}{blood serum measurement: glu}
    \item{\code{Y}}{response}
  }
}

\source{
  \url{http://www-stat.stanford.edu/~hastie/Papers/LARS/diabetes.data}
}
\examples{
data(diabetes)
X <- as.matrix(diabetes[, 1:10])
y <- diabetes$Y
Xy <- standardize(X,y)
Xs <- Xy$X
ys <- Xy$y

costs <- list()
costs[[1]] <- c(0, 5, 10, 20, rep(30, 6))      
costs[[2]] <- list(1, 1, 2, 3, c(4,5), c(4,6), c(4,7),c(4,8), c(4,9), c(4,10))
names(costs[[2]]) <- colnames(X)

result1 <- branchLars(Xs, ys, lambda=10, costs, gamma=1)
unstandardize(result1, Xy$normX)
}
\keyword{datasets}

