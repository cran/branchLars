.packageName <- "branchLars"

"Loss12" <- function(lambda, coef, X, y){
# The first two parts of the loss function
  c(sum((y-X%*%coef)^2), lambda*sum(abs(coef))) 
}

"changeOrderX" <- function(X, cost, Bin, orderLars){
  # if more than one bin, order the variables by cost for different bins 
  # order the variables by lars entry in each bin
  # put zero cost variables at last
  n <- dim(X)[1]       # n observations
  maxB <- floor(max(cost)/Bin)        # maximum bin number 
  if(maxB==0){            # only one bin
     orderX0 <- which(cost == 0)
     orderX <- which((cost >0) & (cost < Bin))
     if(length(orderX)==1){
       X1 <- as.matrix(X[, orderX])
       colnames(X1) <- colnames(X)[orderX]
     }
     else X1 <- X[, orderX]
     if(length(orderX0)==1){
       X2 <- as.matrix(X[, orderX0])
       colnames(X2) <- colnames(X)[orderX0]
     }
     else X2 <- X[, orderX0]
     if(dim(X1)[2]==1) Xnew <- cbind(X1, X2)
     else Xnew <- cbind(X1[, order(orderLars[orderX])], X2)
  } 
  else{ 
    Xnew <- matrix(numeric(), nrow=n)    # the new X whose order will be changed
    for(i in maxB:0){
      orderX <- which((cost >= i*Bin) & (cost < (i+1)*Bin))
      if(length(orderX) > 0){
        if(i==0 && any(cost == 0)){
           orderX0 <- which(cost == 0)
           orderX <- which((cost >0) & (cost < Bin))
           if(length(orderX)==1){
             X1 <- as.matrix(X[, orderX])
             colnames(X1) <- colnames(X)[orderX]
           }
           else X1 <- X[, orderX]
           if(length(orderX0)==1){
             X2 <- as.matrix(X[, orderX0])
             colnames(X2) <- colnames(X)[orderX0]
           }
           else X2 <- X[, orderX0]
           if(dim(X1)[2]==1) Xnew <- cbind(Xnew, X1, X2)
           else Xnew <- cbind(Xnew, X1[, order(orderLars[orderX])], X2)
        }
        else{
           if(length(orderX)==1){
             X1 <- as.matrix(X[, orderX])
             colnames(X1) <- colnames(X)[orderX]
             Xnew <- cbind(Xnew, X1)
           }
           else Xnew <- cbind(Xnew, X[, orderX][, order(orderLars[orderX])])
        }
      }
    }
  }
  Xnew
}

"buildCostFun" <- function(cost){
  function(alpha){
    Name <- names(alpha)
    vIn <- Name[alpha==1]         # selected variables
    vUnknown <- Name[alpha==-1]    # haven't been searched variables
    cIn <- unique(unlist(cost[[2]][vIn]))    # costs of selected variables
    costSel <- sum(cost[[1]][cIn])   # sum of the costs of selected variables

    costNew <- cost[[1]]
    costNew[cIn] <- 0       # let the cost of selected variables be 0
    if(any(is.na(vUnknown))) return(list(costSel, NULL))
    costUnknown <- cost[[2]][vUnknown]     # costs of unknown variables
    if(length(vUnknown)>1){
      costUnknown <- sapply(costUnknown, FUN=function(x){sum(costNew[x])})
    }
    else costUnknown <- sum(costNew[unlist(costUnknown)])
    list(costSel, costUnknown)
  }
}

"standardizeG" <- function(Xs, normX, cost="NULL"){
  nm <- dim(Xs)
  n <- nm[1]    # n obs
  m <- nm[2]    # m variables
  one <- rep(1, n)

  # add the square terms
  Xsquare <- Xs^2
  nameX <- colnames(Xs)
  colnames(Xsquare) <- paste(nameX, "^2", sep="")
  normXsquare <- normX^2
  if(cost[1] != "NULL"){
    cost2 <- cost[[2]]   # add the square term costs 
    names(cost2) <- colnames(Xsquare)
  }

  # Center and scale Xsquare
  meanXsquare <- drop(one%*%Xsquare)/n
  Xsquare <- scale(Xsquare, meanXsquare, FALSE)
  normXsquare2 <- sqrt(drop(one%*%(Xsquare^2)))
  names(normXsquare2) <- NULL
  Xsquare <- scale(Xsquare, FALSE, normXsquare2)

  # add the interaction terms
  Xinter <- as.matrix(Xs[,1]*Xs[, 2])
  normXinter <- normX[1]*normX[2]
  colnames(Xinter) <- paste(nameX[1], ":", nameX[2], sep="")
  if(cost[1] != "NULL"){
    costInter <- list(unique(unlist(cost[[2]][1:2])))  # add the interaction term cost
    names(costInter) <- colnames(Xinter)
  }
  for(i in 2:(m-1)){
    XinterNew <- Xs[,1:i]*Xs[, (i+1)]
    colnames(XinterNew) <- paste(nameX[1:i], ":", nameX[i+1], sep="")
    normXinterNew <- normX[1:i]*normX[i+1]
    Xinter <- cbind(Xinter, XinterNew)
    normXinter <- c(normXinter, normXinterNew)
    if(cost[1] != "NULL"){
      costInterNew <- lapply(apply(cbind(cost[[2]][1:i], cost[[2]][i+1]), 1, unique), unlist)
      costInterNew <- lapply(costInterNew, unique)
      names(costInterNew) <- colnames(XinterNew)
      costInter <- c(costInter, costInterNew)
    }
  }
  
  # Center and scale Xinter
  meanXinter <- drop(one%*%Xinter)/n
  Xinter <- scale(Xinter, meanXinter, FALSE)   
  normXinter2 <- sqrt(drop(one%*%(Xinter^2)))
  names(normXinter2) <- NULL
  Xinter <- scale(Xinter, FALSE, normXinter2)  

  X <- cbind(Xs, Xsquare, Xinter)
  if(cost[1] != "NULL") cost[[2]] <- c(cost[[2]], cost2, costInter)
  normx <- c(normX, normXsquare*normXsquare2, normXinter*normXinter2)
  list(Xs=X, normX=normx, cost=cost)
}


"standardize" <- function(X, y){
  if(is.data.frame(X)) X <- as.matrix(X)
  nm <- dim(X)
  n <- nm[1]     # n obs
  m <- nm[2]     # m variables
  one <- rep(1, n)

  # Center X and y, and scale X, and save the means
  meanx <- drop(one%*%X)/n
  X <- scale(X, meanx, FALSE)   # centers X
  y <- drop(y-mean(y))          # centers y
  normx <- sqrt(drop(one%*%(X^2)))
  names(normx) <- NULL
  X <- scale(X, FALSE, normx)   # scales X

  list(y=y, X=X, normX=normx)
}

"unstandardize" <- function(object, normX){
  coef <- object$coef
  m <- length(normX)
  nameX <- names(coef)
  coef <- scale(t(coef), FALSE, normX)[1:m]  # unstandardize
  names(coef) <- nameX
  t(coef)
}


"branchOpt" <- function(X, y, lambda, M, alpha, step, preSolution, preLarsEntry, 
                      obj0, bestObj, costFun, gamma, preNode, sumC, Bin){
## recursive search function
  if(step==(M+1)) return(bestObj)
  X1 <- X
  X1[, (alpha[1,]==0)] <- 0
  X2 <- X1
  evals <- bestObj$evals 
  tree <- bestObj$tree
  nameOld <- colnames(X)

  tree$Nnode <- tree$Nnode + 2
  tree$nodes <- c(tree$nodes, (tree$Nnode-1), tree$Nnode)
  tree$nodeLabel <- c(tree$nodeLabel, paste(colnames(X)[step],"=",0,sep=""),
                                      paste(colnames(X)[step],"=",1,sep=""))
  tree$edges$from <- c(tree$edges$from, preNode, preNode)
  tree$edges$to <- c(tree$edges$to, (tree$Nnode-1), tree$Nnode)


  ## Solve k(th) Right Relaxation:
  # right bound and loss:
  alphaRight <- alpha
  alphaRightR1 <- alpha[1,]
  if(step < M) alphaRightR1[(step+1):M] <- 0
  costK <- costFun(alphaRightR1)[[2]]
  alphaRight[1, step] <- 1
  costRight <- costFun(alphaRight[1,])
  sumCRight <- costRight[[1]]
  costUnknownVR <- costRight[[2]]
  boundRight <- preSolution$bound + gamma*costK

  alphaRight[2,][(alphaRight[1,])==1] <- 1
  sumCRightTrue <- costFun(alphaRight[2,])[[1]]
  costRight2 <- gamma*sumCRightTrue
  lossRight <- sum(preSolution$loss[2:3]) + costRight2
  alphaRight[2,] <- alpha[2,]

  # right search
  bestObjOld <- bestObj
  treeRight <- tree
  treeRight$path <- c(tree$path, paste(preNode, "~", tree$Nnode, sep=""))
  if(any(boundRight > bestObj$loss[1], sumCRight > sumC, sum(costUnknownVR)==0)){
    # right branch pruned:
    if(any(all(abs(lossRight - bestObj$loss[1])<0.001,
               preSolution$loss[2] < bestObj$loss[2]), 
           (bestObj$loss[1] - lossRight)>=0.001)){
      lossV <- c(lossRight, preSolution$loss[2], preSolution$loss[3], costRight2)
      names(lossV) <- c("Total_Loss", "RSS", "L1_penalty", "Cost")
      bestObj <- list(loss=lossV, coef=preSolution$coef, evals=evals,
                      larsObj=bestObjOld$larsObj)
    }
    bestObjRight <- bestObj
    bestObjRight$tree <- treeRight
  }
  else{
    # no need to change order if step=(M-1):M or the cost for unknown is unchanged
    if((step>=(M-1))|all(costUnknownVR==costFun(alpha[1,])[[2]][names(costUnknownVR)])){
      Xright <- X2 
      coefRight <- preSolution$coef
    }
    else{     # change order 
      orderLars <- preLarsEntry[(step+1):M]
      Xchange <- X2[,(step+1):M]
      Xnew <- changeOrderX(Xchange, costUnknownVR, Bin, orderLars)
      if(step==1){
        Xright <- cbind(X2[,1], Xnew)
        colnames(Xright)[1] <- colnames(X2)[1]
      }
      else Xright <- cbind(X2[,1:step], Xnew)
      nameNew <- colnames(Xright)
      coefRight <- preSolution$coef[nameNew]
      preLarsEntry <- preLarsEntry[nameNew]
      alphaRight <- alphaRight[, nameNew]
    }
    
    lossV <- c(lossRight, preSolution$loss[2], preSolution$loss[3], costRight2)
    names(lossV) <- c("Total_Loss", "RSS", "L1_penalty", "Cost")
    solutionRight <- list(loss=lossV, coef=coefRight, bound=boundRight)
    if(any(all(abs(lossRight - bestObj$loss[1])<0.001, lossV[2] < bestObj$loss[2]),
           (bestObj$loss[1]-lossRight)>=0.001)){
      bestObj <- list(loss=lossV, coef=coefRight, evals=evals,
                      tree=treeRight, larsObj=bestObj$larsObj)
    }
    else  bestObj$tree <- treeRight
    bestObjRight <- branchOpt(Xright,y,lambda,M,alphaRight,step+1,solutionRight,
            preLarsEntry, obj0, bestObj, costFun, gamma, tree$Nnode, sumC, Bin)
  }


  ## Solve k(th) Left Relaxation:
  bestObj <- bestObjRight
  treeLeft <- bestObjRight$tree
  treeLeft$path <- c(tree$path, paste(preNode, "~", (tree$Nnode-1), sep=""))
  evals <- bestObjRight$evals

  alphaLeft <- alpha
  alphaLeft[1,step] <- 0
  # left bound and loss: 
  if(all(alphaLeft[1,]==0)){
    	if(obj0$loss[1] <= bestObj$loss[1]) bestObjLeft <- obj0
    	else bestObjLeft <- bestObj
   		bestObjLeft$tree <- treeLeft
    	bestObjLeft$evals <- evals
  }
  else{
      X1[, step] <- 0
      larsLeft <- lars(X1, y)
      orderLeft <- larsLeft$entry
      names(orderLeft) <- colnames(X1)
      evals <- evals + 1
      coefLeft <- predict(larsLeft, s=lambda, type="coef", mode="lambda")$coef
      alphaLeft[2,] <- 1-(coefLeft==0)
      alphaLeft[2,][(alphaLeft[1,])==1] <- 1
      costLeft <- costFun(alphaLeft[1,])
      sumCLeft <- costLeft[[1]]
      costUnknownVL <- costLeft[[2]]
      costLeft1 <- gamma*sumCLeft
      sumCLeftTrue <- costFun(alphaLeft[2,])[[1]]
      costLeft2 <- gamma*sumCLeftTrue
      lossLeft12 <- Loss12(lambda, coefLeft, X1, y)
      boundLeft <- sum(lossLeft12) + costLeft1
      lossLeft <- sum(lossLeft12) + costLeft2
      alphaLeft[2,] <- 1-(coefLeft==0)

      # next search
      bestObj$evals <- evals
      if(any(boundLeft > bestObj$loss[1], sumCLeft > sumC, sum(costUnknownVL)==0)){
        	# left branch pruned:
        	if(any(all(abs(lossLeft - bestObj$loss[1])<0.001,
                     lossLeft12[1] < bestObj$loss[2]),
               	 (bestObj$loss[1]-lossLeft)>=0.001)){
          		lossV <- c(lossLeft, lossLeft12[1], lossLeft12[2], costLeft2)
          		names(lossV) <- c("Total_Loss", "RSS", "L1_penalty", "Cost")
          		bestObj <- list(loss=lossV, coef=coefLeft, larsObj=larsLeft)
        	}
        	bestObjLeft <- bestObj
        	if(all((alphaLeft[1,][alphaLeft[1,]!=-1])==0) && 
               	 bestObjLeft$loss[1] > obj0$loss[1]){
          		bestObjLeft <- obj0
        	}
        	bestObjLeft$tree <- treeLeft
        	bestObjLeft$evals <- evals
      }
      else{
        # no need to change order if step = (M-1) or M
        if(step>=(M-1))  Xleft <- X1       
        else{             
					  # change order if step < (M-1)
          	orderLars <- orderLeft[(step+1):M]
          	Xchange <- X1[,(step+1):M]
          	Xnew <- changeOrderX(Xchange, costUnknownVL, Bin, orderLars)
          	if(step==1){
            		Xleft <- cbind(X1[,1], Xnew)
            		colnames(Xleft)[1] <- colnames(X1)[1]
          	}
          	else Xleft <- cbind(X1[,1:step], Xnew)
          	nameNew <- colnames(Xleft)
          	coefLeft <- coefLeft[nameNew]
          	orderLeft <- orderLeft[nameNew]  
          	alphaLeft <- alphaLeft[, nameNew]  
        }

        lossV <- c(lossLeft, lossLeft12[1], lossLeft12[2], costLeft2)
        names(lossV) <- c("Total_Loss", "RSS", "L1_penalty", "Cost")
        solutionLeft <- list(loss=lossV, coef=coefLeft, bound=boundLeft)
        if(any(all(abs(lossLeft - bestObj$loss[1])<0.001,
                   lossLeft12[1] < bestObj$loss[2]),
               (bestObj$loss[1]-lossLeft)>=0.001)){
          	bestObj <- list(loss=lossV, coef=coefLeft, evals=evals, 
                          tree=treeLeft, larsObj=larsLeft)
        }
        else bestObj$tree <- treeLeft
        bestObjLeft <- branchOpt(Xleft, y, lambda, M, alphaLeft, step+1, 
                         solutionLeft, orderLeft, obj0, bestObj, costFun, 
                         gamma,(tree$Nnode-1), sumC, Bin)
      }
  }

  # compare left with right
  if((bestObjRight$loss[1] - bestObjLeft$loss[1])>0.001) return(bestObjLeft)
  else if(all(abs(bestObjRight$loss[1] - bestObjLeft$loss[1])<=0.001,
              bestObjRight$loss[2] > bestObjLeft$loss[2])) return(bestObjLeft)
  else{
			if(bestObjRight$coef[step] == 0){
      		bestObjRight$tree <- bestObjLeft$tree
			} 
			else{
      		bestPath <- bestObjRight$tree$path
      		bestObjRight$tree <- bestObjLeft$tree
      		bestObjRight$tree$path <- bestPath
			}
      bestObjRight$evals <- bestObjLeft$evals
      return(bestObjRight)
  }
}




"branchLars" <-
function(X, y, lambda, cost, gamma=1, sumC=0){
## need more than 2 variables (m>2)
## gamma is the user defined parameter used to change weight on cost
  call <- match.call()
  nm <- dim(X)
  n <- nm[1]     # n obs
  m <- nm[2]     # m variables
  gammaOld <- gamma
  gamma <- n*gamma    
  nameOrg <- colnames(X)  
  costFun <- buildCostFun(cost)
  loss0 <- sum(y^2)

  # all variables are in the model (initial lars regression)
  lars1 <- lars(X, y, normalize=FALSE)
  coef1 <- predict(lars1, s=lambda, type="coef", mode="lambda")$coef
  alpha1 <- rep(-1, m)
  names(alpha1) <- nameOrg   
  costUnknownV <- costFun(alpha1)[[2]]
  alpha1 <- rbind(alpha1, (1-(coef1==0)))
  cost1 <- costFun(alpha1[2,])[[1]]
  if(is.na(cost1)) cost1 <- 0
  if(sumC==0) sumC <- cost1
  cost1 <- gamma*cost1
  loss12 <- Loss12(lambda, coef1, X, y)
  bound1 <- sum(loss12)
  loss1 <- bound1 + cost1
  lossV1 <- c(loss1, loss12[1], loss12[2], cost1)
  names(lossV1) <- c("Total_Loss", "RSS", "L1_penalty", "Cost")

  tree <- list(Nnode=1, nodes="1", nodeLabel="root", edges=NULL, path=NULL)
  obj1 <- list(loss=lossV1, coef=coef1, evals=1, tree=tree, larsObj=lars1)
  solution1 <- list(loss=lossV1, coef=coef1, bound=bound1)
 

  if(gamma == 0) result <- obj1
  else{
    # change the order of the variables   
    orderLars1 <- lars1$entry
    names(orderLars1) <- nameOrg
    Bin <- 10/(log(1+gammaOld)*log(1+lambda))*(loss0/(n-1))

    # change the order of the variables 
    Xnew <- changeOrderX(X, costUnknownV, Bin, orderLars1)
    nameNew <- colnames(Xnew)
    orderLars1 <- orderLars1[nameNew]  
    alpha1 <- alpha1[,nameNew]
    coef1 <- coef1[nameNew]

    # no variable in the model:
    lossV0 <- c(loss0, loss0, 0, 0)
    names(lossV0) <- c("Total_Loss", "RSS", "L1_penalty", "Cost")
    coef0 <- rep(0, m)
    names(coef0) <- nameOrg
    obj0 <- list(loss=lossV0, coef=coef0, evals=0, tree=tree, larsObj=NULL)

    result <- branchOpt(Xnew, y, lambda, m, alpha1,1, solution1, orderLars1, obj0, 
                      obj1, costFun, gamma, 1, sumC, Bin)

    result$coef <- result$coef[nameOrg]
  }
  result$sumC <- costFun(1-(result$coef==0))[[1]]
  result$call <- call
  result$sigma2 <- as.numeric(attr(lars1$Cp, "sigma2"))
  result$n <- n
  class(result) <- "branchLars"  
  result
}


"print.branchLars" <- function(x, ...){
  cat("\nCall:\n")
  dput(x$call)
  cat("\nRegression Coefficient:\n")
  print(t(x$coef)) 
}
 
"summary.branchLars" <- function(object, ...){
  cat("\nCall:\n")
  dput(object$call)
  cat("\nOptimal Total Loss and its Components:\n")
  print(object$loss)
}

"drawTree" <- function(object){
  tree <- object$tree
  if(tree$Nnode < 3) return("Less than three nodes")
  Node <- tree$nodes
  nodeLabel <- tree$nodeLabel
  names(nodeLabel) <- Node
  path <- tree$path
  edges <- tree$edges

  g <- new("graphNEL", nodes=Node, edgemode="directed") 
  g <- addEdge(from=as.character(edges$from), to=as.character(edges$to), graph=g) 
  nAttrs <- list()
  nAttrs$label <- nodeLabel
  eAttrs <- list()
  eColor <- rep("red", length(path))
  names(eColor) <- path
  eAttrs$color <- eColor
  attrs <- list(node=list(shape="box", fixedsize=FALSE))
  attrs$node$fontcolor <- "dark blue"
  plot(g, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
}

"predict.branchLars" <- function(object, newdata, ...){
  coef <- object$coef
  if(is.data.frame(newdata)) newdata <- as.matrix(newdata)
  newdata%*%coef
}

"bic" <- function(object){
  RSS <- as.numeric(object$loss[2])
  df <- sum(object$coef != 0)+1
  n <- object$n
  sigma2 <- object$sigma2
  RSS/sigma2 + log(n)*df
}


"cp" <- function(object){
  RSS <- as.numeric(object$loss[2])
  df <- sum(object$coef != 0)+1
  n <- object$n
  sigma2 <- object$sigma2
  RSS/sigma2 - n + 2*df
}

"buildChoose" <- function(X, y, cost, gamma, method){
  function(lambda, sumC=0){
    if(method=="CV"){
      n <- dim(X)[1]
      set.seed(123)
      folds <- split(sample(1:n), rep(1:10, length=n))
      cvs <- c() 
      sumCs <- c()
      for(i in 1:10){
        omit <- folds[[i]]
        bLars <- branchLars(X[-omit,], y[-omit], lambda=lambda,
                            cost=cost, gamma=gamma, sumC=sumC)
        fit <- predict.branchLars(bLars, X[omit,])
        cvs <- c(cvs, mean((y[omit]-fit)^2))
        sumCs <- c(sumCs, bLars$sumC)
      }
      return(c(mean(cvs), max(sumCs)))
    }
    else{
      bLars <- branchLars(X, y, lambda=lambda, cost=cost,
                         gamma=gamma, sumC=sumC)
      if(method=="BIC") return(c(bic(bLars), bLars$sumC))
      if(method=="Cp") return(c(cp(bLars), bLars$sumC))
    }
  }
}

"lambdaOpt" <- function(X, y, cost, gamma=1, lower=0, upper=Inf, method="BIC", 
                       tol=.Machine$double.eps^0.25/3){
  if((lower < 0) | (upper < 0)) return("The lambda value can not be negative")
  if(upper==Inf){                          # upper bound
    lambdas <- lars(X, y)$lambda
    upper <- lambdas[1]
  }

  # objective function
  opt <- buildChoose(X, y, cost, gamma, method) 

  # initialization
  optL <- opt(lower, sumC=0)
  optU <- opt(upper, sumC=optL[2])
  r <- (sqrt(5)-1)/2
  x1 <- lower + (1-r)*(upper-lower)
  x2 <- lower + r*(upper-lower)
  opt1 <- opt(x1)
  opt2 <- opt(x2)
 
  # main loop - golden section search
  while((upper-lower) >= tol){
    if(all(optL[1]==optU[1], optL[1]==opt1[1], optL[1]==opt2[1])){
       if(method=="BIC") 
          return(list(Optimal_Lambda=upper, BIC=optU[1], sumC=optU[2]))
       if(method=="Cp")
          return(list(Optimal_Lambda=upper, Cp=optU[1], sumC=optU[2]))
       return(list(Optimal_Lambda=upper, CV=optU[1], sumC=optU[2]))
    }
    if(opt1[1] <= opt2[1]){
      upper <- x2
      optU <- opt2
      x2 <- x1
      opt2 <- opt1
      x1 <- lower + (1-r)*(upper-lower)
      opt1 <- opt(x1, sumC=optL[2])
    }
    else{
      lower <- x1
      optL <- opt1
      x1 <- x2
      opt1 <- opt2
      x2 <- lower + r*(upper-lower)
      opt2 <- opt(x2, sumC=opt1[2])
    }
  }
  opts <- c(optL[1], opt1[1], opt2[1], optU[1])
  sumCs <- c(optL[2], opt1[2], opt2[2], optU[2])
  Opt <- min(opts)
  x <- c(lower, x1, x2, upper)[which.min(opts)]
  sumC <- sumCs[which.min(opts)]
  if(method=="BIC") return(list(Optimal_Lambda=x, BIC=Opt, sumC=sumC))
  if(method=="Cp") return(list(Optimal_Lambda=x, Cp=Opt, sumC=sumC))
  return(list(Optimal_Lambda=x, CV=Opt, sumC=sumC))
}

