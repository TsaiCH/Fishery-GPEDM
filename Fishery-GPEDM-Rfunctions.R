##########################################################################################################
getYardstickMetrics_fn <- function(results,msg=TRUE) {
  require(yardstick)
  if(!is.data.frame(results)) {
    message("Input is not data.frame format in getYardstickMetrics_fn") 
  }
  numeric_mets <- yardstick::metric_set(mae,mape,rmse,rsq)
  out <- numeric_mets(results,obs,pred)
  return(out)
}
##########################################################################################################
getR2_fn <- function(obs,pred,msg=FALSE) {
  if(msg&&(any(is.na(obs))||any(is.na(pred)))) {
    message("NA in observations or predictions are omitted in getR2_fn") 
  }
  dat <- na.omit(cbind(obs,pred))
  R2 <- 1-sum((dat[,1]-dat[,2])^2)/sum((dat[,1]-mean(dat[,1]))^2)
  return(R2)
}
##########################################################################################################
getKernelDist_fn <- function(X,Y,w=NULL) {
  if(is.null(w)) {
    stop("Error: length-scale parameters are not defined in getKernelDist_fn")
  }
  if(is.null(dim(X))) {
    X <- matrix(X,nrow=1)
  }
  if(is.null(dim(Y))) {
    Y <- matrix(Y,nrow=1)
  }
  D <- matrix(NA,nrow(X),nrow(Y))
  for(i in 0:(nrow(X)-1)) {
    for(j in 0:(nrow(Y)-1)) {
      v1 <- abs(X[i+1,]-Y[j+1,])
      v2 <- (w*v1)*(w*v1)
      D[i+1,j+1] <- sum(v2)
    }
  }
  return(D)
}
##########################################################################################################
getPosteriorGP_example_fn <- function(logest,K=NULL,ysigma2,X,Y,XX,msg=FALSE) {
  if(is.null(dim(XX))) { 
    XX <- matrix(XX,nrow=1)
    if(msg) {
      message("Input XX is 1-D numeric and is converted to a 1-row matrix in getPosteriorGP_fn")
    }
  }
  #
  phi <- exp(logest[1:ncol(X)])
  tau <- (ysigma2-0.001)/(1+exp(-logest[ncol(X)+1]))+0.001
  g <- (ysigma2-0.001)/(1+exp(-logest[ncol(X)+2]))+0.001
  #
  if(is.null(K)) {   
    D <- getKernelDist_fn(X=X,Y=X,w=phi)
    K <- tau*exp(-D)+diag(g,nrow(D))
  }
  Ki <- solve(K)
  DX <- getKernelDist_fn(X=XX,Y=X,w=phi)
  KX <- tau*exp(-DX)
  DXX <- getKernelDist_fn(X=XX,Y=XX,w=phi)
  KXX <- tau*exp(-DXX)
  #
  m <- KX%*%Ki%*%Y
  sigma2 <- KXX-KX%*%Ki%*%t(KX)
  return(list("m"=m,"sigma2"=sigma2))
}
###########################################################################################################
getLOOfromGP_fn <- function(logest,ysigma2,X,y) {
  phi <- exp(logest[1:ncol(X)])
  tau <- (ysigma2-0.001)/(1+exp(-logest[ncol(X)+1]))+0.001
  g <- (ysigma2-0.001)/(1+exp(-logest[ncol(X)+2]))+0.001
  n <- length(y)
  nc <- ncol(X)
  loopred_m <- loopred_sigma2 <- rep(NA,n)
  for(i in 1:n){
    predX <- matrix(X[i,],ncol=nc)
    trainX <- matrix(X[-i,],ncol=nc)
    predY <- y[i]
    trainY <- y[-i]
    D <- getKernelDist_fn(X=trainX,Y=trainX,w=phi)
    K <- tau*exp(-D)+diag(g,nrow(D))
    Ki <- solve(K)
    DX <- getKernelDist_fn(X=predX,Y=trainX,w=phi)
    KX <- tau*exp(-DX)
    DXX <- getKernelDist_fn(X=predX,Y=predX,w=phi)
    KXX <- tau*exp(-DXX)
    m <- KX%*%Ki%*%trainY
    sigma2 <- KXX-KX%*%Ki%*%t(KX)
    loopred_m[i] <- m
    loopred_sigma2[i] <- diag(sigma2)
  }
  return(list("m"=loopred_m,"sigma2"=loopred_sigma2))
}
###########################################################################################################
getSequentialGP_fn <- function(logest,ysigma2,X,y,split=0.8,finalPred=FALSE) {
  phi <- exp(logest[1:ncol(X)])
  tau <- (ysigma2-0.001)/(1+exp(-logest[ncol(X)+1]))+0.001
  g <- (ysigma2-0.001)/(1+exp(-logest[ncol(X)+2]))+0.001
  n <- length(y)
  nc <- ncol(X)
  flag_index <- floor(n*split)
  end_index <- n-floor(n*split)
  seqpred_m <- seqpred_sigma2 <- rep(NA,n+1)
  if(finalPred) {
    X <- matrix(rbind(X,matrix(c(y[length(y)],tail(X,1)[1:(ncol(X)-1)]),nrow=1)),ncol=ncol(X))
    end_index <- end_index+1
  }
  for(i in 1:end_index) {
    trainX <- matrix(X[1:flag_index,],ncol=nc)
    predX <- matrix(X[(flag_index+1),],ncol=nc)
    trainY <- y[1:flag_index]
    predY <- y[flag_index+1]
    D <- getKernelDist_fn(X=trainX,Y=trainX,w=phi)
    K <- tau*exp(-D)+diag(g,nrow(D))
    Ki <- solve(K)
    DX <- getKernelDist_fn(X=predX,Y=trainX,w=phi)
    KX <- tau*exp(-DX)
    DXX <- getKernelDist_fn(X=predX,Y=predX,w=phi)
    KXX <- tau*exp(-DXX)
    m <- KX%*%Ki%*%trainY
    sigma2 <- KXX-KX%*%Ki%*%t(KX)
    seqpred_m[flag_index+1] <- m
    seqpred_sigma2[flag_index+1] <- diag(sigma2)
    flag_index <- flag_index+1
  }
  #out-of-sample indices included in the prediction set
  seqpredIndex <- which(!is.na(seqpred_m[-(n+1)])) 
  if(!finalPred) {
    seqpred_m <- seqpred_m[-(n+1)]
    seqpred_sigma2 <- seqpred_sigma2[-(n+1)]
    return(list("m"=seqpred_m,"sigma2"=seqpred_sigma2,"seqpredIndex"=seqpredIndex))
  } else {
    final_m <- tail(seqpred_m,1)
    final_sigma2 <- tail(seqpred_sigma2,1)
    seqpred_m <- seqpred_m[-(n+1)]
    seqpred_sigma2 <- seqpred_sigma2[-(n+1)]
    return(list("m"=seqpred_m,"sigma2"=seqpred_sigma2,"seqpredIndex"=seqpredIndex,
                "final_m"=final_m,"final_sigma2"=final_sigma2))
  }
}
###########################################################################################################
transformY_fn <- function(cpue,E,op) {
  switch(op,
         "original"=embed(cpue,E)[,1],
         "log"=embed(log(cpue),E)[,1],
         "logdifference"=embed(log(cpue),E)[,1]-embed(log(cpue),E)[,2],
         stop("Unknown operation for transformY")
  )
}
###########################################################################################################
opmethod_fn <- function(obj,op) {
  switch(op,
         "L-BFGS-B"= optim(obj$par,obj$fn,obj$gr,method="L-BFGS-B"),
         "nlminb"= nlminb(obj$par,obj$fn,obj$gr),
         "rprop"= CaDENCE::rprop(obj$par,obj$fn),
         stop("Unknown operation for optimization method")
  )
}
######################################################################################################################
fisheryGPEDMfit_example_fn <- function(cpue=NULL,catch=NULL,E=5,method="L-BFGS-B",split=0.8,transformY="original") {
  xcpue <- as.matrix(embed(cpue,E)[,-1])
  xcatch <- as.matrix(embed(catch,E)[,-1])
  y <- transformY_fn(cpue=cpue,E=E,op=transformY)
  ysigma2 <- var(y)
  DLLname <- "gpEDMtmb_v5_linearscaling"
  dyn.load(TMB::dynlib("gpEDMtmb_v5_linearscaling"))
  data_tmb <- list("xcpue"=as.matrix(xcpue),"xcatch"=as.matrix(xcatch),"y"=as.vector(y),"ysigma2"=as.numeric(ysigma2))
  parms_tmb <- list("logW"=log(rep(1.01,ncol(xcpue))),"logtau"=log(ysigma2*0.1),"logg"=log(ysigma2*0.1),"logb"=log(0.01))
  obj <- TMB::MakeADFun(data=data_tmb,parameters=parms_tmb,DLL=DLLname,hessian=TRUE)
  out <- opmethod_fn(obj=obj,op=method)
  rep <- TMB::sdreport(obj)
  nLL <- out$value 
  logest <- out$par
  Kest <- obj$report()$K
  qEst <- as.numeric(exp(logest["logb"]))
  XX <- X <- as.matrix(xcpue)-qEst*as.matrix(xcatch)
  insamplePred <- getPosteriorGP_example_fn(logest=logest,K=Kest,ysigma2=as.numeric(ysigma2),X=X,Y=y,XX=XX,msg=FALSE)
  obs <- y
  insampleR2 <- getR2_fn(obs=obs,pred=insamplePred$m)
  looPred <- getLOOfromGP_fn(logest=logest,ysigma2=as.numeric(ysigma2),X=X,y=y) 
  looR2 <- getR2_fn(obs=obs,pred=looPred$m)
  looSE <- sqrt(looPred$sigma2)
  seqPred <- getSequentialGP_fn(logest=logest,ysigma2=as.numeric(ysigma2),X=X,y=y,split=split,finalPred=FALSE)
  seqR2 <- getR2_fn(obs=obs[seqPred$seqpredIndex],pred=seqPred$m[seqPred$seqpredIndex])
  seqSE <- sqrt(seqPred$sigma2)
  #
  parmslist <- list()
  parmslist$logest <- logest
  parmslist$K <- Kest
  parmslist$ysigma2 <- as.numeric(ysigma2)
  parmslist$X <- X
  parmslist$y <- y
  parmslist$qEst <- qEst
  parmslist$nLL <- nLL
  parmslist$inssamplepred <- insamplePred
  return(list(rep=rep,obs=obs,looPred=looPred$m,looSE=looSE,seqPred=seqPred$m,
              seqSE=seqSE,insampleR2=insampleR2,looR2=looR2,seqR2=seqR2,parms=parmslist))
}
######################################################################################################################
