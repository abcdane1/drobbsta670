#' Doubly Robust Estimation of ATE 
#'
#' This package provides a doubly robust point estimate and confidence interval for the average treatment effect (ATE) with binary treatment. The outcome models are assumed to be linear, where
#' the outcome for those under treatment 1 and the outcome for those under treatment 0 are regressed against user specified subsets of non-treatment covariates. The propensity score is modeled using logistic regression, 
#' where treatment status is regressed against a subset of non-treatment user-specified covariates. For confidence intervals, the user has the option to specify whether they wish to have the asymptotic interval, 
#' the basic (or empirical) bootstrap, or the percentile bootstrap.
#' 
#'
#' @param trt  A vector of binary treatments.
#' @param x A matrix, data frame, or vector of non-treatment covariates.
#' @param var1 The numbers or names of columns to be included in outcome model for treatment 1.
#' @param var0 The numbers or names of columns to be included in outcome model for treatment 0. If `var0` is not specified, `var1=var0`.
#' @param yout A vector of outcome variables.
#' @param varp The numbers or names of columns to be included in propensity score model. If `varp` is not specified, `var1=varp`. 
#' 
#'
#' @return The return value is an object of...
#'
#' @examples
#'
#' @export

drest_ate<-function(trt,x,var1,yout,var0=var1,varp=var1,ci="asymptotic",level=.95,B=1000){
  
  drest_sub<-function(trt,x,var1,yout,var0,varp){
    #ensuring vectors are matrices for indexing
    x<-as.matrix(x)
    
    #treatment indexing
    trtind1<-which(trt==1)
    trtind0<-which(trt==0)
    
    #outcome models
    x1<-x[trtind1,var1]
    yout1<-yout[trtind1]
    
    x0<-x[trtind0,var0]
    yout0<-yout[trtind0]
    
    #outcome models
    
    out1fit<-lm(yout1~x1)
    out0fit<-lm(yout0~x0)
    
    #intercept 
    int<-rep(1,length(yout))
    fullx1<-cbind(int,x[,var1])
    fullx0<-cbind(int,x[,var0])
    yhatout1<-as.vector(fullx1%*%out1fit$coefficients)
    yhatout0<-as.vector(fullx0%*%out0fit$coefficients)
    
    #propensity score
    ps<-glm(trt~x[,varp],family="binomial")$fitted
    
    #dr1
    dr1<-mean(trt*yout/ps-(trt-ps)/ps*yhatout1)
    #dr0
    dr0<-mean((1-trt)*yout/(1-ps)+(trt-ps)/(1-ps)*yhatout0)
    
    #delta of estimator
    drdelta<-dr1-dr0
    
    return(list(yhatout1,yhatout0,ps,drdelta))
  }
  
  
  #asymptotic variance and CI
  if (ci=="asymptotic"){
    yhatout1<-drest_sub(trt,x,var1,yout,var0,varp)[[1]]
    yhatout0<-drest_sub(trt,x,var1,yout,var0,varp)[[2]]
    ps<-drest_sub(trt,x,var1,yout,var0,varp)[[3]]
    drdelta<-drest_sub(trt,x,var1,yout,var0,varp)[[4]]
    
    varterm1<-(trt*yout-yhatout1*(trt-ps))/ps-((1-trt)*yout+yhatout0*(trt-ps))/(1-ps)
    meansq<-sum((varterm1-drdelta)^2)
    se<-sqrt(meansq/length(yout)^2)
    ci<-drdelta+qnorm((1-level)/2+level)*c(-se,se)
    return(c(drdelta,se,ci))}
  
  #basic bootstrap  
  if (ci=="basic"){
    drdelta<-drest_sub(trt,x,var1,yout,var0,varp)[[4]]
    boot<-rep(0,B)
    for (j in 1:B){
      sampind<-sample(1:length(trt),size=length(trt),replace=TRUE)
      trtnew<-trt[sampind]
      xnew<-as.matrix(x)[sampind,]
      youtnew<-yout[sampind]
      boot[j]<-drest_sub(trtnew,xnew,var1,youtnew,var0,varp)[[4]]
    }
    se<-sd(boot,na.rm=TRUE)
    ci<-drdelta+qnorm((1-level)/2+level)*c(-se,se)
    nacount<-length(which(is.na(boot)))
    if(nacount/B>.1){
      warning(paste("Unreliable interval: model cannot be fit to",nacount,"bootstrap samples, which is more than 10% of the tried samples"))
    }
    return(c(drdelta,se,ci))}
  
  #percentile boot 
  if (ci=="percentile"){
    drdelta<-drest_sub(trt,x,var1,yout,var0,varp)[[4]]
    boot<-rep(0,B)
    for (j in 1:B){
      sampind<-sample(1:length(trt),size=length(trt),replace=TRUE)
      trtnew<-trt[sampind]
      xnew<-as.matrix(x)[sampind,]
      youtnew<-yout[sampind]
      boot[j]<-drest_sub(trtnew,xnew,var1,youtnew,var0,varp)[[4]]
    }
    se<-"NA: Percentile Interval"
    nacount<-length(which(is.na(boot)))
    ci<-quantile(boot,c((1-level)/2,(1-level)/2+level),na.rm=TRUE)
    if(nacount/B>.1){
      warning(paste("Unreliable interval: model cannot be fit to",nacount,"bootstrap samples, which is more than 10% of the tried samples"))
    }
    return(c(drdelta,ci))}
  
  #bca
  if (ci=="bca"){
    drdelta<-drest_sub(trt,x,var1,yout,var0,varp)[[4]]
    boot<-rep(0,B)
    jackrep<-length(trt)
    jacktheta<-rep(0,jackrep)
    for (i in 1:jackrep){
      ind<-setdiff(1:jackrep,i)
      trtj<-trt[ind]
      xj<-as.matrix(x)[ind,]
      youtj<-yout[ind]
      jacktheta[i]<-drest_sub(trtj,xj,var1,youtj,var0,varp)[[4]]}
    for (j in 1:B){
      sampind<-sample(1:length(trt),size=length(trt),replace=TRUE)
      trtnew<-trt[sampind]
      xnew<-as.matrix(x)[sampind,]
      youtnew<-yout[sampind]
      boot[j]<-drest_sub(trtnew,xnew,var1,youtnew,var0,varp)[[4]]
    }
    meanjack<-mean(jacktheta,na.rm=TRUE)
    nacountjack<-length(which(is.na(jacktheta)))
    jacktheta<-jacktheta[!is.na(jacktheta)]
    nacount<-length(which(is.na(boot)))
    alpha<-sum((jacktheta-meanjack)^3)/(6*sum((jacktheta-meanjack)^2)^1.5)
    z0<-qnorm(mean(boot<drdelta,na.rm=TRUE))
    la<-pnorm(z0+(z0+qnorm((1-level)/2))/(1-alpha*(z0+qnorm((1-level)/2))))
    lb<-pnorm(z0+(z0-qnorm((1-level)/2))/(1-alpha*(z0-qnorm((1-level)/2))))
    ci<-quantile(boot,c(la,lb),na.rm=TRUE)
    if(nacountjack/jackrep >.1){
      warning(paste("Unreliable bias correction estimation: jacknife model cannot be fit to",nacountjack,"samples, which is more than 10% of the tried samples"))}
    if(nacount/B>.1){
      warning(paste("Unreliable interval: model cannot be fit to",nacount,"bootstrap samples, which is more than 10% of the tried samples"))
    }
    return(c(drdelta,ci))}
}
