#' Doubly Robust Estimation of ATE
#'
#' This package provides a doubly robust point estimate and confidence interval for the average treatment effect (ATE) with binary treatment. The outcome models are assumed to be linear, where
#' the outcome for those under treatment 1 and the outcome for those under treatment 0 are regressed against user specified subsets of non-treatment covariates. The propensity score is modeled using logistic regression, 
#' where treatment status is regressed against a subset of non-treatment user-specified covariates. For confidence intervals, the user has the option to specify whether they wish to have the asymptotic interval, 
#' or non-parametric bootstrap intervals, which include the basic/empirical, percentile, and bias corrected and accelerated (BCa) intervals.
#' 
#'
#' @param trt  A vector of binary treatments.
#' @param x A matrix, data frame, or vector of non-treatment covariates.
#' @param var1 A vector of numbers or names of columns to be included in outcome model for treatment 1.
#' @param var0 A vector of numbers or names of columns to be included in outcome model for treatment 0. If `var0` is not specified, `var1=var0`.
#' @param yout A vector of outcome variables.
#' @param varp A vector of numbers or names of columns to be included in propensity score model. If `varp` is not specified, `var1=varp`. 
#' @param ci The choice of confidence interval. Possible values are the `asymptotic` (the default) interval, and non-parametric bootstraps `basic`, 
#' `percentile`, and `bca` interval. If `ci=asymptotic`, the se reported is the asymptotic approximation. If `ci=basic`, the se reported is of the bootstrapped sample.
#' @param level A value between 0 and 1 which gives confidence level of the confidence interval - default is .95.
#' @param B A postive integer numbers of bootstrap samples - default is 1000. 
#'
#' @return The return value is an object of...
#'
#'
#' @references Reference
#'
#'
#' @export

drest_ate<-function(trt,x,var1,yout,var0=var1,varp=var1,ci="asymptotic",level=.95,B=1000){
  
  #sub function which fits models and finds dr estimates
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
    
    #linear fits
    out1fit<-lm(yout1~x1)
    out0fit<-lm(yout0~x0)
    
    #intercept 
    int<-rep(1,length(yout))
    fullx1<-cbind(int,x[,var1])
    fullx0<-cbind(int,x[,var0])
    
    #fitted values on all data
    yhatout1<-as.vector(fullx1%*%out1fit$coefficients)
    yhatout0<-as.vector(fullx0%*%out0fit$coefficients)
    
    #propensity score model with fitted
    ps<-glm(trt~x[,varp],family="binomial")$fitted
    
    #dr1
    dr1<-mean(trt*yout/ps-(trt-ps)/ps*yhatout1)
    #dr0
    dr0<-mean((1-trt)*yout/(1-ps)+(trt-ps)/(1-ps)*yhatout0)
    
    #delta of estimator
    drdelta<-dr1-dr0
    
    #outputs estim as well as values for variance 
    return(list(yhatout1,yhatout0,ps,drdelta))
  }
  
  names<-c("Estimate","StdError","Lower","Upper")
  
  #asymptotic variance and CI
  if (ci=="asymptotic"){
    #needed for asymptotic variance calculation
    yhatout1<-drest_sub(trt,x,var1,yout,var0,varp)[[1]]
    yhatout0<-drest_sub(trt,x,var1,yout,var0,varp)[[2]]
    ps<-drest_sub(trt,x,var1,yout,var0,varp)[[3]]
    drdelta<-drest_sub(trt,x,var1,yout,var0,varp)[[4]]
    
    #asymptotic se Lunceford 
    varterm1<-(trt*yout-yhatout1*(trt-ps))/ps-((1-trt)*yout+yhatout0*(trt-ps))/(1-ps)
    meansq<-sum((varterm1-drdelta)^2)
    se<-sqrt(meansq/length(yout)^2)
    
    #asymptotic normal interval 
    ci<-drdelta+qnorm((1-level)/2+level)*c(-se,se)
    output<-c(drdelta,se,ci)
    names(output)<-names
    return(output)
  }
  
  #basic bootstrap  
  if (ci=="basic"){
    drdelta<-drest_sub(trt,x,var1,yout,var0,varp)[[4]]
    #bootstrap samples
    boot<-rep(0,B)
    for (j in 1:B){
      sampind<-sample(1:length(trt),size=length(trt),replace=TRUE)
      trtnew<-trt[sampind]
      xnew<-as.matrix(x)[sampind,]
      youtnew<-yout[sampind]
      boot[j]<-drest_sub(trtnew,xnew,var1,youtnew,var0,varp)[[4]]
    }
    
    #standard deviation of bootstrap
    se<-sd(boot,na.rm=TRUE)
    
    #basic normal ci 
    ci<-drdelta+qnorm((1-level)/2+level)*c(-se,se)
    
    #spits out warning if two many bootstraps estimates cannot be obtained
    nacount<-length(which(is.na(boot)))
    if(nacount/B>.1){
      warning(paste("Unreliable interval: model cannot be fit to",nacount,"bootstrap samples, which is more than 10% of the tried samples"))
    }
    
    #outputs dr estim, se and ci 
    output<-c(drdelta,se,ci)
    names(output)<-names
    return(output)}
  
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
    #no se for percentile interval
    se<-"NA: Percentile Interval"
    
    #keeping track of nas
    nacount<-length(which(is.na(boot)))
    
    #2.5 and 97.5% quantile interval
    ci<-quantile(boot,c((1-level)/2,(1-level)/2+level),na.rm=TRUE)
    
    #spits out warning if two many bootstraps estimates cannot be obtained
    if(nacount/B>.1){
      warning(paste("Unreliable interval: model cannot be fit to",nacount,"bootstrap samples, which is more than 10% of the tried samples"))
    }
    
    #outputs dr estim, se (as NA) and ci 
    output<-c(drdelta,NA,ci)
    names(output)<-names
    return(output)}
  
  #bca interval
  if (ci=="bca"){
    drdelta<-drest_sub(trt,x,var1,yout,var0,varp)[[4]]
    boot<-rep(0,B)
    jackrep<-length(trt)
    jacktheta<-rep(0,jackrep)
    #jackknife estimator for acceleration factor 
    for (i in 1:jackrep){
      #leaving one out 
      ind<-setdiff(1:jackrep,i)
      trtj<-trt[ind]
      xj<-as.matrix(x)[ind,]
      youtj<-yout[ind]
      #fitting model to leave one out
      jacktheta[i]<-drest_sub(trtj,xj,var1,youtj,var0,varp)[[4]]}
    #bootstrap samples
    for (j in 1:B){
      sampind<-sample(1:length(trt),size=length(trt),replace=TRUE)
      trtnew<-trt[sampind]
      xnew<-as.matrix(x)[sampind,]
      youtnew<-yout[sampind]
      boot[j]<-drest_sub(trtnew,xnew,var1,youtnew,var0,varp)[[4]]
    }
    #acceleration factor
    #mean of jackknife
    meanjack<-mean(jacktheta,na.rm=TRUE)
    nacountjack<-length(which(is.na(jacktheta)))
    jacktheta<-jacktheta[!is.na(jacktheta)]
    #keeping track of nas
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
    #outputs standard dr estim, se (as NA) and ci
    output<-c(drdelta,NA,ci)
    names(output)<-names
    return(output)}
}
