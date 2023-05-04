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


drest_atev0<-function(trt,x,var1,yout,var0=var1,varp=var1,ci="asymptotic",level=.95){
  
  #fitted models 
  
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
  
  #asymptotic variance and CI
  if (ci=="asymptotic"){
  varterm1<-(trt*yout-yhatout1*(trt-ps))/ps-((1-trt)*yout+yhatout0*(trt-ps))/(1-ps)
  meansq<-sum((varterm1-drdelta)^2)
  se<-sqrt(meansq/length(yout)^2)
  ci<-drdelta+qnorm((1-level)/2+level)*c(-se,se)}
  
  return(c(drdelta,se,ci))
  
}

