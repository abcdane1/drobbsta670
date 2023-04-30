#' Double Robust Estimation of ATE 
#'
#' This package provides a doubly robust point estimate and confidence interval for the average treatment effect (ATE) with binary treatment. The outcome models are assumed to be linear, where
#' the outcomes for those under treatment 1 and the outcomes for those under treatment 0 are regressed against user specified subsets of non-treatment covariates. The propensity score is modeled using logistic regression, 
#' where treatment status is regressed against a subset of non-treatment user-specified covariates. For confidence intervals, the user has the option to specify whether they wish to have the asymptotic interval, 
#' the basic (or empirical) bootstrap, or the percentile bootstrap.
#' 
#'
#' @param trt  A vector of binary treatments
#' @param xout A matrix, data frame, or vector of non-treatment covariates 
#' @param var1 The numbers or names of columns to be included in outcome model for treatment 1
#' @param var0 The numbers or names of columns to be included in outcome model for treatment 0
#' @param yout A vector of outcome variables 
#' @param varp The numbers or names of columns to be included in propensity score model 
#' 
#'
#' @return The return value is an object of...
#'
#' @examples
#'
#' @export


drest_ate<-function(trt,xout,var1,var0,yout,varp){
  
  #fitted models 
  
  #ensuring vectors are matrices for indexing
  xout<-as.matrix(xout)
  
  #treatment indexing
  trtind1<-which(trt==1)
  trtind0<-which(trt==0)
  
  #outcome models
  xout1<-xout[trtind1,var1]
  yout1<-yout[trtind1]
  
  xout0<-xout[trtind0,var0]
  yout0<-yout[trtind0]
  
  #outcome models
  
  out1fit<-lm(yout1~xout1)
  out0fit<-lm(yout0~xout0)
  
  #intercept 
  int<-rep(1,length(yout))
  fullxout1<-cbind(int,xout[,var1])
  fullxout0<-cbind(int,xout[,var0])
  yhatout1<-as.vector(fullxout1%*%out1fit$coefficients)
  yhatout0<-as.vector(fullxout0%*%out0fit$coefficients)
  
  #propensity score
  ps<-glm(trt~xout[,varp],family="binomial")$fitted
  
  #dr1
  dr1<-mean(trt*yout/ps-(trt-ps)/ps*yhatout1)
  #dr0
  dr0<-mean((1-trt)*yout/(1-ps)+(trt-ps)/(1-ps)*yhatout0)
  
  #delta of estimator
  drdelta<-dr1-dr0
  
  #asymtpotic variance 
  varterm1<-(trt*yout-yhatout1*(trt-ps))/ps-((1-trt)*yout+yhatout0*(trt-ps))/(1-ps)
  meansq<-sum((varterm1-drdelta)^2)
  se<-sqrt(meansq/length(yout)^2)
  
  return(c(drdelta,se))
  
}
