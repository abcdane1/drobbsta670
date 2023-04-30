#' Linear Regression
#'
#' Provides a Summary Output of Regression Model
#'
#' @param y dependent variable vector "y"
#' @param x independent variable matrix "x" 
#'
#' @return lm summary
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
