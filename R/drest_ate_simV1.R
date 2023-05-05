#' Doubly Robust Estimation of ATE: Simulation
#'
#' This function provides a simulation study to assess coverage of confidence intervals for the doubly robust estimator of ATE under various sample sizes 
#' and single model misspecification options (discussed below). The confidence intervals assessed are the asymptotic confidence interval and 
#' non-parametric bootstrap confidence intervals: basic/empirical, percentile, and BCa. The first models with the BCa interval excluded reproduce the results in (paper)
#'  
#'
#' @param model   A value 1-7 representing the choice of misspecification model. The models are discussed in details.
#' @param sample  A positive integer sample size (greater than 100 recommended) for the simulated data.
#' @param iterations The number of times that the simulation is run. Default is 1000.
#' @param level A value between 0 and 1 which gives confidence level of the confidence interval - default is .95. 
#' @param boot If TRUE return full bootstrap table for simulated data. 
#' @param B A postive integer numbers of bootstrap samples - default is 1000
#'
#' @return The return value table of the form 
#' 
#' @references Reference
#'
#' @export

drest_ate_simV1<-function(model,sample,iterations=1000,level=.95,boot=FALSE,B=1000){
  set.seed(010590)
  if(model==1){
    #model 1 (both correct)
    xsimm<-c(1,3)
    pssim<-c(1,3)}
  #model 2 (ps misspecified) 
  if (model==2|model==4){
    xsimm<-c(1,3)
    pssim<-1}
  #model 3 (outcome misspecified)
  if (model==3|model==5){
    xsimm<-1
    pssim<-c(1,3)}
  
  if(model==6){
    xsimm<-c(1,3,5)
    pssim<-4
  }
  
  if(model==7){
    xsimm<-4
    pssim<-c(1,3,5)
  }
  
  
  drest<-rep(0,iterations)
  acm<-rep(0,iterations)
  sd<-rep(0,iterations)
  cov_as<-rep(0,iterations)
  cov_bas<-rep(0,iterations)
  cov_perc<-rep(0,iterations)
  cov_bca<-rep(0,iterations)
  
  if (model==1|model==2|model==3){
    for (i in 1:iterations){
      x1<-rnorm(sample)
      x2<-rnorm(sample)
      probx3<-runif(sample,0,1)
      x3<-ifelse(probx3<=.3,1,0)
      x4<-rnorm(sample) #error 
      lcpssim<-1.5+x1-2*x2+x3
      psvsim<-1/(1+exp(-lcpssim))
      rsim<-runif(sample,0,1)
      trtsim<-ifelse(psvsim+rsim<.91,1,0)
      outsim<-x1+x3+2*x4
      xmatsim<-cbind(x1,x2,x3)
      estim_as<-drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,level=level)
      drest[i]<-estim_as[1]
      acm[i]<-estim_as[2]
      cov_as[i]<-ifelse(0>=estim_as[3] & 0<=estim_as[4],1,0)
      
      if(boot==TRUE){
        estim_bas<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="basic",level=level,B=B))
        estim_perc<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="percentile",level=level,B=B))
        estim_bca<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="bca",level=level,B=B))
        cov_bas[i]<-ifelse(0>=estim_bas[3] & 0<=estim_bas[4],1,0)
        cov_perc[i]<-ifelse(0>=estim_perc[3] & 0<=estim_perc[4],1,0)
        cov_bca[i]<-ifelse(0>=estim_bca[3] & 0<=estim_bca[4],1,0)}
    }}
  
  if (model==4|model==5){
    for (i in 1:iterations){
      x1<-rexp(sample)
      x2<-rnorm(sample)
      probx3<-runif(sample,0,1)
      x3<-ifelse(probx3<=.3,1,0)
      x4<-rnorm(sample) #error 
      lcpssim<-1.5+x1-2*x2+x3
      psvsim<-1/(1+exp(-lcpssim))
      rsim<-runif(sample,0,1)
      trtsim<-ifelse(psvsim+rsim<.91,1,0)
      outsim<-x1+x3+2*x4
      xmatsim<-cbind(x1,x2,x3)
      estim_as<-drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,level=level)
      drest[i]<-estim_as[1]
      acm[i]<-estim_as[2]
      cov_as[i]<-ifelse(0>=estim_as[3] & 0<=estim_as[4],1,0)
      
      if(boot==TRUE){
        estim_bas<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="basic",level=level,B=B))
        estim_perc<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="percentile",level=level,B=B))
        estim_bca<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="bca",level=level,B=B))
        cov_bas[i]<-ifelse(0>=estim_bas[3] & 0<=estim_bas[4],1,0)
        cov_perc[i]<-ifelse(0>=estim_perc[3] & 0<=estim_perc[4],1,0)
        cov_bca[i]<-ifelse(0>=estim_bca[3] & 0<=estim_bca[4],1,0)}
      
    }}
  
  
  #add variable, with x5 only, vs sin(x_5)
  if (model==6|model==7){
    for (i in 1:iterations){
      x1<-rnorm(sample)
      x2<-rnorm(sample)
      probx3<-runif(sample,0,1)
      x3<-ifelse(probx3<=.3,1,0)
      x4<-rnorm(sample) #error 
      #added term
      x5<-rnorm(sample)
      lcpssim<-1.5+x1-2*x2+x3
      psvsim<-1/(1+exp(-lcpssim))
      rsim<-runif(sample,0,1)
      trtsim<-ifelse(psvsim+rsim<.91,1,0)
      outsim<-x1+x3+2*x4-1.5*sin(x5)
      xmatsim<-cbind(x1,x2,x3,x5,sin(x5))
      estim_as<-drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,level=level)
      drest[i]<-estim_as[1]
      acm[i]<-estim_as[2]
      cov_as[i]<-ifelse(0>=estim_as[3] & 0<=estim_as[4],1,0)
      
      if(boot==TRUE){
        estim_bas<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="basic",level=level,B=B))
        estim_perc<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="percentile",level=level,B=B))
        estim_bca<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="bca",level=level,B=B))
        cov_bas[i]<-ifelse(0>=estim_bas[3] & 0<=estim_bas[4],1,0)
        cov_perc[i]<-ifelse(0>=estim_perc[3] & 0<=estim_perc[4],1,0)
        cov_bca[i]<-ifelse(0>=estim_bca[3] & 0<=estim_bca[4],1,0)}
      
    }}
  
  bias<-mean(drest,na.rm=TRUE)
  seacm<-mean(acm,na.rm=TRUE)
  sdout<-sd(drest,na.rm=TRUE)
  sesd<-seacm/sdout
  pcov_as<-100*mean(cov_as,na.rm=TRUE)
  cov_as<-cov_as[!is.na(cov_as)]
  ci_as<-100*prop.test(sum(cov_as),length(cov_as),correct=FALSE)$conf.int
  
  if(boot==TRUE){
    pcov_bas<-100*mean(cov_bas,na.rm=TRUE)
    pcov_perc<-100*mean(cov_perc,na.rm=TRUE)
    pcov_bca<-100*mean(cov_bca,na.rm=TRUE)
    
    cov_bas<-cov_bas[!is.na(cov_bas)]
    cov_perc<-cov_perc[!is.na(cov_perc)]
    cov_bca<-cov_bca[!is.na(cov_bca)]
    
    ci_bas<-100*prop.test(sum(cov_bas),length(cov_bas),correct=FALSE)$conf.int
    ci_perc<-100*prop.test(sum(cov_perc),length(cov_perc),correct=FALSE)$conf.int
    ci_bca<-100*prop.test(sum(cov_bca),length(cov_bca),correct=FALSE)$conf.int}
  
  #}
  #return(c(model,sample,bias,seacm,sdout,sesd,pcov_as,pcov_bas,pcov_perc))
  names<-c("Model","SampleSize","Bias","SE_ACM","SD","SE_ACM/SD","CovSE_ACM","CovSE_ACMLB","CovSE_ACMUB",
           "CovBas","CovBasLB","CovBasUB","CovPerc","CovPercLB","CovPercUB","CovBCa","CovBCaLB","CovBCaUB")
  results<-c(model,sample,bias,seacm,sdout,sesd,pcov_as,ci_as)
  names(results)<-names[1:9]
  if(boot==TRUE){
    results<-c(model,sample,bias,seacm,sdout,sesd,pcov_as,ci_as,pcov_bas,ci_bas,
               pcov_perc,ci_perc,pcov_bca,ci_bca)
    names(results)<-names}
  
  return(round(results,3))
}