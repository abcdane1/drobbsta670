#' Doubly Robust Estimation of ATE: Simulation
#'
#' This function provides a simulation study to assess coverage of confidence intervals for the doubly robust estimator of ATE under various sample sizes 
#' and single model misspecification options (discussed below). The confidence intervals assessed are the asymptotic confidence interval and 
#' non-parametric bootstrap confidence intervals: basic/empirical, percentile, and BCa. The first models with the BCa interval excluded reproduce the results in (paper)
#' This function allows for parallelization, where you can specify the number of parallel tasks and number of iterations for each task. 
#'
#' @param model   A value 1-7 representing the choice of misspecification model. The models are discussed in details.
#' @param sample  A positive integer sample size (greater than 100 recommended) for the simulated data.
#' @param iterations The number of times that each parallel operation is run. Default is 100. 
#' @param rounds The number of parallel operations to be performed. If `round=1` as default, there will be no parallelization (not recommended, especially for `boot=TRUE`). 
#' The total number of simulations will be `iterations` times `rounds`.
#' @param level A value between 0 and 1 which gives confidence level of the confidence interval - default is .95. 
#' @param boot If TRUE return full bootstrap table for simulated data. 
#' @param B A positive integer numbers of bootstrap samples - default is 1000. 
#' @param nc Number of cores for parallelization - default is 4. If `nc=1`, there will be no parallelization.
#'
#' @return The return value table of the form 
#' 
#' @importFrom stats lm glm qnorm pnorm prop.test quantile rexp runif sd rnorm
#' @importFrom parallel mcmapply
#' 
#' @references Reference
#' 
#' 
#' @export
#' 
#' 
drest_ate_sim<-function(model,sample,iterations=100,rounds=1,level=.95,boot=FALSE,B=1000,nc=4){

#subfunction to choose model misspecification
drest_ate_simsub<-function(model,sample,iterations,rounds,level,boot,B){
  if(model==1){
    #model 1 (both correct)
    xsimm<-c(1,3)
    pssim<-c(1,3)}
  #model 2 (ps misspecified) 
  if (model==2|model==4){
    xsimm<-c(1,3)
    pssim<-1}
  #model 3,5 (outcome misspecified, X1 now exp)
  if (model==3|model==5){
    xsimm<-1
    pssim<-c(1,3)}
  #model 6 (ps misspecified, more missing and non-linear transform)
  if(model==6){
    xsimm<-c(1,3,5)
    pssim<-4
  }
  
  #model 6 (outcome misspecified, more missing and non-linear transform)
  if(model==7){
    xsimm<-4
    pssim<-c(1,3,5)
  }
  
  #empty vectors to fill
  drest<-rep(0,iterations)
  acm<-rep(0,iterations)
  sd<-rep(0,iterations)
  cov_as<-rep(0,iterations)
  cov_bas<-rep(0,iterations)
  cov_perc<-rep(0,iterations)
  cov_bca<-rep(0,iterations)
  
  #data generation process for first three models
  if (model==1|model==2|model==3){
    for (i in 1:iterations){
      #x1, x2, x4 norm
      x1<-rnorm(sample)
      x2<-rnorm(sample)
      #x3 bin(.3)
      probx3<-runif(sample,0,1)
      x3<-ifelse(probx3<=.3,1,0)
      x4<-rnorm(sample) 
      #true ps model
      lcpssim<-1.5+x1-2*x2+x3
      psvsim<-1/(1+exp(-lcpssim))
      rsim<-runif(sample,0,1)
      trtsim<-ifelse(psvsim+rsim<.91,1,0)
      #true outcome model
      outsim<-x1+x3+2*x4
      #data
      xmatsim<-cbind(x1,x2,x3)
      #estimates as in drest_ate
      estim_as<-drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,level=level)
      #doubly robust
      drest[i]<-estim_as[1]
      #asymptotic variance
      acm[i]<-estim_as[2]
      #asymptotic variance, coverage of ci
      cov_as[i]<-ifelse(0>=estim_as[3] & 0<=estim_as[4],1,0)
      
      #adds basic, percentile, and bootstrap intervals
      if(boot==TRUE){
        #estimation of intervals
        #try must be specified since some calulations are difficult
        estim_bas<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="basic",level=level,B=B))
        estim_perc<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="percentile",level=level,B=B))
        estim_bca<-try(drobbsta670::drest_ate(trtsim,xmatsim,xsimm,outsim,varp=pssim,ci="bca",level=level,B=B))
        #coverage of intervals
        cov_bas[i]<-ifelse(0>=estim_bas[3] & 0<=estim_bas[4],1,0)
        cov_perc[i]<-ifelse(0>=estim_perc[3] & 0<=estim_perc[4],1,0)
        cov_bca[i]<-ifelse(0>=estim_bca[3] & 0<=estim_bca[4],1,0)}}}
  
  #same as above but changes x1 to rexp, always misspecify one 
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
        cov_bca[i]<-ifelse(0>=estim_bca[3] & 0<=estim_bca[4],1,0)}}}
  
  
  #same as above adds variable sin(x5) (so removes more) and always misspecify one
  if (model==6|model==7){
    for (i in 1:iterations){
      x1<-rnorm(sample)
      x2<-rnorm(sample)
      probx3<-runif(sample,0,1)
      x3<-ifelse(probx3<=.3,1,0)
      x4<-rnorm(sample) #error 
      #added term
      x5<-rnorm(sample)
      lcpssim<-1.5+x1-2*x2+x3-1.5*sin(x5)
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
        cov_bca[i]<-ifelse(0>=estim_bca[3] & 0<=estim_bca[4],1,0)}}}
  
  #if boot=T, adds more information
  if(boot==TRUE){
    results<-list(drest,acm,cov_as,cov_bas,cov_perc,cov_bca)
    return(results)}
  
  #standard information (default)
  results<-list(drest,acm,cov_as)
  return(results)}

#parallel operations 
RNGkind("L'Ecuyer-CMRG")
set.seed(010590) #fixed output if fixed number of cores  
l<-parallel::mcmapply(drest_ate_simsub,rounds=1:rounds,MoreArgs=list(model=model,sample=sample,iterations=iterations,
                                                                           boot=boot,level=level,B=B),mc.cores=nc)

#final listings
drestfin<-unlist(l[1,])
acmfin<-unlist(l[2,])
covfin_as<-unlist(l[3,])

#final output without boot
bias<-mean(drestfin,na.rm=TRUE)
seacm<-mean(acmfin,na.rm=TRUE)
sdout<-sd(drestfin,na.rm=TRUE)
sesd<-seacm/sdout
pcov_as<-100*mean(covfin_as,na.rm=TRUE)
#removes potential na values
covfin_as<-covfin_as[!is.na(covfin_as)]
#wilson score method ci without continuity correction for coverage ci
ci_as<-100*prop.test(sum(covfin_as),length(covfin_as),correct=FALSE)$conf.int

#added to output without boot
#repeats interval for coverage and coverage ci (na removed)
if(boot==TRUE){
  covfin_bas<-unlist(l[4,])
  covfin_perc<-unlist(l[5,])
  covfin_bca<-unlist(l[6,])
  
  pcov_bas<-100*mean(covfin_bas,na.rm=TRUE)
  pcov_perc<-100*mean(covfin_perc,na.rm=TRUE)
  pcov_bca<-100*mean(covfin_bca,na.rm=TRUE)
  
  cov_bas<-covfin_bas[!is.na(covfin_bas)]
  cov_perc<-covfin_perc[!is.na(covfin_perc)]
  cov_bca<-covfin_bca[!is.na(covfin_bca)]
  
  ci_bas<-100*prop.test(sum(cov_bas),length(cov_bas),correct=FALSE)$conf.int
  ci_perc<-100*prop.test(sum(cov_perc),length(cov_perc),correct=FALSE)$conf.int
  ci_bca<-100*prop.test(sum(cov_bca),length(cov_bca),correct=FALSE)$conf.int}

#names 
names<-c("Model","SampleSize","Bias","SE_ACM","SD","SE_ACM/SD","CovSE_ACM","CovSE_ACMLB","CovSE_ACMUB",
         "CovBas","CovBasLB","CovBasUB","CovPerc","CovPercLB","CovPercUB","CovBCa","CovBCaLB","CovBCaUB")

#results and names without boot
results<-c(model,sample,bias,seacm,sdout,sesd,pcov_as,ci_as)
names(results)<-names[1:9]

#results and names with boot
if(boot==TRUE){
  results<-c(model,sample,bias,seacm,sdout,sesd,pcov_as,ci_as,pcov_bas,ci_bas,
             pcov_perc,ci_perc,pcov_bca,ci_bca)
  names(results)<-names}

return(round(results,3))
}
