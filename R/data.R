#' Instance of Simulated Data for Doubly Robust Estimator of ATE
#'
#' This is simply an example of a simulated data set as per Scenarios 1-3 as described for Table 3 of Funk et al. 2011. \cr
#' \cr
#' #data generation \cr
#' `sample<-1000`
#' `x1<-rnorm(sample); x2<-rnorm(sample); probx3<-runif(sample,0,1)` \cr
#' `x3<-ifelse(probx3<=.3,1,0); x4<-rnorm(sample)` \cr
#' #probability of treatment model \cr
#' `lcpssim<-1.5+x1-2*x2+x3; psvsim<-1/(1+exp(-lcpssim)); rsim<-runif(sample,0,1)`\cr
#' #treatment, outcome, covariates \cr
#' `trtsim<-ifelse(psvsim+rsim<.91,1,0); outsim<-x1+x3+2*x4; xmatsim<-cbind(x1,x2,x3)` \cr
#' #covariates to be used in outcome and ps models \cr
#' `s2500df<-data.frame(trtsim,xmatsim,outsim)`
#'
#' @docType data
#'
#' @usage data(s2500df)
#'
#' @format `s2500df` is a data frame with 1000 individuals (rows) and 5 variables where the first column is the treatment variable,
#' the second, third, and fourth columns are covariates, and the last column is the outcome variable.
#'
#' @keywords datasets
#'
#' @references Michele Jonsson Funk, Daniel Westreich, Chris Wiesen, Til Sturmer, M Alan Brookhart, and Marie Davidian. Doubly robust estimation of causal effects. American journal of epidemiology, 173(7):761–767, 2011.
#'
"s2500df"

