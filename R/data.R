#' Instance of Simulation Data for Doubly Robust Estimator of ATE
#'
#' This is simply an example of one simulated data set as per Scenarios 1-3 as described in the () paper.\cr
#' `sample<-1000`
#' `x1<-rnorm(sample);x2<-rnorm(sample);probx3<-runif(sample,0,1)`\cr
#' `x3<-ifelse(probx3<=.3,1,0); x4<-rnorm(sample);lcpssim<-1.5+x1-2*x2+x3`\cr
#' `psvsim<-1/(1+exp(-lcpssim)); rsim<-runif(sample,0,1)`\cr
#' #treatment, outcome, covaraiates \cr
#' `trtsim<-ifelse(psvsim+rsim<.91,1,0); outsim<-x1+x3+2*x4; xmatsim<-cbind(x1,x2,x3)` \cr 
#' #covariates to be used in outcome and ps models \cr 
#' `s2500df<-data.frame(trtsim,xmatsim,outsim)`
#'
#' @docType data
#'
#' @usage data(s2500df)
#'
#' @format \code{"iris"} is a data frame with 150 cases (rows) and 5 variables
#'  (columns) named Sepal.Length, Sepal.Width, Petal.Length, Petal.Width, and Species.
#'
#' @keywords datasets
#'
#' @references tbd
#'
#' @source \href{https://phenome.jax.org/projects/Moore1b}{QTL Archive}
#' 
#' #' @examples
#' data(s2500df)
"s2500df"

