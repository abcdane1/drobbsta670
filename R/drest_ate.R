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


lm_sum<-function(y,x){
  fit<-lm(y~x^2)
  return(summary(fit))
}

