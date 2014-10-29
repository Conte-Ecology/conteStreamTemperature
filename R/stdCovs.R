#' @title Standardize variables using the mean and standard deviation from the analysis
#'
#' @description
#' \code{stdCovs} returns a dataframe of the standardized values of all the variables in x.
#' @param x Dataframe of covariates for prediction or validation
#' @param y Dataframe used for model fitting
#' @param var.names Character vector of covariate names common in x and y to standardize
#' @details
#' This function standardizes the list of covariates in the new dataframe (x) using the means and standard deviations from the dataframe (y) that was used to fit (calibrate) the model
#' @export
stdCovs <- function(x, y, var.names){
  xStd <- as.data.frame(matrix(NA, dim(x)[1], length(var.names)))
  names(xStd) <- var.names
  for(i in 1:length(var.names)){
    xStd[ , var.names[i]] <- (x[ , var.names[i]] - mean(y[ , var.names[i]], na.rm=T)) / sd(y[ , var.names[i]], na.rm=T)
  }
  return(xStd)
}