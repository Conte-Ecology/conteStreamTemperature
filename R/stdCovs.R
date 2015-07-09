#' @title Standardize variables using the mean and standard deviation from the analysis
#'
#' @description
#' \code{stdCovs} returns a dataframe of the standardized values of all the variables in x.
#' @param x Dataframe of covariates for prediction or validation
#' @param y Dataframe of means and standard deviations
#' @param var.names Character vector of covariate names common in x and y to standardize
#' @details
#' This function standardizes the list of covariates in the new dataframe (x) using the means and standard deviations from the dataframe (y) that was used to fit (calibrate) the model
#' @export
stdCovs <- function(x, y, var.names){
  for(i in 1:length(var.names)){
    x[ , var.names[i]] <- (x[ , var.names[i]] - y[which(var.names == var.names[i]), "means"]) / y[which(var.names == var.names[i]), "stdevs"]
  }
  return(x)
}


#' @title Standardize covariates for model fitting
#'
#' @description
#' \code{stdFitCovs} returns a dataframe of the standardized values of all the variables in x.
#' @param x Dataframe of covariates for prediction or validation
#' @param var.names Character vector of covariate names common in x and y to standardize
#' @details
#' This function standardizes the list of covariates in the new dataframe (x) using the means and standard deviations from the dataframe (y) that was used to fit (calibrate) the model
#' @export
stdFitCovs <- function(x, var.names){
  for(i in 1:length(var.names)){
    x[ , var.names[i]] <- (x[ , var.names[i]] - mean(x[ , var.names[i]], na.rm = T)) / sd(x[ , var.names[i]], na.rm = T)
  }
  return(x)
}





