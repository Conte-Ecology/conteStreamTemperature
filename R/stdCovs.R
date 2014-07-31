#' @title Standardize variables using the mean and standard deviation from the analysis
#'
#' @description
#' \code{stdCovs} returns the standardized values of all the variables in x.
#'
#' @details
#' x: data frame of covariates for prediction or validation
#' y: original data frame of covariates used in the analysis
stdCovs <- function(x, y, varNames){
  xStd <- as.data.frame(matrix(NA, dim(x)[1], length(varNames)))
  names(xStd) <- varNames
  for(i in 1:length(varNames)){
    xStd[ , varNames[i]] <- (x[ , varNames[i]] - mean(y[ , varNames[i]], na.rm=T)) / sd(y[ , varNames[i]], na.rm=T)
  }
  return(xStd)
}