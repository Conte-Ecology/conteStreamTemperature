#' @title indexDeployments: Create an index for each unique temperature logger deployment
#'
#' @description
#' \code{indexDeployments} returns the input dataframe with a deployment ID added as a column
#'
#' @param data Dataframe of temperature data including date and unique site ID
#' @param regional Optional logical if true then also sort the dataframe by HUC8
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
indexDeployments <- function(data, regional = FALSE) {
  tbl_df(data)
  if(regional) {
  data <- arrange(data, HUC8, site, date)
  } else {
    data <- arrange(data, site, date)
  }
  
  # Within each site, if column above doesn't equal date + 1 then start new deployment ID
  # SLOW - consider rewriting or moving to C++
  deployID <- NA
  deployID[1] <- 1
  for(i in 2:nrow(data)) {
    if(data$site[i] == data$site[i-1] & 
         data$date[i] == data$date[i -1] + 1) {
      data$deployID[i] <- data$deployID[i-1]
    } else {
      data$deployID[i] <- data$deployID[i-1] + 1
    }
  } # rewrite by shifting one df down by one and comparing if equal
  return(data)
}


#' @title createDeployRows: Create rows to loop through for autoregressive function
#'
#' @description
#' \code{createDeployRows} returns a list of two dataframes each with a rowNum column for looping
#'
#' @param data Dataframe for analysis with date and a deployID row (can generate with indexDeployments function)
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
createDeployRows <- function(data) {
  data$rowNum <- 1:dim(data)[1]
  
  firstObsRows <- data %>%
    group_by(deployID) %>%
    filter(date == min(date) | is.na(date)) %>%
    select(rowNum)
  
  evalRows <- data %>%
    group_by(deployID) %>%
    filter(date != min(date) & !is.na(date)) %>%
    select(rowNum)
  
  return(list(firstObsRows, evalRows)) # this can be a list or 1 dataframe with different columns
}

