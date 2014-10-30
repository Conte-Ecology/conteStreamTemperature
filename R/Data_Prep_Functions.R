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
  data$sitef <- as.numeric(as.factor(data$site))
  data$rowNum <- 1:nrow(data)
  
  if(regional) {
  data <- arrange(data, HUC8, site, date)
  } else {
    data <- arrange(data, site, date)
  }

#  test  
#  data1=data.frame(site=rep(1:4,each=4),date=rep(1:4))
#  data=rbind(data1,data1[15:16,])
  
  data <- 
  data %>%
      mutate( siteShift = c( 1,sitef[ 1:(nrow(data)-1) ] ),
              dateShift = c( 1,date[ 1:(nrow(data)-1) ] ),
              newSite = sitef == siteShift + 1,
              newDate = date != dateShift + 1,
              newDeploy = (newSite | newDate) * 1,              
              deployID= cumsum(newDeploy) )
  
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
createFirstRows <- function(data) {
  #data$rowNum <- 1:dim(data)[1] This is created in indexDeployments
  
  firstObsRows <- data %>%
    group_by(deployID) %>%
    filter(date == min(date) | is.na(date)) %>%
    select(rowNum)
  
  return( firstObsRows$rowNum ) # this can be a list or 1 dataframe with different columns. can't be df - diff # of rows
}

#' @title createEvalRows: Create rows to loop through for autoregressive function
#'
#' @description
#' \code{createDeployRows} returns a list of two dataframes each with a rowNum column for looping
#'
#' @param data Dataframe for analysis with date and a deployID row (can generate with indexDeployments function)
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
createEvalRows <- function(data) {
  #data$rowNum <- 1:dim(data)[1]
  evalRows <- data %>%
    group_by(deployID) %>%
    filter(date != min(date) & !is.na(temp)) %>%
    select(rowNum)
  
  return( evalRows$rowNum ) # this can be a list or 1 dataframe with different columns. can't be df - diff # of rows
}



#' @title addStreamMuResid 
#'
#' @description
#' \code{addStreamMuResid} get means for stream.mu and resid from the jags output. append means to tempDataSyncS
#'
#' @param M,tempDataSyncS
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
addStreamMuResid <- function(M.wb,tempDataSyncS){
    
  getMeansStreamMuResid <- 
    M.wb %>%
    as.matrix(.) %>%
    melt(.) %>%
    group_by(Var2) %>%
    summarise( mean = mean(value) ) %>%
    mutate( rowNum = as.numeric( substring(Var2,regexpr("\\[",Var2)+1,regexpr("\\]",Var2)-1) ) ) %>%
    filter( row_number() %in% grep("stream.mu",x=Var2) | 
              row_number() %in% grep("residuals",x=Var2)) 
  
  m1 <- getMeansStreamMuResid %>%
    filter(row_number() %in% grep("residuals",x=Var2) ) %>%
    mutate( resid.wb = mean ) 
  
  tempDataSyncS$resid.wb <- NULL
  tempDataSyncS <- left_join( tempDataSyncS,m1[,c('rowNum','resid.wb')], by = 'rowNum')
  
  m2 <- getMeansStreamMuResid %>%
    filter(row_number() %in% grep("stream.mu",x=Var2) ) %>%
    mutate( pred.wb = mean ) 
  
  tempDataSyncS$pred.wb <- NULL
  tempDataSyncS <- left_join( tempDataSyncS,m2[,c('rowNum','pred.wb')], by = 'rowNum')
  
  return( tempDataSyncS )
  
}

# this could go in another file. Ben stuck it here for now

#' @title not in: 
#'
#' @description
#' \code{!in} antii %in%
#'
#' @param x,table
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0