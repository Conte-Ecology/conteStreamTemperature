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
  data$rowNum <- 1:dim(data)[1]
  
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
  data$rowNum <- 1:dim(data)[1]
  evalRows <- data %>%
    group_by(deployID) %>%
    filter(date != min(date) & !is.na(temp)) %>%
    select(rowNum)
  
  return( evalRows$rowNum ) # this can be a list or 1 dataframe with different columns. can't be df - diff # of rows
}


#' @title not in: 
#'
#' @description
#' \code{!in} antii %in%
#'
#' @param x
#' @param table
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
"%!in%" <- function(x, table) match(x, table, nomatch = 0) == 0


#' @title prepDataWrapper: Wrapper to prepare data for analysis or predictions
#'
#' @description
#' \code{prepDataWrapper} Wrapper to prepare data for analysis or predictions
#'
#' @param var.names Character vector naming the variables to standardize for analysis
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
prepDataWrapper <- function(data.fit = NULL, var.names, dataInDir, dataOutDir, predict.daymet = FALSE, validate = FALSE, validateFrac = NULL, filter.area = NULL) {
  
  tempData <- readStreamTempData(timeSeries=TRUE, covariates=TRUE, dataSourceList=dataSource, fieldListTS=fields, fieldListCD='ALL', directory=dataInDir)
  
  springFallBPs$site <- as.character(springFallBPs$site)
  
  if(predict.daymet) {
    ########## How to add BP for years without data and clip data to the sync period ??? #######
    # Join with break points
    covariateDataBP <- left_join(covariateData, springFallBPs, by=c('site', 'year'))
    # rm(covariateData)
    
    # temp hack
    climateData$site <- as.character(climateData$site)
    tempDataBP <- left_join(climateData, covariateDataBP, by=c('site'))
    
    # Clip to syncronized season
    # tempFullSync <- filter(tempDataBP, dOY >= finalSpringBP & dOY <= finalFallBP)
    
    # temp hack
    tempDataSync <- filter(tempDataBP, dOY >= 50 & dOY <= 350)
    tempDataSync$Latitude <- tempDataSync$Latitude.x
    tempDataSync$Longitude <- tempDataSync$Longitude.x
    ##################
  } else {
    # Join with break points
    tempDataBP <- left_join(tempData, springFallBPs, by=c('site', 'year'))
    rm(tempData) # save some memory
    
    # Clip to syncronized season
    tempDataSync <- filter(tempDataBP, dOY >= finalSpringBP & dOY <= finalFallBP)
  }
  
  # Order by group and date
  tempDataSync <- tempDataSync[order(tempDataSync$site,tempDataSync$year,tempDataSync$dOY),]
  
  # For checking the order of tempDataSync
  tempDataSync$count <- 1:length(tempDataSync$year)
  
  tempDataSync <- tempDataSync[order(tempDataSync$count),] # just to make sure tempDataSync is ordered for the slide function
  
  # airTemp
  tempDataSync <- slide(tempDataSync, Var = "airTemp", GroupVar = "site", slideBy = -1, NewVar='airTempLagged1')
  tempDataSync <- slide(tempDataSync, Var = "airTemp", GroupVar = "site", slideBy = -2, NewVar='airTempLagged2')
  
  # prcp
  tempDataSync <- slide(tempDataSync, Var = "prcp", GroupVar = "site", slideBy = -1, NewVar='prcpLagged1')
  tempDataSync <- slide(tempDataSync, Var = "prcp", GroupVar = "site", slideBy = -2, NewVar='prcpLagged2')
  tempDataSync <- slide(tempDataSync, Var = "prcp", GroupVar = "site", slideBy = -3, NewVar='prcpLagged3')
  
  # Make dataframe with just variables for modeling and order before standardizing
  tempDataSync <- tempDataSync[ , c("agency", "date", "AgencyID", "year", "site", "date", "finalSpringBP", "finalFallBP", "FEATUREID", "HUC4", "HUC8", "HUC12", "temp", "Latitude", "Longitude", "airTemp", "airTempLagged1", "airTempLagged2", "prcp", "prcpLagged1", "prcpLagged2", "prcpLagged3", "dOY", "Forest", "Herbacious", "Agriculture", "Developed", "TotDASqKM", "ReachElevationM", "ImpoundmentsAllSqKM", "HydrologicGroupAB", "SurficialCoarseC", "CONUSWetland", "ReachSlopePCNT", "srad", "dayl", "swe")] #  
  
  if(class(filter.area) == "numeric") tempDataSync <- filter(tempDataSync, filter = TotDASqKM <= filter.area)
  
  ### Separate data for fitting (training) and validation
  
  #If validating:
  if(validate) {
    n.fit <- floor(length(unique(tempDataSync$site)) * (1 - validateFrac))
    
    set.seed(2346)
    site.fit <- sample(unique(tempDataSync$site), n.fit, replace = FALSE) # select sites to hold back for testing 
    tempDataSyncValid <- subset(tempDataSync, !site %in% site.fit) # data for validation
    tempDataSync <- subset(tempDataSync, site %in% site.fit)    # data for model fitting (calibration)
    
    tempDataSyncValidS <- stdCovs(x = tempDataSyncValid, y = tempDataSync, var.names = var.names)
    
    tempDataSyncValidS <- indexDeployments(tempDataSyncValidS, regional = TRUE)
    firstObsRowsValid <- createFirstRows(tempDataSyncValidS)
    evalRowsValid <-createEvalRows(tempDataSyncValidS)
    
  } else {
    tempDataSyncValid <- NULL
  }
  
  # Standardize for Analysis
  
  tempDataSyncS <- stdFitCovs(x = tempDataSync, var.names = var.names)
  
  ### Add unique deployment column and create vector to loop through for each unique site-deployment
  
  # Data for fitting
  tempDataSyncS <- indexDeployments(tempDataSyncS, regional = TRUE)
  firstObsRows <- createFirstRows(tempDataSyncS)
  evalRows <-createEvalRows(tempDataSyncS)
  
  if(validate) {
    save(tempDataSync, tempDataSyncS, tempDataSyncValid, tempDataSyncValidS, firstObsRows, evalRows, firstObsRowsValid, evalRowsValid, file = paste0(dataOutDir, 'tempDataSync.RData'))
  } else {
    save(tempDataSync, tempDataSyncS, firstObsRows, evalRows, file = paste0(dataOutDir, 'tempDataSync.RData'))
  }
}