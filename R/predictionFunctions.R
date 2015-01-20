#' @title predictTemp: Predict conditional daily temperatures for observed or unobserved sites
#'
#' @description
#' \code{predictTemp} Predict daily stream temperatures
#'
#' @param data Data frame of covariates for prediction
#' @param data.fit Data frame used for fitting (calibrating) the JAGS model
#' @param cov.list List of covariates used in the model
#' @param coef.list List of coefficient values estimated from the model
#' @return Numeric vector of predicted daily stream temperatures
#' @details
#' The predictions are all conditional on site, HUC8, and year when the data is available and predicts the mean values when any component was not used in the fitting (not in the calibration dataset)
#' @examples
#' 
#' \dontrun{
#' Predictions <- predictTemp(data = tempDataSyncValidS, data.fit = tempDataSyncS, cov.list = cov.list, coef.list = coef.list))
#' }
#' @export
predictTemp <- function(data, data.fit = tempDataSyncS, coef.list, cov.list) {
  
  B.site <- prepConditionalCoef(coef.list = coef.list, cov.list = cov.list, var.name = "site")
  B.huc <- prepConditionalCoef(coef.list = coef.list, cov.list = cov.list, var.name = "huc")
  B.year <- prepConditionalCoef(coef.list = coef.list, cov.list = cov.list, var.name = "year")
  B.ar1 <- prepConditionalCoef(coef.list = coef.list, cov.list = cov.list, var.name = "ar1")
  
  df <- prepPredictDF(data = data, coef.list = coef.list, cov.list = cov.list, var.name = "site")
  df <- prepPredictDF(data = df, coef.list = coef.list, cov.list = cov.list, var.name = "huc")
  df <- prepPredictDF(data = df, coef.list = coef.list, cov.list = cov.list, var.name = "year")
  df <- prepPredictDF(data = df, coef.list = coef.list, cov.list = cov.list, var.name = "ar1")
  
  
  df$trend <- NA
  df$trend <- as.vector(coef.list$B.fixed$mean %*% t(as.matrix(select(data, one_of(cov.list$fixed.ef))))) +
    rowSums(as.matrix(select(df, one_of(cov.list$site.ef))) * as.matrix(select(df, one_of(names(B.site[-1]))))) +
    rowSums(as.matrix(select(df, one_of(cov.list$huc.ef))) * as.matrix(select(df, one_of(names(B.huc[-1]))))) +
    rowSums(as.matrix(select(df, one_of(cov.list$year.ef))) * as.matrix(select(df, one_of(names(B.year[-1])))))
  
  # Add B.ar1 to predictions
  df <- mutate(df, prev.temp = c(NA, temp[(2:(nrow(data))) -1]))
  df <- mutate(df, prev.trend = c(NA, trend[(2:nrow(data)) - 1]))
  df <- mutate(df, prev.err = prev.temp - prev.trend)
  df <- mutate(df, tempPredicted = trend)
  df[which(!is.na(df$prev.err)), ]$tempPredicted <- df[which(!is.na(df$prev.err)), ]$trend + df[which(!is.na(df$prev.err)), ]$B.ar1 * df[which(!is.na(df$prev.err)), ]$prev.err
  
  return(df)
}

#' @title prepPredictDF
#'
#' @description
#' \code{prepPredictDF} Helper function to prepare dataframe for predictions
#'
#' @param data Dataframe for which predictions will be calculated
#' @param cov.list List of covariates used in the model
#' @param coef.list List of coefficient values estimated from the model
#' 
#' @return Returns Dataframe of covariates, coefficients from a fitted model, and observed temperature when available
#' @details
#' Used within the predictTemp function
#'  
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
prepPredictDF <- function(data, coef.list, cov.list, var.name) {
  B <- prepConditionalCoef(coef.list = coef.list, cov.list = cov.list, var.name = var.name)
  if(var.name == "ar1") {
    var.name <- "site"
    data[ , var.name] <- as.factor(data[ , var.name])
    B <- dplyr::select(B, site = site, B.ar1 = mean)
    df <- left_join(data, B, by = var.name)
    df[ , names(B[-1])][is.na(df[ , names(B[-1])])] <- colMeans(B[-1]) # replace NA with mean
  } else {
    B[ , var.name] <- as.character(B[ , var.name])
    data[ , var.name] <- as.character(data[ , var.name])
    df <- left_join(data, B, by = var.name) # merge so can apply/mutate by rows without a slow for loop
    for(i in 2:length(names(B))) {
      df[ , names(B[i])][is.na(df[ , names(B[i])])] <- colMeans(B[i])
    }
    #df[ , names(B[-1])][is.na(df[ , names(B[-1])])] <- colMeans(B[-1]) # replace NA with mean
  }
  return(df)
}

#' @title prepConditionalCoef
#'
#' @description
#' \code{prepConditionalCoef} Helper function to prepare dataframe for predictions
#'
#' @param cov.list List of covariates used in the model
#' @param coef.list List of coefficient values estimated from the model
#' @param var.name Name of variable for which coefficients need to be organized
#' 
#' @return Returns Dataframe of coefficients from a fitted model associated with a variable (var.name)
#' @details
#' Used within the prepPredictDF function
#'  
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
prepConditionalCoef <- function(coef.list, cov.list, var.name) {
  if(var.name == "ar1") {
    B <- coef.list[[paste0("B.", var.name)]]
  } else {
    f <- paste0(var.name, " ~ coef")
    B <- dcast(coef.list[[paste0("B.", var.name)]], formula = as.formula(f), value.var = "mean") # conver long to wide
    B <- dplyr::select(B, one_of(c(var.name, cov.list[[paste0(var.name, ".ef")]]))) # recorder to match for matrix multiplcation
    names(B) <- c(names(B[1]), paste0(names(B[-1]), ".B.", var.name)) # rename so can differentiate coeficients from variables
  }
  
  return(B)
}


#' @title Plot observed and predicted values of stream temperature
#'
#' @description
#' \code{plotPredicted} plot observed and predicted values of stream temperature along with daily air temperature
#'
#' @param observed Dataframe of observed stream temperatures with at least columns "temp", "sites", and "years"
#' @param predicted Dataframe of predicted stream temperatures with at least columns "tempPredicted", "sites", "years", and "airTemp"
#' @param siteList Optional character vector of site names for predictions. Default "ALL" will predict to all sites in the predicted dataframe
#' @param yearList Optional character vector of years for predictions. Default "ALL" will predict for all years in the predicted dataframe
#' @param dir Directory where the files will be saved. Defaults to the working directory
#' @return Saves a png file for every site in the siteList to the dir with all years in the yearList 
#' @details
#' blah, blah, blah
#' @examples
#' 
#' \dontrun{
#' plotPredict(observed = tempDataSync, predicted = tempFull, siteList = "ALL", yearList = "ALL", dir = paste0(dataLocalDir,'/', 'plots/fullRecord/'))
#' }
#' @export
plotPredict <- function(observed, predicted, siteList = "ALL", yearList = "ALL", dir = getwd(), display = FALSE){ # add option to not include predicted or make similar function that makes observation plots
  if(siteList == "ALL" & display == TRUE) stop("If plots wanted for all sites change display to FALSE and set output directory for saved plots")
  if(siteList == "ALL"){
    sites <- unique(as.character(predicted$site))
  } else {
    sites <- siteList
  }
  if(yearList == "ALL"){
    years <- unique(as.character(predicted$year))
  } else {
    years <- yearList # add check if character vector, change if numeric, check if match any in predicted, stop with error else
  }
  
  ###### Need to convert back to original scale or join with original DF
  #predicted.origin.scale <- left_join(observed, predicted[ , c("site", "date", "tempPredicted")], by = c("site", "date"))
  
  for(i in 1:length(unique(sites))){
    dataSite <- dplyr::filter(predicted, filter = site == sites[i] & year %in% years)
    #dataSiteObs <- dplyr::filter(observed, filter = site == sites[i] & year %in% years)
    foo <- ggplot(dataSite, aes(dOY, tempPredicted)) + 
      coord_cartesian(xlim = c(100, 300), ylim = c(0, 35)) + 
      geom_point(data=dataSite, aes(dOY, temp), colour='blue') +
      geom_point(colour = 'red') + 
      geom_line(colour = 'red') + 
      geom_point(aes(dOY, airTemp), colour = 'black') + 
      ggtitle(dataSite$site[i]) + 
      facet_wrap(~year) + 
      xlab(label = 'Day of the year') + ylab('Temperature (C)') + 
      theme(axis.text.x = element_text(angle = 45), 
            axis.text.y = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'),
            legend.key = element_rect(colour = "grey80"), 
            panel.background = element_rect(fill = "white", colour = NA), 
            panel.border = element_rect(fill = NA, colour = "grey50"), 
            panel.grid.major = element_line(colour = "grey", size = 0.2), 
            panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
            strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2))
    if(display) {
      print(foo)
    } else {
      ggsave(filename=paste0(dir, dataSite$site[i], '.png'), plot=foo, dpi=300, width=12, height=8, units='in' )
    }
  } # surprisingly fast but wouldn't do for all catchments
}


#' @title predCubic
#'
#' @description
#' \code{predCubic} Get predictions of the cubic random temperature effect across years
#'
#' @param v is a data frame of values for each paramter in the mcmc temperature model. Required columns are 'Parameter' and 'mean'
#' 
#' @return Returns predCubic a data frame with the predictions by day of year and bY which contains the paramter estimates and predicted temperatures from the random effect component of the temperature model
#' @details
#' This function extracts the among-year temperature random effects parameters from the model output, does predictions of the effects by day of year and prevodes predictions at 4 day of year values (z-score values)
#'  
#' @examples
#' 
#' \dontrun{
#' preds <- predCubic( values )
#' }
#' @export
predCubic <- function( v ){
  
  #find B.year values and do predictions
  
  # rows that start with B.year
  bYearRows <- grep("^B.year",x=v$Parameter) # the ^ excludes anything before 'B.year
  bYear <- v[bYearRows,] 
  
  # Need to pull out the 'year' and 'term' variables from B.year[year,term]
  #location in text sequence
  firstBracket <- regexpr("\\[",bYear$Parameter)
  comma <- regexpr("\\,",bYear$Parameter)
  lastBracket <- regexpr("\\]",bYear$Parameter)
  
  # pull out values
  bYear$year <- substring(bYear$Parameter,firstBracket+1,comma-1)
  bYear$p <- substring(bYear$Parameter,comma+1,lastBracket-1)
  
  # cast bYear to wide format
  bYear$m <- bYear$mean #can't use 'mean' in dcast
  bY <- dcast(bYear, year~p,value=m)
  bY$year <- as.numeric(bY$year)
  names(bY) <- c('year','int','q1','q2','q3')
  
  # set up data frame for predictions for each year, dOY combo
  predCubic1 <- data.frame(dOY=rep(unique(tempDataSyncS$dOY),dim(bY)[1]),
                           year=rep(1:nrow(bY), each=length(unique(tempDataSyncS$dOY)))
  )
  
  predCubic <- left_join(predCubic1,bY)
  predCubic$pred <- predCubic$int + predCubic$q1 * predCubic$dOY + 
    predCubic$q2 * predCubic$dOY^2 +
    predCubic$q3 * predCubic$dOY^3
  
  
  # get maximum for each year
  # need to multiply the q's by the derivative multiplier
  bY <- bY %>% group_by(year) %>% mutate( maxX=as.numeric(polyroot(c(q1,2*q2,3*q3))[1]) )
  bY <- bY %>% mutate( maxY = int + q1*maxX + q2*maxX^2 + q3*maxX^3,
                       dOYMinus2 =  int + q1*-2 + q2*(-2)^2 + q3*(-2)^3,
                       dOYMinus1 =  int + q1*-1 + q2*(-1)^2 + q3*(-1)^3,
                       dOYPlus1  =  int + q1* 1 + q2*( 1)^2 + q3*( 1)^3,
                       dOYPlus2  =  int + q1 *2 + q2*   2^2 + q3*   2^3
                     )  
  
  return(list(predCubic=predCubic,bY=bY))
   
} 



#' @title calcThresholdDays
#'
#' @description
#' \code{calcThresholdDays} Calculate derived metrics from predicted stream temperatures
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param temp.threshold Optional numeric temperature threshold value in degrees C
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' blah, blah, blah, something, something
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calcThresholdDays <- function(grouped.df, derived.df, temp.threshold, summer = FALSE) {
  var.name <- paste0("meanDays.", temp.threshold)
  if(summer) {
    meanDays <- grouped.df %>%
      mutate(month = as.numeric(format(date, "%m"))) %>%
      filter(month >= 6 & month <= 8) %>%
      filter(tempPredicted > temp.threshold) %>%
      dplyr::summarise(days = n()) %>%
      dplyr::summarise(meanDays = round(median(days, na.rm = T))) %>%
      dplyr::select(-month) %>%
      setNames(c("featureid", var.name))
  } else { 
    meanDays <- grouped.df %>%
      filter(tempPredicted > temp.threshold) %>%
      dplyr::summarise(days = n()) %>%
      dplyr::summarise(meanDays = round(median(days, na.rm = T))) %>%
      setNames(c("featureid", var.name))
  }
  derived.df <- left_join(derived.df, meanDays, by = "featureid")
  derived.df[ , var.name][is.na(derived.df[ , var.name])] <- 0
  rm(meanDays)
  
  return(derived.df)
}

#' @title calcYearsMaxTemp
#'
#' @description
#' \code{calcYearsMaxTemp} Calculate derived metrics from predicted stream temperatures
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param temp.threshold Optional numeric temperature threshold value in degrees C
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' blah, blah, blah, something, something
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calcYearsMaxTemp <- function(grouped.df, derived.df, temp.threshold, summer = FALSE) {
  library(lazyeval)
  var.name1 <- paste0("yearsMaxTemp.", temp.threshold)
  var.name2 <- paste0("freqMax.", temp.threshold)
  dots <- list(~n)
  if(summer) {
    yearsMaxTemp <- grouped.df %>%
      mutate(month = as.numeric(format(date, "%m"))) %>%
      filter(month >= 6 & month <= 8) %>%
      dplyr::summarise(maxTemp = max(tempPredicted, na.rm = T)) %>%
      filter(maxTemp > temp.threshold) %>%
      dplyr::summarise(yearsMaxTemp = n()) %>%
      mutate(freqMax = yearsMaxTemp / length(unique(grouped.df$year))) %>%
      dplyr::select(-month)
      setNames(c("featureid", var.name1, var.name2))
  }
  yearsMaxTemp <- grouped.df %>%
    dplyr::summarise(maxTemp = max(tempPredicted, na.rm = T)) %>%
    filter(maxTemp > temp.threshold) %>%
    dplyr::summarise(yearsMaxTemp = n()) %>%
    mutate(freqMax = yearsMaxTemp / length(unique(grouped.df$year))) %>%
    setNames(c("featureid", var.name1, var.name2))
  derived.df <- left_join(derived.df, yearsMaxTemp, by = "featureid")
  derived.df[ , var.name1][is.na(derived.df[ , var.name1])] <- 0
  derived.df[ , var.name2][is.na(derived.df[ , var.name2])] <- 0
  rm(yearsMaxTemp)
  
  return(derived.df)
}

#' @title calcYearsCold
#'
#' @description
#' \code{calcYearsCold} Calculate derived metrics from predicted stream temperatures
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param state Character string for specific State 
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' Different states have different classifications of cold, cool, and warm streams for regulatory purposes. This function takes the state abreviation and returns the number and frequency of years that any stream reach is expected to fall within each category.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calcYearsCold <- function(grouped.df, derived.df, temp.threshold, states) {
  library(lazyeval)
  if(class(states) != "character") stop("calcYearsCold function expect state to be a character string of 2 letter abreviations")
  avail.states <- c("CT")
  if(!all(states %in% avail.states)) stop(paste("In calcYearsCold function, water temperature classifications only available for ", avail.states))
  
  if(any(states == "CT")) {
    # Cold
    var.name1 <- paste0("years.cold.", CT)
    var.name2 <- paste0("freq.cold.", CT)
    dots <- list(~n)
    yearsCold <- grouped.df %>%
      mutate(month = as.numeric(format(date, "%m"))) %>%
      filter(month >= 6 & month <= 8) %>%
      dplyr::summarise(meanTemp = mean(tempPredicted, na.rm = T)) %>%
      filter(meanTemp < 18.29) %>%
      dplyr::summarise(yearsCold = n()) %>%
      mutate(freqCold = yearsCold / length(unique(grouped.df$year))) %>%
      dplyr::select(-month)
    setNames(c("featureid", var.name1, var.name2))
    derived.df <- left_join(derived.df, yearsCold, by = "featureid")
    derived.df[ , var.name1][is.na(derived.df[ , var.name1])] <- 0
    derived.df[ , var.name2][is.na(derived.df[ , var.name2])] <- 0
    rm(yearsCold)
    # Cool
    var.name1 <- paste0("years.cool.", CT)
    var.name2 <- paste0("freq.cool.", CT)
    dots <- list(~n)
    yearsCool <- grouped.df %>%
      mutate(month = as.numeric(format(date, "%m"))) %>%
      filter(month >= 6 & month <= 8) %>%
      dplyr::summarise(meanTemp = mean(tempPredicted, na.rm = T)) %>%
      filter(meanTemp >= 18.29 & meanTemp <= 21.70) %>%
      dplyr::summarise(yearsCool = n()) %>%
      mutate(freqCold = yearsCool / length(unique(grouped.df$year))) %>%
      dplyr::select(-month)
    setNames(c("featureid", var.name1, var.name2))
    derived.df <- left_join(derived.df, yearsCool, by = "featureid")
    derived.df[ , var.name1][is.na(derived.df[ , var.name1])] <- 0
    derived.df[ , var.name2][is.na(derived.df[ , var.name2])] <- 0
    rm(yearsCool)
    # Warm
    var.name1 <- paste0("years.warm.", CT)
    var.name2 <- paste0("freq.warm.", CT)
    dots <- list(~n)
    yearsCool <- grouped.df %>%
      mutate(month = as.numeric(format(date, "%m"))) %>%
      filter(month >= 6 & month <= 8) %>%
      dplyr::summarise(meanTemp = mean(tempPredicted, na.rm = T)) %>%
      filter(meanTemp > 21.70) %>%
      dplyr::summarise(yearsCool = n()) %>%
      mutate(freqCold = yearsCool / length(unique(grouped.df$year))) %>%
      dplyr::select(-month)
    setNames(c("featureid", var.name1, var.name2))
    derived.df <- left_join(derived.df, yearsCool, by = "featureid")
    derived.df[ , var.name1][is.na(derived.df[ , var.name1])] <- 0
    derived.df[ , var.name2][is.na(derived.df[ , var.name2])] <- 0
    rm(yearsCool)
  }
  
  return(derived.df)
}

#' @title deriveMetrics
#'
#' @description
#' \code{deriveMetrics} Calculate derived metrics from predicted stream temperatures
#'
#' @param data Dataframe that includes at least featureid, year, temp (observed), and tempPredicted
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' This function calculates the following metrics based on predicted stream temperature across years
#' 
#' * total observations (days with data) per featureid
#' * Mean maximum daily mean temperature by featureid (over years)
#' * Maximum max daily mean temperature
#' * Number of days with stream temp > 18, 20, 22 C and optional user-defined temperature
#' * Number of years with mean maximum over 18, 20, 22 C and optional user-defined temperature
#' * frequency of years with a mean max over 18, 20, 22 C and optional user-defined temperature
#' * Number of days with summer stream temp > 18.29, 21.70 C for CT DEEP
#' * Number of years with summer maximum over 18.29 & 21.7 C for CT DEEP
#' * frequency of years with a summer max over 18.29 & 21.7 C for CT DEEP
#' * Mean resistance to peak air temperature (difference between observed air and predicted stream temperatures during the hottest part of the year)
#' * Mean RMSE for each featureid
#' * Flag based on RMSE > 95%. These featureids should probably be checked for unrecorded impoundments, restrictive culverts, or large groundwater influences before making management or policy decisions
#' 
#'  
#' @examples
#' 
#' \dontrun{
#' df <- predictTemp(data = df, data.fit = tempDataSyncS, cov.list = cov.list, coef.list = coef.list)
#' derived.metrics <- deriveMetrics(data = df)
#' 
#' }
#' @export
deriveMetrics <- function(data, threshold = NULL) {
  library(dplyr)
  library(zoo)
  
  byfeatureid <- group_by(data, featureid)
  byfeatureid <- group_by(fullDataSync, featureid)
  byfeatureidYear <- group_by(byfeatureid, year, add = TRUE)
  #(maxTempfeatureid <- dplyr::dplyr::summarise(byfeatureid, max(tempPredicted, na.rm = T)))
  
  # Mean maximum daily mean temperature by featureid (over years)
  meanMaxTemp <- byfeatureidYear %>%
    dplyr::summarise(maxTempPredicted = max(tempPredicted, na.rm = T)) %>%
    dplyr::summarise(meanMaxTemp = mean(maxTempPredicted))
  derivedfeatureidMetrics <- meanMaxTemp
  rm(meanMaxTemp)
  
  # total observations (days with data) per featureid
  totalObs <- byfeatureidYear %>%
    filter(!is.na(temp)) %>%
    dplyr::summarise(Obs = n()) %>%
    dplyr::summarise(totalObs = sum(Obs)) %>%
    dplyr::select(featureid, totalObs)
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, totalObs, by = "featureid") %>%
    mutate(totalObs = ifelse(is.na(totalObs), 0, totalObs))
  derivedfeatureidMetrics <- dplyr::select(derivedfeatureidMetrics, featureid, totalObs, meanMaxTemp)
  
  # Maximum max daily mean temperature
  maxMaxTemp <- byfeatureidYear %>%
    dplyr::summarise(maxTemp = max(tempPredicted, na.rm = T)) %>%
    dplyr::summarise(maxMaxTemp = max(maxTemp))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, maxMaxTemp, by = "featureid")
  rm(maxMaxTemp)
  
  # Number of days with stream temp > threshold
  derivedfeatureidMetrics <- calcThresholdDays(byfeatureidYear, derivedfeatureidMetrics, 18)
  derivedfeatureidMetrics <- calcThresholdDays(byfeatureidYear, derivedfeatureidMetrics, 20)
  derivedfeatureidMetrics <- calcThresholdDays(byfeatureidYear, derivedfeatureidMetrics, 22)
  if(class(threshold) == "numeric") derivedfeatureidMetrics <- calcThresholdDays(byfeatureidYear, derivedfeatureidMetrics, threshold)
  
  # CT DEEP Thresholds
  #derivedfeatureidMetrics <- calcYearsCold(byfeatureidYear, derivedfeatureidMetrics, states = c("CT"))
  
  # Number and frequency of years with mean max over threshold
  derivedfeatureidMetrics <- calcYearsMaxTemp(grouped.df = byfeatureidYear, derived.df = derivedfeatureidMetrics, temp.threshold = 18)
  derivedfeatureidMetrics <- calcYearsMaxTemp(grouped.df = byfeatureidYear, derived.df = derivedfeatureidMetrics, temp.threshold = 20)
  derivedfeatureidMetrics <- calcYearsMaxTemp(grouped.df = byfeatureidYear, derived.df = derivedfeatureidMetrics, temp.threshold = 22)
  if(class(threshold) == "numeric") derivedfeatureidMetrics <- calcYearsMaxTemp(byfeatureidYear, derivedfeatureidMetrics, threshold)
  
  # Resistance to peak air temperature
  ## This probably makes the most sense during minimum flow periods but we don't have a sufficient flow model
  ## 60 or 90 days max air temp?
  # error if use standardized values rather than original scale
 # dOY.max.warm <- byfeatureidYear %>%
#    mutate(warm.90 = rollsum(x = airTemp, 90, align = "right", fill = NA))
#  dOY.max.warm <- dOY.max.warm %>%
 #   group_by(featureid, year, add = TRUE) %>%
  #  filter(warm.90 == max(warm.90, na.rm = T)) %>%
   # select(dOY)
  meanResist <- byfeatureidYear %>%
    #filter(dOY > dOY.max.warm$dOY - 90 & dOY <= dOY.max.warm$dOY) %>%
    filter(dOY >= 152 & dOY < 244) %>% # clip to summer
    mutate(absResid = abs(airTemp - tempPredicted)) %>%
    dplyr::summarise(resistance = sum(absResid, na.rm = T)) %>%
    dplyr::summarise(meanResist = mean(resistance, na.rm = T))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanResist, by = "featureid")
  rm(meanResist)
  
   # RMSE for each featureid (flag highest)
foo <- byfeatureidYear %>%
  filter(!(is.na(temp))) %>%
  mutate(error2 = (temp - tempPredicted)^2)

if(dim(foo)[1] > 0) {
  meanRMSE <- foo %>%
    dplyr::summarise(RMSE = sqrt(mean(error2))) %>%
    dplyr::summarise(meanRMSE = mean(RMSE))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanRMSE, by = "featureid")
} else {
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, rename(foo, meanRMSE = error2), by = "featureid")
}
  
  # Flag based on RMSE > 95%
# add an ifelse for whether data is present or not #######################
  derivedfeatureidMetrics <- mutate(derivedfeatureidMetrics, flag = ifelse(meanRMSE > quantile(derivedfeatureidMetrics$meanRMSE, probs = c(0.95), na.rm=TRUE), "Flag", ""))
  
  gc(verbose = FALSE)
  
  return(derivedfeatureidMetrics)
}
