#' @title predictTemp: Predict conditional daily temperatures for observed or unobserved sites
#'
#' @description
#' \code{predictTemp} Predict daily stream temperatures
#'
#' @param data Data frame of covariates for prediction
#' @param data.fit Data frame used for fitting (calibrating) the JAGS model
#' @param cov.list List of covariates used in the model
#' @param coef.list List of coefficient values estimated from the model
#' @param featureid_site dataframe with featureid and sites linked
#' @param validate Logical if TRUE then the trend is used for the predictions and not adjusted for autocorrelation in the error. Defaults to FALSE
#' @return Numeric vector of predicted daily stream temperatures
#' @details
#' The predictions are all conditional on site, HUC8, and year when the data is available and predicts the mean values when any component was not used in the fitting (not in the calibration dataset)
#' @examples
#' 
#' \dontrun{
#' Predictions <- predictTemp(data = tempDataSyncValidS, data.fit = tempDataSyncS, cov.list = cov.list, coef.list = coef.list, featureid_site = featureid_site))
#' }
#' @export
predictTemp <- function(fullDataSyncS, coef.list, rand_ids, Random_AR1 = FALSE) {
  
  if(!exists("fullDataSyncS$sitef")) {
    fullDataSyncS <- fullDataSyncS %>%
      left_join(rand_ids$df_site)
  }
  if(!exists("fullDataSyncS$hucf")) {
    fullDataSyncS <- fullDataSyncS %>%
      left_join(rand_ids$df_huc)
  }
  if(!exists("fullDataSyncS$yearf")) {
    fullDataSyncS <- fullDataSyncS %>%
      left_join(rand_ids$df_year)
  }
  
  ############# Predictions ##############
  #fullDataSyncS <- predictTemp(data = fullDataSyncS, coef.list = coef.list, cov.list = cov.list, featureid_site = featureid_site)
  
  fixed.ef <- as.numeric(coef.list$B.fixed$mean) # separate out the iteration or do for mean/median
  
  # add specific random effects to the dataframe
  fullDataSyncS <- left_join(fullDataSyncS, coef.list$B.site)
  fullDataSyncS <- left_join(fullDataSyncS, coef.list$B.huc) # problem with validation data, need to use the mean when huc don't match
  fullDataSyncS <- left_join(fullDataSyncS, coef.list$B.year)
  
  
  for (j in 2:length(names(coef.list$B.site))) {
    fullDataSyncS[, names(coef.list$B.site[j])][is.na(fullDataSyncS[, names(coef.list$B.site[j])])] <- colMeans(coef.list$B.site[j])
  }
  for (j in 2:length(names(coef.list$B.huc))) {
    fullDataSyncS[, names(coef.list$B.huc[j])][is.na(fullDataSyncS[, names(coef.list$B.huc[j])])] <- colMeans(coef.list$B.huc[j])
  }
  for (j in 2:length(names(coef.list$B.year))) {
    fullDataSyncS[, names(coef.list$B.year[j])][is.na(fullDataSyncS[, names(coef.list$B.year[j])])] <- colMeans(coef.list$B.year[j])
  }
  
  fullDataSyncS$fixed.ef <- as.vector(fixed.ef %*% t(as.matrix(as.data.frame(unclass(select(ungroup(fullDataSyncS), one_of(cov.list$fixed.ef)))))))
  fullDataSyncS$site.ef <- rowSums(as.matrix(select(fullDataSyncS, one_of(cov.list$site.ef))) * as.matrix(select(fullDataSyncS, starts_with("B.site"))))
  fullDataSyncS$huc.ef <- rowSums(as.matrix(select(fullDataSyncS, one_of(cov.list$huc.ef))) * as.matrix(select(fullDataSyncS, starts_with("B.huc"))))
  fullDataSyncS$year.ef <- rowSums(as.matrix(select(fullDataSyncS, one_of(cov.list$year.ef))) * as.matrix(select(fullDataSyncS, starts_with("B.year"))))
  
  # fullDataSyncS$trend <- rowSums(as.matrix(dplyr::select(fullDataSyncS, one_of(c("fixed.ef", "site.ef", "huc.ef", "year.ef")))))
  # FAILS
  
  #fullDataSyncS <- fullDataSyncS %>%
  # dplyr::mutate(trend = fixed.ef + site.ef + huc.ef + year.ef) # works fine outside foreach
  
  fullDataSyncS$trend <- fullDataSyncS$fixed.ef + fullDataSyncS$site.ef + fullDataSyncS$huc.ef + fullDataSyncS$year.ef # WORKS!!!
  
  # Add B.ar1 to predictions
  #fullDataSyncS <- group_by(fullDataSyncS, sitef)
  fullDataSyncS <- mutate(fullDataSyncS, prev.temp = c(NA, fullDataSyncS$temp[(2:(nrow(fullDataSyncS))) -1]),
                          prev.trend = c(NA, fullDataSyncS$trend[(2:nrow(fullDataSyncS)) - 1]),
                          prev.err = prev.temp - prev.trend,
                          tempPredicted = trend,
                          prev.temp = ifelse(newDeploy == 1, NA, prev.temp),
                          prev.err = ifelse(newDeploy == 1, NA, prev.err))
  
  if(Random_AR1) {
    #B.ar1.sub <- data.frame(sitef = rand_ids$df_site$sitef)
    B.ar1.sub <- coef.list$B.ar1 %>%
      dplyr::select(sitef, mean) %>%
      dplyr::rename(B.ar1 = mean)
    fullDataSyncS <- left_join(fullDataSyncS, B.ar1.sub, by = c("sitef"))
    fullDataSyncS <- fullDataSyncS %>%
      dplyr::mutate(B.ar1 = ifelse(is.na(B.ar1), mean(B.ar1.sub$B.ar1, na.rm = T), B.ar1)) %>%
      dplyr::arrange(featureid, date) 
  } else {
    fullDataSyncS$B.ar1 <- coef.list$B.ar1$mean
    fullDataSyncS <- fullDataSyncS %>% dplyr::arrange(featureid, date) 
  }
  
  
  fullDataSyncS[which(!is.na(fullDataSyncS$prev.err)), ]$tempPredicted <- fullDataSyncS[which(!is.na(fullDataSyncS$prev.err)), ]$trend + fullDataSyncS[which(!is.na(fullDataSyncS$prev.err)), ]$B.ar1 * fullDataSyncS[which(!is.na(fullDataSyncS$prev.err)), ]$prev.err
  
  # unofficial warning message
  #   mean.pred <- mean(fullDataSync$tempPredicted, na.rm = T)
  #   if(mean.pred == "NaN") {
  #     warning(paste0(i, " of ", n.loops, " loops has no predicted temperatures"))
  #   } 
  
  fullDataSyncS <- dplyr::arrange(fullDataSyncS, rowNum)
  
  fullDataSyncS <- data.frame(unclass(fullDataSyncS), stringsAsFactors = FALSE)
  
  return(fullDataSyncS)
  
}



#' @title prepData
#'
#' @description
#' \code{prepData} Helper function to prepare dataframe for predictions
#'
#' @param data Dataframe for which predictions will be calculated
#' @param cov.list List of covariates used in the model
#' @param coef.list List of coefficient values estimated from the model
#' @param featureid_site dataframe with featureid and site as 2 columns
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
prepData <- function(catches_string, springFallBPs, df_covariates_upstream, tempDataSync, featureid_lat_lon, featureid_huc8, rand_ids) {
  
  ########## pull daymet records ##########
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(drv, dbname='sheds', host='osensei.cns.umass.edu', user=options('SHEDS_USERNAME'), password=options('SHEDS_PASSWORD'))
  
  qry_daymet <- paste0("SELECT featureid, date, tmax, tmin, prcp, dayl, srad, vp, swe, (tmax + tmin) / 2.0 AS airtemp FROM daymet WHERE featureid IN (", catches_string, ") ;")
  rs <- dbSendQuery(con, statement = qry_daymet)
  climateData <- dbFetch(rs, n=-1)
  
  dbClearResult(rs)
  dbDisconnect(con)
  #dbUnloadDriver(drv)
  
  ########## Assign synchronized period ##########
  mean.spring.bp <- mean(dplyr::filter(springFallBPs, finalSpringBP != "Inf")$finalSpringBP, na.rm = T)
  mean.fall.bp <- mean(dplyr::filter(springFallBPs, finalFallBP != "Inf")$finalFallBP, na.rm = T)
  
  foo <- springFallBPs %>%
    dplyr::mutate(site = as.character(site),
                  featureid = as.integer(site),
                  finalSpringBP = ifelse(finalSpringBP == "Inf" | is.na(finalSpringBP), mean.spring.bp, finalSpringBP),
                  finalFallBP = ifelse(finalFallBP == "Inf" | is.na(finalFallBP), mean.fall.bp, finalFallBP))
  
  ########## Combine Datat ##########
  fullDataSync <- climateData %>%
    left_join(df_covariates_upstream, by=c('featureid')) %>%
    left_join(dplyr::select(tempDataSync, date, featureid, site, temp), by = c('date', 'featureid')) %>%
    left_join(featureid_huc8, by = c('featureid')) %>%
    left_join(featureid_lat_lon, by = c('featureid')) %>%
    dplyr::mutate(year = as.numeric(format(date, "%Y")),
                  airTemp = (tmax + tmin)/2) %>%
    left_join(dplyr::select(foo, -site), by = c('featureid', 'year')) %>%
    dplyr::mutate(finalSpringBP = ifelse(finalSpringBP == "Inf" | is.na(finalSpringBP), mean.spring.bp, finalSpringBP),
                  finalFallBP = ifelse(finalFallBP == "Inf" | is.na(finalFallBP), mean.fall.bp, finalFallBP)) %>%
    dplyr::mutate(dOY = yday(date)) %>%
    dplyr::filter(AreaSqKM < 200)
  # dplyr::filter(AreaSqKM >= 1 & AreaSqKM < 200 & allonnet < 70) # changed so don't deal with problematically small drainage areas (somre were 0.00006 sq km) - for loop didn't like this!!!!!!!!!
  
  ################### PROBLEM ################
  # if allonnet is very large (maybe > 75% of drainage area) the predictions are probably non-sense 
  ##########################
  
  ################### PROBLEM #################
  # 2-day precip as large as 210 - not sure if this is realistic and if so it might be outside the scope of our predictions
  ##################################
  
  
  rm(climateData)
  gc()
  
  # Order by group and date
  fullDataSync <- fullDataSync[order(fullDataSync$featureid, fullDataSync$year, fullDataSync$dOY),]
  
  # For checking the order of fullDataSync
  #fullDataSync$count <- 1:length(fullDataSync$year)
  
  #fullDataSync <- fullDataSync[order(fullDataSync$count),] # just to make sure fullDataSync is ordered for the slide function
  
  # moving means instead of lagged terms in the future
  fullDataSync <- fullDataSync %>%
    group_by(featureid, year) %>%
    arrange(featureid, year, dOY) %>%
    mutate(impoundArea = AreaSqKM * allonnet,
           airTempLagged1 = lag(airTemp, n = 1, fill = NA),
           temp5p = rollapply(data = airTempLagged1, 
                              width = 5, 
                              FUN = mean, 
                              align = "right", 
                              fill = NA, 
                              na.rm = T),
           temp7p = rollapply(data = airTempLagged1, 
                              width = 7, 
                              FUN = mean, 
                              align = "right", 
                              fill = NA, 
                              na.rm = T),
           prcp2 = rollsum(x = prcp, 2, align = "right", fill = NA),
           prcp7 = rollsum(x = prcp, 7, align = "right", fill = NA),
           prcp30 = rollsum(x = prcp, 30, align = "right", fill = NA))
  
  # clip to synchronized period of the year
  #dplyr::mutate(huc = huc8) %>%
  
  ############ PROBLEM ############################
  # not assigning breakpoints properly - mostly NA 
  #############
  
  fullDataSync <- fullDataSync %>%
    dplyr::filter(dOY >= finalSpringBP & dOY <= finalFallBP | is.na(finalSpringBP) | is.na(finalFallBP & finalSpringBP != "Inf" & finalFallBP != "Inf")) %>%
    dplyr::filter(dOY >= mean.spring.bp & dOY <= mean.fall.bp)
  
  var.names <- c("airTemp", 
                 "temp7p",
                 "prcp", 
                 "prcp2",
                 "prcp7",
                 "prcp30",
                 "dOY", 
                 "forest", 
                 "herbaceous", 
                 "agriculture", 
                 "devel_hi", 
                 "developed",
                 "AreaSqKM",  
                 "allonnet",
                 "alloffnet",
                 "surfcoarse", 
                 "srad", 
                 "dayl", 
                 "swe",
                 "impoundArea")
  
  fullDataSync <- fullDataSync %>%
    mutate(HUC8 = as.character(HUC8),
           huc8 = as.character(HUC8),
           huc = as.character(HUC12),
           site = as.numeric(as.factor(featureid))) 
  
  ################# PROBLEM ###############################
  # HUGE VALUES FOR STANDARDIZED COVARIATES 
  # consider adjusting those for agriculture and herbaceous
  # maybe also for precip
  ###########
  
  fullDataSyncS <- stdCovs(x = fullDataSync, y = df_stds, var.names = var.names)
  
  fullDataSyncS <- addInteractions(fullDataSyncS)
  
  fullDataSyncS <- dplyr::arrange(fullDataSyncS, featureid, date)
  
  fullDataSyncS$site <- as.character(fullDataSyncS$featureid)
  
  if(exists("fullDataSyncS$sitef")) {
    fullDataSyncS <- dplyr::select(fullDataSyncS, -sitef)
  }
  
  fullDataSyncS <- fullDataSyncS %>%
    left_join(rand_ids$df_site) %>%
    left_join(rand_ids$df_huc) %>%
    left_join(rand_ids$df_year)
  
  fullDataSyncS <- indexDeployments(fullDataSyncS, regional = TRUE)
  
  return(list(fullDataSync = fullDataSync, fullDataSyncS = fullDataSyncS))
}


#' @title prepPredictDF
#'
#' @description
#' \code{prepPredictDF} Helper function to prepare dataframe for predictions
#'
#' @param data Dataframe for which predictions will be calculated
#' @param cov.list List of covariates used in the model
#' @param coef.list List of coefficient values estimated from the model
#' @param featureid_site dataframe with featureid and site as 2 columns
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
prepPredictDF <- function(data, coef.list, cov.list, var.name, featureid_site) {
  B <- prepConditionalCoef(coef.list = coef.list, cov.list = cov.list, var.name = var.name)
  if(var.name == "site" | var.name == "ar1") {
    B$site <- as.factor(B$site)
    featureid_site$site <- as.factor(featureid_site$site)
    B <- dplyr::left_join(B, featureid_site, by = c("site"))
    if(var.name == "ar1") {
      var.name <- "site"
      data[ , var.name] <- as.character(data[ , var.name])
      B <- dplyr::select(B, site = site, featureid = featureid, B.ar1 = mean) # 
      #df <- left_join(data, B, by = c("featureid"))
      df <- merge(data, B, by = c("featureid"), all.x = T)
      df[ , names(B[-1])][is.na(df[ , names(B[-1])])] <- colMeans(dplyr::select(B, B.ar1)) # replace NA with mean
    } else {
      #B[ , var.name] <- as.character(B[ , var.name])
      #data[ , var.name] <- as.character(data[ , var.name])
      #df <- left_join(data, dplyr::select(B, -site), by = c("featureid")) # merge so can apply/mutate by rows without a slow for loop
      df <- merge(data, dplyr::select(B, -site), by = c("featureid"), all.x = T)
      for(i in 2:length(names(B))) {
        df[ , names(B[i])][is.na(df[ , names(B[i])])] <- colMeans(B[i])
      }
    }
  } else {
    B[ , var.name] <- as.character(B[ , var.name])
    data <- as.data.frame(unclass(data))
    data[ , var.name] <- as.character(data[ , var.name])
    df <- merge(data, B, by = var.name, all.x = T) # merge so can apply/mutate by rows without a slow for loop - merge with dplyr causes error if no overlap in by list
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
    B$site <- as.character(B$site)
  } else {
    f <- paste0(var.name, " ~ coef")
    B <- dcast(coef.list[[paste0("B.", var.name)]], formula = as.formula(f), value.var = "mean") # convert long to wide
    B <- dplyr::select(B, one_of(c(var.name, cov.list[[paste0(var.name, ".ef")]]))) # recorder to match for matrix multiplcation
    if(var.name == "site" | var.name == "ar1") B$site <- as.character(B$site)
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



