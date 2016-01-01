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
predictTemp <- function (data, coef.list, cov.list, rand_ids, Random_AR1 = FALSE) {
  if (!exists("data$sitef")) {
    data <- data %>% left_join(rand_ids$df_site)
  }
  if (!exists("data$hucf")) {
    data <- data %>% left_join(rand_ids$df_huc)
  }
  if (!exists("data$yearf")) {
    data <- data %>% left_join(rand_ids$df_year)
  }
  fixed.ef <- as.numeric(coef.list$B.fixed$mean)
  data <- left_join(data, coef.list$B.site)
  data <- left_join(data, coef.list$B.huc)
  data <- left_join(data, coef.list$B.year)
  for (j in 2:length(names(coef.list$B.site))) {
    data[, names(coef.list$B.site[j])][is.na(data[,names(coef.list$B.site[j])])] <-
      colMeans(coef.list$B.site[j])
  }
  for (j in 2:length(names(coef.list$B.huc))) {
    data[, names(coef.list$B.huc[j])][is.na(data[,
                                                 names(coef.list$B.huc[j])])] <-
      colMeans(coef.list$B.huc[j])
  }
  for (j in 2:length(names(coef.list$B.year))) {
    data[, names(coef.list$B.year[j])][is.na(data[,
                                                  names(coef.list$B.year[j])])] <-
      colMeans(coef.list$B.year[j])
  }
  data$fixed.ef <-
    as.vector(fixed.ef %*% t(as.matrix(as.data.frame(unclass(
      select(ungroup(data),
             one_of(cov.list$fixed.ef))
    )))))
  data$site.ef <- rowSums(as.matrix(select(
    data,
    one_of(cov.list$site.ef)
  )) * as.matrix(select(data,
                        starts_with("B.site"))))
  data$huc.ef <- rowSums(as.matrix(select(
    data,
    one_of(cov.list$huc.ef)
  )) * as.matrix(select(data,
                        starts_with("B.huc"))))
  data$year.ef <- rowSums(as.matrix(select(
    data,
    one_of(cov.list$year.ef)
  )) * as.matrix(select(data,
                        starts_with("B.year"))))
  data$trend <-
    data$fixed.ef + data$site.ef +
    data$huc.ef + data$year.ef
  data <- mutate(
    data, prev.temp = c(NA,
                        data$temp[(2:(nrow(data))) - 1]), prev.trend = c(NA,
                                                                         data$trend[(2:nrow(data)) - 1]), prev.err = prev.temp -
      prev.trend, tempPredicted = trend, prev.temp = ifelse(newDeploy ==
                                                              1, NA, prev.temp), prev.err = ifelse(newDeploy == 1,
                                                                                                   NA, prev.err)
  )
  if (Random_AR1) {
    B.ar1.sub <- coef.list$B.ar1 %>% dplyr::select(sitef,
                                                   mean) %>% dplyr::rename(B.ar1 = mean)
    data <- left_join(data, B.ar1.sub,
                      by = c("sitef"))
    data <-
      data %>% dplyr::mutate(B.ar1 = ifelse(is.na(B.ar1),
                                            mean(B.ar1.sub$B.ar1, na.rm = T), B.ar1)) %>% dplyr::arrange(featureid,
                                                                                                         date)
  }
  else {
    data$B.ar1 <- coef.list$B.ar1$mean
    data <- data %>% dplyr::arrange(featureid,
                                    date)
  }
  data[which(!is.na(data$prev.err)),]$tempPredicted <-
    data[which(!is.na(data$prev.err)),]$trend + data[which(!is.na(data$prev.err)),]$B.ar1 * data[which(!is.na(data$prev.err)),]$prev.err
  data <- dplyr::arrange(data, rowNum)
  data <-
    data.frame(unclass(data), stringsAsFactors = FALSE)
  return(data)
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



