#' @title predictTemp: Predict conditional daily temperatures for observed or unobserved sites
#'
#' @description
#' \code{predictTemp} Predict daily stream temperatures
#'
#' @param data Data frame of covariates for prediction
#' @param data.fit Data frame used for fitting (calibrating) the JAGS model
#' @param firstObsRows Dataframe with the rowNum column indicating the rows where a logger was first deployed
#' @param evalRows Dataframe with the rowNum column indicating the rows during a deployment after the first day for use in the autoregressive
#' @param cov.list List of covariates used in the model
#' @param coef.list List of coefficient values estimated from the model
#' @return Numeric vector of predicted daily stream temperatures
#' @details
#' The predictions are all conditional on site, HUC8, and year when the data is available and predicts the mean values when any component was not used in the fitting (not in the calibration dataset)
#' @examples
#' 
#' \dontrun{
#' Predictions <- predictTemp(data = tempDataSyncValidS, data.fit = tempDataSyncS, firstObsRows = firstObsRows, evalRows = evalRows, cov.list = cov.list, coef.list = coef.list))
#' }
#' @export
predictTemp <- function(data, data.fit = tempDataSyncS, coef.list, cov.list, firstObsRows, evalRows, observed = TRUE) {
  
  B.fixed = coef.list$B.fixed
  B.site = coef.list$B.site
  B.huc = coef.list$B.huc
  B.year = coef.list$B.year
  B.ar1 = coef.list$B.ar1
  
  df <- prepDF(data, covars = cov.list)
  
  if(!(identical(data, data.fit))) {
    site.new <- data %>%
      dplyr::filter(!(site %in% levels(as.factor(data.fit$site)))) %>%
      dplyr::select(site)
    B.site.new <- expand.grid(site = levels(as.factor(site.new$site)), coef = levels(B.site$coef))
    B.site.new$mean <- 0
    B.site <- rbind(B.site[ , c("site", "coef", "mean")], B.site.new[ , c("site", "coef", "mean")])
    
    huc.new <- data %>%
      dplyr::filter(!(HUC8 %in% levels(as.factor(data.fit$HUC8)))) %>%
      dplyr::select(huc = HUC8) 
    B.huc.new <- expand.grid(huc = levels(as.factor(huc.new$huc)), coef = levels(B.huc$coef))
    mu.huc <- dplyr::filter(coef.list$fix.ef, grepl('^mu.huc', Parameter))
    mu.huc$coef <- colnames(df$data.random.sites)
    B.huc.new <- dplyr::left_join(mu.huc, B.huc.new, by = "coef")
    B.huc <- rbind(B.huc[ , c("huc", "coef", "mean")], B.huc.new[ , c("huc", "coef", "mean")])
    
    year.new <- data %>%
      dplyr::filter(!(year %in% levels(as.factor(data.fit$year)))) %>%
      dplyr::select(year = year) 
    B.year.new <- expand.grid(year = levels(as.factor(year.new$year)), coef = levels(B.year$coef))
    mu.year <- dplyr::filter(coef.list$fix.ef, grepl('^mu.year', Parameter))
    mu.year$coef <- colnames(df$data.random.years)
    B.year.new <- dplyr::left_join(mu.year, B.year.new, by = "coef")
    B.year <- rbind(B.year[ , c("year", "coef", "mean")], B.year.new[ , c("year", "coef", "mean")])
    
    ar1.new <- data %>%
      dplyr::filter(!(site %in% levels(as.factor(data.fit$site)))) %>%
      dplyr::select(site)
    B.ar1.new <- data.frame(site = levels(as.factor(ar1.new$site)))
    B.ar1.new$mean <- mean(B.ar1$mean)
    B.ar1 <- rbind(B.ar1[ , c("site", "mean")], B.ar1.new)
  }
  
  Pred <- NA
  
  if(observed) {
    trend <- NA
    for(i in 1:length(firstObsRows)) {
      trend[firstObsRows[i]] <- as.numeric(
        B.fixed$mean %*% t(df$data.fixed[firstObsRows[i], ]) + 
          dplyr::filter(B.huc, huc == as.character(data$HUC8[firstObsRows[i]]))$mean %*% t(df$data.random.sites[firstObsRows[i], ]) + 
          dplyr::filter(B.site, site == as.character(data$site[firstObsRows[i]]))$mean %*% t(df$data.random.sites[firstObsRows[i], ]) + 
          dplyr::filter(B.year, year == as.character(data$year[firstObsRows[i]]))$mean %*% t(df$data.random.years[firstObsRows[i], ])
      )
      
      Pred[firstObsRows[i]] <- trend[firstObsRows[i]]
    }
    
    for(i in 1:length(evalRows)) {
      trend[evalRows[i]] <- as.numeric(
        B.fixed$mean %*% t(df$data.fixed[evalRows[i], ]) + 
          dplyr::filter(B.huc, huc == as.character(data$HUC8[evalRows[i]]))$mean %*% t(df$data.random.sites[evalRows[i], ]) + 
          dplyr::filter(B.site, site == as.character(data$site[evalRows[i]]))$mean %*% t(df$data.random.sites[evalRows[i], ]) + 
          dplyr::filter(B.year, year == as.character(data$year[evalRows[i]]))$mean %*% t(df$data.random.years[evalRows[i], ]) 
      )
      
      if(is.na(data$temp[evalRows[i]-1]) | is.null(data$temp[evalRows[i]-1])) {
        Pred[evalRows[i]] <- trend[evalRows[i]]
      } else {
        Pred[evalRows[i]] <- trend[evalRows[i]] + 
          dplyr::filter(B.ar1, site == as.character(data$site[evalRows[i]]))$mean * (data$temp[evalRows[i]-1] - trend[evalRows[i]-1])
      }
      #as.numeric(B.ar1$mean * (data$temp[evalRows[i]-1] - trend[evalRows[i]-1]))
    }
    
  } else {
    for(i in 1:dim(data)[1]) {
      Pred[i] <- as.numeric(
        B.fixed$mean %*% t(df$data.fixed[i, ]) + 
          dplyr::filter(B.huc, huc == as.character(data$HUC8[i]))$mean %*% t(df$data.random.sites[i, ]) + 
          dplyr::filter(B.site, site == as.character(data$site[i]))$mean %*% t(df$data.random.sites[i, ]) + 
          dplyr::filter(B.year, year == as.character(data$year[i]))$mean %*% t(df$data.random.years[i, ])
      )
      
    }
  }
  
  return(Pred)
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
plotPredict <- function(observed, predicted, siteList = "ALL", yearList = "ALL", dir = getwd()){ # add option to not include predicted or make similar function that makes observation plots
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
  predicted.origin.scale <- left_join(observed, predicted[ , c("site", "date", "tempPredicted")], by = c("site", "date"))
  
  for(i in 1:length(unique(sites))){
    dataSite <- dplyr::filter(predicted.origin.scale, filter = site == sites[i] & year %in% years)
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
      theme(axis.text.x = element_text(angle = 45))
    ggsave(filename=paste0(dir, dataSite$site[i], '.png'), plot=foo, dpi=300 , width=12,height=8, units='in' )
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


