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
  
  for(i in 1:length(unique(sites))){
    dataSite <- filter(predicted, filter = site == sites[i], year == years)
    dataSiteObs <- filter(observed, filter = site == sites[i], year == years)
    foo <- ggplot(dataSite, aes(dOY, tempPredicted)) + 
      coord_cartesian(xlim = c(100, 300), ylim = c(0, 35)) + 
      geom_point(data=dataSiteObs, aes(dOY, temp), colour='black') +
      geom_point(colour = 'blue') + 
      geom_line(colour = 'blue') + 
      geom_point(aes(dOY, airTemp), colour = 'red') + 
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


