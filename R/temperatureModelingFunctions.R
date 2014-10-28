#' @title riverLabeller: Label Rivers
#'
#' @description
#' \code{riverLabeller} returns something of great importance
#'
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
riverLabeller <- function(var, value){
  value <- as.character(value)
  if (var=="site") { 
    value[value=="1"] <- "WB"
    value[value=="2"] <- "OL"
    value[value=="3"] <- "OS"
    value[value=="4"] <- "Is"
  }
  return(value)
}

riverLabeller2 <- function(var, value){
  value <- as.character(value)
  if (var=="site") { 
    value[value=="WEST BROOK"] <- "WB"
    value[value=="WB JIMMY"] <- "OL"
    value[value=="WB MITCHELL"] <- "OS"
    value[value=="WB OBEAR"] <- "Is"
  }
  return(value)
}



#' @title validateModel: Leave 1 site out model validation
#'
#' @description
#' \code{validateModel} returns something of great importance
#'
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
validateModel <- function ( model, d, vd ) {
  f=model$call
  siteList <- unique(vd$site)
  siteYearList <- unique( vd[,c('site','year')])
  
  v <- data.frame(s=NA,y=NA,n=NA,rmse=NA,biasMean=NA,biasSD=NA,
                  sigma=NA,lat=NA,lon=NA)
  g <- list()
  i=0
  for(s in siteList){
    
    # run model without site = s
    withOut <- d[d$site != s,]    
    m <- lm( f, data=withOut )
    
    yearList <- unique( vd$year[vd$site == s ] )
    
    for(y in yearList){
     
      # data for the site/year combo
      with <- vd[vd$site == s & vd$year == y,]
      # predict temp for missing site/year based on all the other sites
      p <- cbind(with,pred=predict(m,with))
      
      i=i+1
      g[[i]] <- ggplot(p[,c('temp','pred')],aes(temp,pred))+
        geom_point()+
        geom_smooth(method='lm') +
        geom_abline(intercept=0,slope=1) +
        scale_x_continuous('Observed temperature') +
        scale_y_continuous('Predicted temperature') +
        ggtitle(paste(as.character(s), as.character(y),sep='_'))
      
      p$bias <- p$pred-p$temp
# #      v <- rbind(v,c(as.character(s), as.character(y), nrow(p), 
#                      sqrt(sum(p$bias^2, na.rm=T)/nrow(p)),
#                      mean(p$bias,na.rm=T),
#                      sd(p$bias,na.rm=T),
#                      summary(m)$sigma,
#                      unique(p$Latitude),unique(p$Longitude))) 
      v <- rbind(v,data.frame(s=as.character(s), 
                              y=y, n=nrow(p), 
                              rmse=sqrt(sum(p$bias^2, na.rm=T)/nrow(p)),
                              biasMean=mean(p$bias,na.rm=T),
                              biasSD=sd(p$bias,na.rm=T),
                              sigma=summary(m)$sigma,
                              lat=unique(p$Latitude),lon=unique(p$Longitude)))
      
     print(paste0(s," ", i, ' out of ', nrow(siteYearList),"   rmse: ", v$rmse[i+1]  ))
    }
  }
  v <- v[-1,] #get rid of first dummy row
  v$biasMeanAbs <- abs(v$biasMean)
  v$biasMeanDir <- ifelse( v$biasMean>0,1,-1 )

  meanRMSE <- mean((v$rmse),na.rm=T)
  meanBiasMean <- mean((v$biasMean),na.rm=T)
  meanBiasSD <- mean((v$biasSD),na.rm=T)

  return( list(v=v,g=g,means=list(meanRMSE=meanRMSE,meanBiasMean=meanBiasMean,meanBiasSD=meanBiasSD) ) )
}


#' @title makeBiaseMap
#'
#' @description
#' \code{makeBiaseMap} returns something of great importance
#'
#' @details
#' requires data produced from validateModel()
#' var: blah, blah, blah
#' value: something, something
#' @export
makeBiasMap <- function (d) {
  library(ggmap)
  #map.center <- geocode("Hartford, CT")
  map.center <- geocode("Springfield, MA")
  baseMap <- qmap(c(lon=map.center$lon, lat=map.center$lat), source="google", zoom=8)
  map <- baseMap + geom_point(aes(x=lon, y=lat, 
                                  size=(biasMeanAbs),
                                  color=factor(biasMeanDir)), 
                              data=d$v) +
    scale_color_manual(values=c('darkred','darkgreen')) 
    #ggtitle( as.character(y)))
  return(map)
}



#' @title indexUpstreamDaymetVariablesForObservedSites
#'
#' @description
#' \code{indexUpstreamDaymetVariablesForObservedSites} This function takes a dataframe of time series data and lists of variables to either scale or log-scale.
#'
#' @param dataFrame dataframe
#' @param logVars somthing
#' @param scaleVars something
#' 
#' @return Returns something
#' @details
#' This function takes a dataframe of time series data and lists of variables to either scale or log-scale.
#' @export
loggAndScaleDaymet <- function(dataFrame, logVars, scaleVars){
  
  # Don't use "scale" because it creates a data type with attributes that make is difficult to use predict
  
  # Log the specified time-series data:
  for ( i in logVars ){
    dataFrame[,paste0(i, 'L')] <- log(dataFrame[,i]  + 0.001)
  }
  
  #Split the dataframe into segments:
  dataFrameSeg2 <- dataFrame[dataFrame$segment %in% 2,]
  dataFrameSeg3 <- dataFrame[dataFrame$segment %in% 3,]
  
  # Scale the specified time-series data:
  for ( j in scaleVars ){
    
    if( j %in% logVars ) {dataFrameSeg2[,paste0(j, 'LS')]  <- (dataFrameSeg2[,paste0(j, 'L')]  - mean(dataFrameSeg2[,paste0(j, 'L')],na.rm=T) ) / sd(dataFrameSeg2[,paste0(j, 'L')], na.rm=T )} 
    else (dataFrameSeg2[,paste0(j, 'S')]  <- (dataFrameSeg2[,j]  - mean(dataFrameSeg2[,j],na.rm=T) ) / sd(dataFrameSeg2[,j], na.rm=T ))
    
    if( j %in% logVars ) {dataFrameSeg3[,paste0(j, 'LS')]  <- (dataFrameSeg3[,paste0(j, 'L')]  - mean(dataFrameSeg3[,paste0(j, 'L')],na.rm=T) ) / sd(dataFrameSeg3[,paste0(j, 'L')], na.rm=T )} 
    else (dataFrameSeg3[,paste0(j, 'S')]  <- (dataFrameSeg3[,j]  - mean(dataFrameSeg3[,j],na.rm=T) ) / sd(dataFrameSeg3[,j], na.rm=T ))  
  }
  
  return(list(dataFrame, dataFrameSeg2, dataFrameSeg3))
}



#' @title prepDF: Prepares dataframe for use or prediction
#'
#' @description
#' \code{prepDF} Prepares dataframe for use or prediction
#'
#' @param data Dataframe of potential covariates to model or predict daily stream temperature
#' @param form Named list of formulae each created with the formula function. Names must be data.fixed, data.random.sites, and data.random.years.
#' @return List of named data frames each created using model.matrix and the input formulae
#' @details
#' blah, blah, blah
#' @examples
#' 
#' \dontrun{
#' data.list <- prepDF(data = df, form = formulae)
#' }
#' @export
prepDF <- function(data, form = formulae) {
  ## Change to use model matrix for simplicity and consistency?
  #data.fixed <- data.frame(intercept = 1,
  #                   lat = data$Latitude,
  #                  lon = data$Longitude,
  #                 drainage = data$TotDASqKM,
  #                forest = data$Forest,
  #               elevation = data$ReachElevationM,
  #              coarseness = data$SurficialCoarseC,
  #             wetland = data$CONUSWetland,
  #            impoundments = data$ImpoundmentsAllSqKM,
  #           airT.drainage = data$airTemp * data$TotDASqKM)
  
  #data.fixed <- data %>%
  # mutate(intercept = 1, 
  #       airT.drainage = airTemp * TotDASqKM) %>%
  #select(intercept,
  #      Latitude,
  #     Longitude,
  #    TotDASqKM,
  #   Forest,
  #  ReachElevationM,
  # SurficialCoarseC,
  # CONUSWetland,
  # ImpoundmentsAllSqKM,
  # airT.drainage)
  
  # Benefit of using model matrix?
  
  data.fixed <- model.matrix(form$fixed.form, data)
  # str(data.fixed)
  # head(data.fixed)
  colnames(data.fixed)
  
  
  data.random.sites <- model.matrix(form$site.form, data)
  
  #data.random.sites <- data.frame(intercept.site = 1, 
  #                    airTemp = data$airTemp, 
  #airTempLag1 = data$airTempLagged1,
  #                   airTempLag2 = data$airTempLagged2,
  #                  precip = data$prcp,
  #                 precipLag1 = data$prcpLagged1,
  #precipLag2 = data$prcpLagged2,
  # effect of interaction with random and fixed effects?
  #                airT.precip = data$airTemp * data$prcp) # airT.precip2 = data$airTemp * data$prcpLagged1 - doesn't converge
  
  
  data.random.years <- model.matrix(form$year.form, data)
  
  #data.random.years <- data.frame(intercept.year = 1, 
  #                    dOY = data$dOY, 
  #                    dOY2 = data$dOY^2,
  #                   dOY3 = data$dOY^3)
  
  return(list(data.fixed = as.data.frame(data.fixed), data.random.sites = as.data.frame(data.random.sites), data.random.years = as.data.frame(data.random.years)))
}
