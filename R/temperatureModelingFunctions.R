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



#' @title rmse: root mean squared error
#'
#' @description
#' \code{rmse} returns root mean squared error
#'
#'@param error Vector of residual error from a model
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
rmse <- function(error) {
  sqrt(mean(error^2))
}

#' @title mae: root mean squared error
#'
#' @description
#' \code{mae} returns mean absolute error
#'
#'@param error Vector of residual error from a model
#' @details
#' var: blah, blah, blah
#' value: something, something
#' @export
mae <- function(error) {
  mean(abs(error))
}

