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

# this works for modes with no interactions, but not for those with them
# deosn't seems all that useful anyway because a linear model is fine for thses data
predSim <- function( nSim=100, model, data ){
  
  # predictive simulation for model checking
  # simulated parameter estimates
  
  s <- sim( model,nSim ) # gets 100 parameter estimates and sigmas from the model
  sC <- coef(s)
  sSigma <- s@sigma #sigma method doesn't work
  
  p <- array(NA,c(nrow(data),nSim))
  pWError <- p 
  etArray <- data.matrix(cbind(rep(1,nrow(data)),data[,fixCov])) # this call to data works well for modes with no interactions, not sure how to get the interactions in the data
  for(i in 1:nrow(data)){
    
    p[i,] <- sC %*% etArray[i,]
    pWError[i,]  <- rnorm(nSim,p[i,],sSigma)
    
  }
  pRowMeans <- apply(p,1,mean)
  pRowSDs <- apply(p,1,sd)
  
  pWErrorRowMeans <- apply(pWError,1,mean)
  pWErrorRowSDs <- apply(pWError,1,sd)
  
  gP <- 
    ggplot(cbind(pRowMeans,data), aes(pRowMeans,temp)) +
    geom_point()+
    geom_abline(intercept=0, slope=1,color='white')
    #ggtitle('data')
  
  gPWError <- 
    ggplot(cbind(pWErrorRowMeans,data), aes(pWErrorRowMeans,temp)) +
    geom_point()+
    geom_abline(intercept=0, slope=1,color='white')
  
  return( list(pRowMeans=pRowMeans,pRowSDs=pRowSDs,
               pWErrorRowMeans=pWErrorRowMeans,pWErrorRowSDs=pWErrorRowSDs,
               gP=gP,gPWError=gPWError) )
  
}


# leave out one site validation
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

# requires data produced from validateModel()

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






#======================================
# Scale time series data (from Daymet):
#======================================
# This function takes a dataframe of time series data and lists of variables to either scale or log-scale.


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



#========================================
# Predict the maximum summer temperature:
#========================================
# This function takes the two dataframes with prediction data for the segmented models and a dataframe of all the prediction outputs.
# It runs the respective models on each segment dataframe and joins them together for one output.
# It also computes the confidence interval for this prediction.
# Predicted values are added to the dataframe of prediction output.

predictSummerTMax <- function(seg2DF, seg3DF, outputDF){
  
  # Index the summer BPs:
  x2 <- ddply( seg2DF, .(site), summarise, dOY = dOY[dOY == round(bp2)])
  x3 <- ddply( seg3DF, .(site), summarise, dOY = dOY[dOY == round(bp2)])
  
  # Pull the prediction data for these values:
  y2 <- merge(x2, seg2DF, by = c('site', 'dOY'), all.x = T, sort = F)
  y3 <- merge(x3, seg3DF, by = c('site', 'dOY'), all.x = T, sort = F)
  
  # Predict the temp and CI for just the summer BP:
  z2 <-  data.frame(y2$site, predict(finalModelS2, y2, interval = "confidence"))
  z3 <-  data.frame(y3$site, predict(finalModelS2, y3, interval = "confidence"))
  
  # Names :
  names(z2) <- c('site', 'summerMax', 'lwr', 'upr')
  names(z3) <- c('site', 'summerMax', 'lwr', 'upr')
  
  # Join the segments:
  comb <- rbind(z2, z3)
  
  # Calculate the Confidence Interval:
  comb$summerMaxCI <- comb$upr - comb$lwr
  
  # Merge with other predictions:
  out <- merge(outputDF, comb[,c('site', 'summerMax', 'summerMaxCI')], by = 'site', all.x = T)
  
  return(out)
}



#===================================
# Predict rising and falling slopes:
#===================================
# This function takes the one dataframe that already contains predicted temperatures and one dataframe of prediction outputs.
# It runs a linear regression (by site and segment) on these predictions vs. the observed air temperature.
# The slopes are returned and added to the prediction dataframe output.

# Note: We can't do this directly in the big regressions because other daily covariates besides airTemp are included (swe, dayl, srad)

predictSlopes <- function(inputDF, outputDF){
  
  # Predict rising and falling slopes
  # ---------------------------------
  slopes <- ddply( inputDF, .(site,segment), summarize, slope=coef(lm(predTemp ~ airTemp))[2])
  
  # Manipulate the output to join to prediction dataframe:
  slopesMelt <- melt(slopes,id.vars=c('site','segment'))
  slopesCast <- cast(slopesMelt,site~segment)
  names(slopesCast) <- c('site','slopeSeg2','slopeSeg3')
  
  # Merge:    
  out1 <- merge( x=outputDF, y=slopesCast, all.x=T, by = c('site') )
  
  # Get the confidence interval for the slopes
  # ------------------------------------------
  slopesCI <- ddply( inputDF, .(site,segment), summarize, range=confint(lm(predTemp ~ airTemp), 'airTemp', level = 0.95))      
  slopesCI$slopeCI <- slopesCI$range[,2] - slopesCI$range[,1]
  
  # Names:
  slopesCI <- slopesCI[,c('site', 'segment', 'slopeCI')]
  
  # Manipulate output:
  slopesMeltCI <- melt(slopesCI,id.vars=c('site','segment'))
  slopesCastCI <- cast(slopesMeltCI,site~segment)
  names(slopesCastCI) <- c('site','slopeSeg2CI','slopeSeg3CI')
  
  # Merge:
  out <- merge( x=out1, y=slopesCastCI, all.x=T, by = c('site') )
  
  return(out)
}

# Function that returns Root Mean Squared Error
rmse <- function(error) {
  sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error) {
  mean(abs(error))
}

