indexUpstreamDaymetVariablesForObservedSites <- function(record, variables, tiles, catchmentShapefile, covariateData, delineatedCatchmentsList, daymetDirectory){
  
  library(ncdf)
  library(rgeos)
  library(sp)
  
  start.time <- proc.time()[3]
  
  numTiles <- length(tiles)
  
  #------------------------------------------------------------------------------
  # Create a master list of the daymet coords by looping through all of the tiles
  #------------------------------------------------------------------------------
  for (i in 1:numTiles){
    
    # Open the NetCDF 
    NCDF <- open.ncdf(paste0(daymetDirectory, tiles[i], '_2010/prcp.nc'))    #netcdf   
    
    # Dimension limits of each of the variables
    start1 = c(1,1)
    latcount <- c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2])
    loncount <- c(NCDF$var$lon$varsize[1], NCDF$var$lon$varsize[2])
    
    # Read in variables
    lat = get.var.ncdf ( nc=NCDF, varid="lat", start = start1, count = latcount )
    lon = get.var.ncdf ( nc=NCDF, varid="lon", start = start1, count = loncount )
    close.ncdf(NCDF)
    
    # Join coordinate lists
    tempCoords <- cbind( as.vector(lon), as.vector(lat))
    colnames(tempCoords) <- c("Longitude", "Latitude")
    
    # Generate list
    if (i ==1) {masterCoords <- tempCoords}
    if (i > 1) {masterCoords <- rbind(masterCoords, tempCoords)}
  }
  
  masterCoordsMatrix <- masterCoords
  masterCoords <- as.data.frame(masterCoords)
  
  #----------------------------------------------
  # Loop through Sites and NetCDFs, getting data.
  #----------------------------------------------  
  
  sites <- unique(masterData$site)
  
  masterLength <- length(delineatedCatchmentsList)
  
  for ( i in 1:length(sites)){
    
    print(paste0(round(i/length(sites), digits = 3)*100, '% done.'))
    
    # Define the catchment polygon:
    #------------------------------
    featureID <- covariateData$FEATUREID[covariateData$site %in% sites[i]]
    features <- delineatedCatchmentsList[[which(sapply(c(1:masterLength),FUN = function(x){delineatedCatchmentsList[[x]][1]==featureID})==TRUE)]]
    
    catchmentShape <- catchmentShapefile[catchmentShapefile$FEATUREID %in% features,]
    basinShape     <- gUnaryUnion(catchmentShape) #dissolve individual catchments
    a <- SpatialPoints(masterCoords, proj4string=CRS(proj4.NHD))
    
    inside <- as.data.frame(a[!is.na(over(a, basinShape)),])
    
    
    #If no point falls within the catchment, find the nearest one:
    #-------------------------------------------------------------
    if(nrow(inside) == 0 ){
      
      tempLat <- covariateData$Latitude [covariateData$site %in% sites[i]]
      tempLon <- covariateData$Longitude[covariateData$site %in% sites[i]]
      
      distances <- spDistsN1(masterCoordsMatrix, c(tempLon, tempLat), longlat = TRUE)
      minDist <- min(distances)
      distpos <- which(distances == minDist)[1]
      
      nearLon  <- masterCoords[distpos, 1]
      nearLat  <- masterCoords[distpos, 2]
      
      inside[1,1] <- nearLon
      inside[1,2] <- nearLat
    }
    
    
    #Pull the Tiles for the points within the catchment:
    #---------------------------------------------------
    for(k in 1:length(inside[,1])){
      
      siteLon <- inside[k,1]
      siteLat <- inside[k,2]
      
      #Index the tile by site location:
      tile <- indexDaymetTileByLatLon(siteLat,siteLon)
      
      temp <- data.frame('Longitude' = siteLon, 'Latitude' = siteLat, 'Tile' = tile)
      
      if (k ==1) spatialLocs <- temp
      if (k > 1) spatialLocs <- rbind(spatialLocs, temp)
    }
    rm(siteLat, siteLon)
    
    spatialLocs <- spatialLocs[ order(spatialLocs$Tile), ]
    subTiles <- unique(spatialLocs$Tile)
    
    
    #Link record with catchment area:
    #--------------------------------
    curRecord <- masterData[which(masterData$site == sites[i]),]
    curRecord <- curRecord[which(curRecord$year < 2014),]
    
    begYr <- min(curRecord$year)
    
    for (j in 1:length(variables)){
      
      for ( year in unique(curRecord$year) ){
        
        for ( t in 1:length(subTiles)){
          
          NCDF <- open.ncdf(paste0(daymetDir, subTiles[t], '_', year,'/', variables[j], '.nc'))    #netcdf
          
          # Dimension limits of each of the variables we'll use:
          # ----------------------------------------------------
          start1 = c(1,1)
          latcount <- c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2])
          loncount <- c(NCDF$var$lon$varsize[1], NCDF$var$lon$varsize[2])
          YDcount  <- NCDF$var$yearday$varsize      
          
          start2 = c(1, 1, 1)
          varcount = c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2], NCDF$var$yearday$varsize)
          
          # Read in variables:
          # ------------------
          lat = get.var.ncdf ( nc=NCDF, varid="lat",                 start = start1, count = latcount )
          lon = get.var.ncdf ( nc=NCDF, varid="lon",                 start = start1, count = loncount )
          dOY = get.var.ncdf ( nc=NCDF, varid="yearday",             start = 1,      count = YDcount  )
          var = get.var.ncdf ( nc=NCDF, varid= paste0(variables[j]), start = start2, count = varcount )
          
          close.ncdf(NCDF)
          
          # Correction for Daymet doy which starts at 0.
          dOY <- dOY + 1  
          
          tileCoords <- as.data.frame(cbind( as.vector(lon), as.vector(lat)))
          names(tileCoords) <- c('Lon', 'Lat')
          
          xx <- spatialLocs[which(spatialLocs$Tile == subTiles[t]),]
          
          for (m in 1:length(xx[,1])){
            position <- which(lon == xx$Longitude[m] & lat == xx$Latitude[m], arr.in = TRUE)
            
            if ( t == 1 & m == 1) {tempVar <- data.frame(year, dOY, var[position[1], position[2], 1:365])} else(tempVar <- cbind(tempVar, var[position[1], position[2], 1:365]))
            
          }
        }
        
        
        # Switch to na.rm = T
        ifelse( ncol(tempVar) > 3, R <- rowMeans(tempVar[,-c(1,2)], na.rm = FALSE, dims = 1),  R <- tempVar[,-c(1,2)] )
        
        tempVar <- data.frame(tempVar[,c(1,2)], R)
        names(tempVar) <- c("year", "dOY", paste0(variables[j]))
        
        if ( year == begYr ) ( mainVar <- tempVar)
        if ( year >  begYr ) ( mainVar <- rbind(mainVar, tempVar))
        
        rm(tempVar, R)
      }
      
      if (j == 1) {allVars <- mainVar} else {allVars <- merge(allVars, mainVar, by = c('year','dOY'), all.x = T)}  
    }
    
    #Add data into main dataframe
    allVars$site <- sites[i]
    tempRecord <- merge(curRecord, allVars, by = c("site", "year", "dOY"), all.x = T, all.y = F, sort = F)
    if (i == 1) {fullRecord <- tempRecord} else {fullRecord <- rbind(fullRecord, tempRecord)}
  }
  
  # How long it takes to run
  end.time   <- proc.time()[3]
  print(paste0((end.time-start.time)/3600, " hours"))
  
  return(fullRecord)
  
}
#=========================================================================================================
