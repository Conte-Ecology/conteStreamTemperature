# Functions used in working with the Daymet climate data


#=========================================================================================================
# This function indexes values from the master list of covariates for observed stream temperature sites.
# It takes the following:
#    1) The stream temperature record (unique site ID, latitude, longitude, columns)
#    2) A dataframe of the covariates for the catchments (FEATUREIDs source)
#    3) A master catchments shapefile
#    4) A CRS string of the spatial data projection
#    5) A string of variables to pull from the covariates list
#
# It returns a dataframe with the site name, lat/lon, FEATUREID, and the select covariate values.
#=========================================================================================================
indexCovariateData <- function(record, masterCovariates, catchmentShapefile, projectionString, fields){
  start.time <- proc.time()[3]
  
  library(sp)
  
  # Select fields
  selectUpstreamStats <- UpstreamStats[,names(masterCovariates) %in% fields]
  
  # Generate a list of sites
  sites <- unique(masterData[,c('site', 'Latitude', 'Longitude')])
  
  # Loop through all sites in the record
  for ( i in 1:nrow(sites)){  
    
    # Select a site
    curSite <- sites[i,]
    
    # Make the site a SpatialPoints object
    point <- SpatialPoints(matrix(data=c(curSite$Longitude,curSite$Latitude),ncol=2,byrow=T),  proj4string=CRS(proj4.NHD))
    
    # Select the catchment that contains the point
    featureID <- over(point,catchmentShapefile)$FEATUREID
    
    # Index the covariate data and join it to the site info
    tempCovs <- data.frame( curSite, selectUpstreamStats[selectUpstreamStats$FEATUREID == featureID, names(selectUpstreamStats) %in% fields])
    
    # Store data from each iteration
    if ( i == 1 ) { newCovs <- tempCovs} else( newCovs <- rbind(newCovs, tempCovs) )
    
    # Track progress
    print(paste0(round(i/nrow(sites), digits = 3)*100, '% done.'))
  }
  
  # Return the covariates list
  return(outCovs)
  
  # How long it takes to run
  end.time   <- proc.time()[3]
  print(paste0((end.time-start.time)/3600, " hours"))
}
#=========================================================================================================



#=========================================================================================================
# This function indexes Daymet tiles by the latitude and longitude of a site.
# It takes a latitude and longitude value.
# It returns the Daymet tile that holds data for the coordinates given.
#=========================================================================================================
indexDaymetTileByLatLon <- function(SiteLat, SiteLon){

  Tile <- ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -74 & SiteLon < -72, 11754, #**
          ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -72 & SiteLon < -70, 11755, #**
          ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -70 & SiteLon < -68, 11756, #**   
          ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -74 & SiteLon < -72, 11934, #**
          ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -72 & SiteLon < -70, 11935, #**
          ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -70 & SiteLon < -68, 11936,
          ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -74 & SiteLon < -72, 12114, #**      
          ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -72 & SiteLon < -70, 12115, #** 
          ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -70 & SiteLon < -68, 12116,
          ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -68 & SiteLon < -66, 12117,
          ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -72 & SiteLon < -70, 12295,     
          ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -70 & SiteLon < -68, 12296,     
          ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -68 & SiteLon < -66, 12297,   
          print("Tile Error"))))))))))))))

  return(Tile)
}
# ** Explicitly checked with mapping software (ArcGIS).
#=========================================================================================================


#=========================================================================================================
# This function pulls the Daymet variables for a time series at one location nearest to the site location.
# It takes the stream temperature record and the string of variables to pull from Daymet.
#     At minimum, the record needs a unique site ID, latitude, longitude, year, and dOY columns.
# It returns the original dataframe with new columns for the Daymet variables.
#=========================================================================================================
indexLocalDaymetVariablesForObservedSites <- function(record, variables){
  start.time <- proc.time()[3]
  
  sites <- unique(record$site)
  
  for ( i in 1:length(sites)){
    
    print(paste0(round(i/length(sites), digits = 3)*100, '% done.'))
    
    # Select site
    curRecord <- record[which(record$site %in% sites[i]),]
    curRecord <- curRecord[which(curRecord$year < 2014),]
    
    # Site coordinates
    siteLon <- unique(curRecord$Longitude)
    siteLat <- unique(curRecord$Latitude)
    
    #Index the tile by site location:
    tile <- indexDaymetTileByLatLon(siteLat,siteLon)
    
    # Site time range
    begYr <- min(curRecord$year)
    
    # Loop through the variables and years in NetCDF files
    for (j in 1:length(variables)){
      
      for ( year in unique(curRecord$year) ){
        
        #Open the NetCDF with the known location:
        #--------------------------------------------
        NCDF <- open.ncdf(paste0(daymetDir, tile, '_', year,'/', variables[j], '.nc'))    #netcdf
        
        #Dimension limits of each of the variables we'll use:
        #----------------------------------------------------
        start1 = c(1,1)
        latcount <- c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2])
        loncount <- c(NCDF$var$lon$varsize[1], NCDF$var$lon$varsize[2])
        YDcount  <- NCDF$var$yearday$varsize      
        
        start2 = c(1, 1, 1)
        varcount = c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2], NCDF$var$yearday$varsize)
        
        #Read in variables:
        #------------------
        lat = get.var.ncdf ( nc=NCDF, varid="lat",                 start = start1, count = latcount )
        lon = get.var.ncdf ( nc=NCDF, varid="lon",                 start = start1, count = loncount )
        dOY = get.var.ncdf ( nc=NCDF, varid="yearday",             start = 1,      count = YDcount  )
        var = get.var.ncdf ( nc=NCDF, varid= paste0(variables[j]), start = start2, count = varcount )
        
        close.ncdf(NCDF)
        
        dOY <- dOY + 1  #Daymet doy starts at 0.
        
        if (year == begYr){
          coords <- cbind( as.vector(lon), as.vector(lat))
          
          distances <- spDistsN1(coords, c(siteLon, siteLat), longlat = TRUE)
          
          minDist <- min(distances)
          
          distpos <- which(distances == minDist)[1]
          
          coords <- as.data.frame(coords) # for indexing
          varLon  <- coords[distpos, 1]
          varLat  <- coords[distpos, 2]
          
          position <- which(lat == varLat & lon == varLon, arr.in = TRUE)  
        }
        
        tempVar <- data.frame(year, dOY, var[position[1], position[2], 1:365])
        names(tempVar) <- c("year", "dOY", paste0(variables[j]))
        
        if ( year == begYr ) ( mainVar <- tempVar)
        if ( year >  begYr ) ( mainVar <- rbind(mainVar, tempVar))
        
        rm(tempVar)
        
      }# End year loop
      
      # Join the variables into one dataframe
      if (j == 1) {allVars <- mainVar} else {allVars <- merge(allVars, mainVar, by = c('year','dOY'), all.x = T)}  
      
    }# End variable loop
    allVars$site <- sites[i]
    tempRecord <- merge(curRecord, allVars, by = c("site", "year", "dOY"), all.x = T, all.y = F, sort = F)
    
    if (i == 1) {fullRecord <- tempRecord} else {fullRecord <- rbind(fullRecord, tempRecord)}
  }
  
  # How long it takes to run
  end.time   <- proc.time()[3]
  print(paste0((end.time-start.time)/3600, " hours"))
  
  fullRecord$airTemp <- (fullRecord$tmin + fullRecord$tmax)/2
  
  return(fullRecord)
}
#=========================================================================================================


#=========================================================================================================
# This function calculates the spatial average of the Daymet variables within a watershed.
# It takes the following:
#    1) The stream temperature record (Site names, latitude, longitude, year, and dOY columns)
#    2) A string of variables to pull from Daymet
#    3) A list of daymet tiles covered by the watersheds
#    4) A master catchments shapefile
#    5) A dataframe of the covariates for the catchments (FEATUREIDs source)
#    6) A list of catchment delineations for the region
# It returns the original dataframe with new columns for the Daymet variables.
#=========================================================================================================
indexUpstreamDaymetVariablesForObservedSites <- function(record, variables, tiles, catchmentShapefile, covariateData, delineatedCatchmentsList){
  
  start.time <- proc.time()[3]
  
  numTiles <- length(tiles)
  
  #------------------------------------------------------------------------------
  # Create a master list of the daymet coords by looping through all of the tiles
  #------------------------------------------------------------------------------
  for (i in 1:numTiles){
    
    # Open the NetCDF 
    NCDF <- open.ncdf(paste0(daymetDir, tiles[i], '_2010/prcp.nc'))    #netcdf   
    
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
    curRecord <- curRecord[which(curRecord$year < 2013),]
    
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





