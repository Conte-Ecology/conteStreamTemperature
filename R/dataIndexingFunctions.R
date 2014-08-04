# Functions used in working with the different data sources.


#=========================================================================================================
# Description: 
#   This function reads in stream temp timeseries and/or the respective covariate data for different data 
#   sources, joins them together, and outputs a dataframe.
# Usage:
#   readStreamTempData(timeSeries = TRUE, covariates = TRUE, dataSourceList = c('CTDEP', 'MAFW', ...), 
#                      fieldListTS = c('site', 'date', 'temp', ...), 
#                      fieldListCS = 'ALL', 
#                      directory = '.../temperatureProjects/dataIn/)
#
# Arguments:
#    1) timeSeries      A TRUE/FALSE statement of whether to read in the timeseries data.
#    2) covariates      A TRUE/FALSE statement of whether to read in the covariate data.
#    3) dataSourceList  A character vector of the agency abbreviations of the data sources.
#    4) fieldListTS     A character vector of the common fields to be in the output dataframe.
#    5) fieldListCD     A character vector of the common fields to be in the output dataframe. If set 
#                         set equal to 'ALL' then all fields are selected.
#    6) directory       A character vector of the parent dataframe of where the data is stored (dataIn).
#
# It returns a dataframe with the site name, lat/lon, FEATUREID, and the select covariate values.
#=========================================================================================================
readStreamTempData <- function(timeSeries, covariates, dataSourceList, fieldListTS, fieldListCD, directory){
  
  # Loop through the agencies
  for ( i in 1:length(dataSourceList)){
    
    # Read in timeseries data
    # -----------------------
    if( timeSeries == T ){
      
      # Individual files
      load(paste0(directory, dataSourceList[i], '/observedStreamTempAndClimateData_', dataSourceList[i], '.RData'))
      
      # Join files
      if ( i == 1 ) { tS <- masterData[,fieldListTS]} else ( tS <- rbind(tS, masterData[,fieldListTS]) )
    }
    
    # Read in covariates
    # ------------------
    if( covariates == T ){
      
      # Individual files
      load(paste0(directory, dataSourceList[i], '/covariateData_', dataSourceList[i], '.RData'))

      # Index the fields
      if( 'ALL' %in% fieldListCD ) {covs <- covariateData} else( covs <- covariateData[, fieldListCD])
      
      # Join files
      if ( i == 1 ) { cD <- covs} else ( cD <- rbind(cD, covs) )
    }
    
  }
  
  # Prepare output
  ifelse(  exists('tS') &  exists ('cD'),  dataOut <- merge(tS, cD, by = 'site', all.x = T, all.F = F, sort = F), 
           ifelse(  exists('tS') & !exists ('cD'),  dataOut <- tS,
                    ifelse( !exists('tS') &  exists ('cD'),  dataOut <- cD, 
                            print("'timeSeries' or 'covariates' must be TRUE."))))
  
  # Output
  return(dataOut)
}
#=========================================================================================================


#=========================================================================================================
# This function indexes values from the master list of covariates for observed stream temperature sites.
# It takes the following:
#    1) The stream temperature record (unique site ID, latitude, and longitude columns )
#    2) A dataframe of the covariates for the catchments (FEATUREIDs source)
#    3) A master catchments shapefile
#    4) A CRS string of the spatial data projection
#    5) A string of variables to pull from the covariates list
#
# It returns a dataframe with the site name, lat/lon, FEATUREID, and the select covariate values.
#=========================================================================================================
indexCovariateData <- function(record, masterCovariates, catchmentShapefile, projectionString, fields){
  start.time <- proc.time()[3]
  
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
  
  newCovs$site <- paste(newCovs$site)
  
  # Return the covariates list
  return(newCovs)
  
  # How long it takes to run
  end.time   <- proc.time()[3]
  print(paste0((end.time-start.time)/3600, " hours"))
}
#=========================================================================================================



#=========================================================================================================
# This function corrects the covariate data file after the site locations have been manually checked.
# It takes the following:
#    1) The original covariate data file
#    2) The "siteChanges" CSV file
#    3) The master dataframes of both local and upstream covariate statistics
#    4) A list of layers that get NAs assigned for local catchment values. (This will hopefully change with an updated layer)
# It returns the same covariate dataframe with corrected values and a column indicating whether or not values
#   for that site changed.
#=========================================================================================================
correctCovariateData <- function(covariateData, siteChanges, LocalStats, UpstreamStats, impoundmentLayers){
  
  d <- covariateData
  s <- siteChanges
  
  # Correct the factor problem
  s$site <- as.character(s$site)
  d$site <- as.character(d$site)
  
  # Select sites that get changed
  s1 <- s[which(s$correctFeatureID > 1 | s$localCatchment > 0 ),]
  
  # Edit the "correctFeatureID" column to contain the correct FEATUREIDs (both changes and ones that will stay the same).
  s1$correctFeatureID[s1$correctFeatureID == 1] <- s1$currentFeatureID[s1$correctFeatureID == 1]
  
  # Don't want to replace Lat/Lon of the site with Lat/Lon of catchment centroid, so remove these from the master list
  up  <- UpstreamStats[, - which(names(UpstreamStats) %in% c('Latitude', 'Longitude'))]
  loc <- LocalStats   [, - which(names(LocalStats) %in% c('Latitude', 'Longitude'))]
  
  # Remove covariate data for sites that will be replaced
  d1 <- d[!(d$site %in% s1$site), ]
  
  #Replace sites that need new Upstream covariate data:
  #====================================================
  
  # Sites that get upstream data:
  # -----------------------------
  # First, check if any sites get upstream data
  if(length(which(s1$localCatchment == 0)) > 0 ){
  
    nU <- s1[s1$localCatchment == 0, c('site', 'correctFeatureID')]
    colnames(nU) <- c('site', 'FEATUREID')
    nU1 <- merge(nU, d[,c('site', 'Latitude', 'Longitude')], by = 'site', all.x = T, sort = F)
    
    # Merge in new covariate data
    nU2 <- merge(nU1, up, by = 'FEATUREID', all.x = T, sort = F)
    
    # Remove unused columns
    nU3 <- nU2[,names(nU2) %in% names(d1)]
  }
  
  # Sites that get local data:
  # --------------------------
  # First, check if any sites get local data
  if(length(which(s1$localCatchment == 1)) > 0 ){
  
    nL <- s1[s1$localCatchment == 1, c('site', 'correctFeatureID')]
    colnames(nL) <- c('site', 'FEATUREID')
    nL1 <- merge(nL, d[,c('site', 'Latitude', 'Longitude')], by = 'site', all.x = T, sort = F)
    
    # Merge in new covariate data
    nL2 <- merge(nL1, loc, by = 'FEATUREID', all.x = T, sort = F)
    
    # Remove unused columns
    nL3 <- nL2[,names(nL2) %in% names(d1)]
    
    # Replace local stats for impoundments with NAs because these values are not applicable to the local catchments (they are only used for upstream calculation)
    nL3[, names(nL3) %in% impoundmentLayers] <- NA
  } 
  
  # Join new data together:
  # -----------------------
  ifelse(exists("nU3") & exists("nL3"), newData <- rbind(nU3, nL3),
  ifelse(exists("nU3") & !exists("nL3"), newData <- nU3,
  ifelse(!exists("nU3") & exists("nL3"), newData <- nL3, 
                noSites <- "Y")))
  
  if( exists("noSites") ) { print( "There are no sites to replace.")}
  
  if( !exists("noSites")){
    # Add a column so we know which sites have had their data changed
    newData$locationChange <- 'Y'
    d1$locationChange <- 'N'
    
    # Join new to existing keeper data
    dataOut <- rbind(d1, newData)
    
    return(dataOut)} else( return(covariateData))
  
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
          "Tile Error")))))))))))))

  return(Tile)
}
# ** Explicitly checked with mapping software (ArcGIS).
#=========================================================================================================


#=========================================================================================================
#
# Find Nearest Daymet Point
#
# Description: 
#   This function takes the latitude and longitude of a site and calculates the nearest Daymet point based
#     on the latitude/longitude arrays read in from Daymet. 
# 
# Usage:
#   findNearestDaymetPoint(siteLat = 44.0 , siteLon = -72.5, 
#                           dayLat = lat, dayLon = lon, 
#                           currentTile = 11934)
#
# Arguments:
#    1) siteLat         A single latitude value of the site.
#    2) siteLon         A single longitude value of the site.
#    3) dayLat          An array of latitude points as read in from Daymet NetCDF file.
#    4) dayLon          An array of longitude points as read in from Daymet NetCDF file.
#    5) currentTile     The Daymet tile that is currently being indexed.
#
# It returns the array position (row, col) of the nearest Daymet point to the site. This is used to index
#   the climate records from the variable array as read in from Daymet.
#=========================================================================================================
findNearestDaymetPoint <- function(siteLat, siteLon, dayLat, dayLon, currentTile){
  
  library(sp)
  
  # Calculate distances between site and Daymet point
  coords <- cbind( as.vector(dayLon), as.vector(dayLat))
  coords <- data.frame(Lon = coords[,1], Lat = coords[,2], dist = spDistsN1(coords, c(siteLon, siteLat), longlat = TRUE))
  
  # Index lat/lon of closest points
  nearPts <-  coords[coords$dist == min(coords$dist), c('Lat', 'Lon')]
  
  # Sometimes duplicate points get selected. Pick only one and make sure it's in your current tile to avoid NAs
  nearPts$Tile <- indexDaymetTileByLatLon(nearPts$Lat, nearPts$Lon)
  
  # Final coords
  varLat <- nearPts$Lat[nearPts$Tile == currentTile][1]
  varLon <- nearPts$Lon[nearPts$Tile == currentTile][1]
  
  # Determine the position of the record in the Daymet array
  position <- which(dayLat == varLat & dayLon == varLon, arr.in = TRUE)
  
  # Return the position in the Daymet array
  return(position)
}# End function
#=========================================================================================================


#=========================================================================================================
#
# Pair Daymet Climate Data for Local Points
#
# Description: 
#   This function takes a period of record for a location and returns the corresponding time series of 
#     climate data from Daymet. This is done by indexing the nearest Daymet point to the site location.
# 
# Usage:
#   indexLocalDaymetVariablesForObservedSites(record = masterData , 
#       variables = c('prcp', 'tmin', 'tmax'), 
#       daymetDirectory = '//IGSAGBEBWS-MJO7/projects/dataIn/environmental/climate/daymet/unzipped/Daily')
#
# Arguments:
#    1) record          A dataframe indicating the period of record ( limited to 1980 - 2013) over which to 
#                         pull climate data. Minimum dataframe column requirements and format:
#                         ---------------------------------------------------------------------------------
#                         'data.frame':  3228940 obs. of  5 variables:
#                          $ site     : chr  "MADEP_W1466_M1" "MADEP_W1466_M1" "MADEP_W1466_M1" ...
#                          $ year     : num  1980 1980 1980 1980 1980 1980 1980 1980 1980 ...
#                          $ dOY      : num  1 2 3 4 5 6 7 8 9 ...
#                          $ Latitude : num  42.5 42.5 42.5 42.5 ...
#                          $ Longitude: num  -72.9 -72.9 -72.9 -72.9 ...
#                         ---------------------------------------------------------------------------------
#    2) variables       A string listing the climate variables (as named by Daymet) to be indexed
#    3) daymetDirectory A string indicating the folder where Daymet files reside
#
# It returns the record with added columns for Daymet variable records.
#=========================================================================================================
indexLocalDaymetVariablesForObservedSites <- function(record, variables, daymetDirectory){
  
  library(ncdf)
  library(sp)
  
  # Add tile assignments
  record$tile <- indexDaymetTileByLatLon(record$Latitude, record$Longitude)
  
  # Generate tile list
  tiles <- unique(record$tile)
  
  # Tile loop
  for (tile in tiles){
    
    # Index records of sites in tile
    tileRecord <- record[record$tile == tile,]
    
    # Read in a sample NetCDF of the tile to determine site positions within the Daymet array
    # ---------------------------------------------------------------------------------------
    # Open NetCDF
    locNCDF <- open.ncdf(paste0(daymetDirectory, tile, '_1980/prcp.nc'))    #netcdf
    
    #Dimension limits of each of the variables we'll use:
    start1 = c(1,1)
    latcount <- c(locNCDF$var$lat$varsize[1], locNCDF$var$lat$varsize[2])
    loncount <- c(locNCDF$var$lon$varsize[1], locNCDF$var$lon$varsize[2])
    
    #Read in variables:
    lat = get.var.ncdf ( nc=locNCDF, varid="lat", start = start1, count = latcount )
    lon = get.var.ncdf ( nc=locNCDF, varid="lon", start = start1, count = loncount )
    
    # Close connection
    close.ncdf(locNCDF)
    
    # Create dataframe to store location info
    meta <- unique(tileRecord[,c('site', 'Latitude', 'Longitude')])
    meta$row <- NA; meta$col <- NA
    
    # Determine the position in Daymet array of the site's record.
    #   This should hold true for all variables and years within the tile.
    for( i in 1:nrow(meta)){
      
      position <- findNearestDaymetPoint(meta$Latitude[i], meta$Longitude[i], lat, lon, tile)
      
      meta$row[i] <- position[,'row']
      meta$col[i] <- position[,'col']  
    }
    
    # Index years in current record
    years <- unique(tileRecord$year)[order(unique(tileRecord$year))]
    
    # Year loop
    for (year in years){
      
      # Clip the record to the sites within the current tile
      tileYearRecord <- tileRecord[tileRecord$year == year,]
      
      # Print status
      print(paste(tile, '-', year))
      
      # List of current sites
      tileYearSites <- unique(tileYearRecord$site)
      
      # Pre-allocate columns for Daymet Variables
      tileYearDaymet <- as.data.frame(matrix(data = NA, nrow = 365 * length(tileYearSites), ncol = 3 + length(variables)))
      names(tileYearDaymet) <- c('site', 'year', 'dOY', variables)
      
      # Fill year column
      tileYearDaymet$year <- year
      
      # Varialble loop
      for (variable in variables){
        
        # Open the NetCDF of the current tile, year, and variable
        NCDF <- open.ncdf(paste0(daymetDirectory, tile, '_', year,'/', variable, '.nc'))    #netcdf
        
        # Dimension limits of each of the objects used
        start1 = c(1,1)
        YDcount  <- NCDF$var$yearday$varsize
        
        start2 = c(1, 1, 1)
        varcount = c(NCDF$var$lat$varsize[1], NCDF$var$lat$varsize[2], NCDF$var$yearday$varsize)
        
        # Read in objects
        dOY = get.var.ncdf ( nc=NCDF, varid="yearday",             start = 1,      count = YDcount  )
        var = get.var.ncdf ( nc=NCDF, varid= variable, start = start2, count = varcount )
        
        # Close the NetCDF
        close.ncdf(NCDF)
        
        #Daymet doy starts at 0
        dOY <- dOY + 1  
        
        # Loop through sites within tile and assign record
        for ( s in 1:length(tileYearSites)){
          
          # Index array position
          row <- meta$row[meta$site == tileYearSites[s]]
          col <- meta$col[meta$site == tileYearSites[s]]
          
          # Fill sites, dOY, and variable columns
          tileYearDaymet$site[dOY+(s-1)*365] <- tileYearSites[s]
          tileYearDaymet$dOY[dOY+(s-1)*365]  <- dOY
          tileYearDaymet[dOY+(s-1)*365, which(names(tileYearDaymet) == variable) ] <- var[row, col, 1:365]
        }# End sites loop
      }# End variable loop
      
      if( ! exists("tileDaymet") ){ tileDaymet <- tileYearDaymet } else(tileDaymet <- rbind(tileDaymet, tileYearDaymet) )
    }# End year loop
    
    if( ! exists("daymet") ){ daymet <- tileDaymet } else(daymet <- rbind(daymet, tileDaymet) )
    
    rm(tileDaymet); gc()
    
    print( (proc.time()[3] - start.time)/3600 )
  }# End tile loop
  
  outRecord <- merge(record, daymet, by = c('site', 'year', 'dOY'), all.x = T, all.y = F, sort = F)
  
  return(outRecord)
}# End function
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





