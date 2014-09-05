#' @title readStreamTempData
#'
#' @description
#' \code{readStreamTempData} reads in stream temp timeseries and/or the respective covariate data for different data sources, joins them together, and outputs a dataframe
#'
#' @param timeSeries A TRUE/FALSE statement of whether to read in the timeseries data.
#' @param covariates A TRUE/FALSE statement of whether to read in the covariate data.
#' @param dataSourceList A character vector of the agency abbreviations of the data sources.
#' @param fieldListTS A character vector of the common fields to be in the output dataframe.
#' @param fieldListCD A character vector of the common fields to be in the output dataframe. If set equal to 'ALL' then all fields are selected.
#' @param directory A character vector of the parent dataframe of where the data is stored (dataIn).
#' @return Returns a dataframe with the site name, lat/lon, FEATUREID, and the select covariate values
#' @details
#' This function reads in stream temp timeseries and/or the respective covariate data for different data sources, joins them together, and outputs a dataframe.
#' @export
readStreamTempData <- function(timeSeries, covariates, dataSourceList, fieldListTS, fieldListCD, directory){
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


#' @title assignCatchments
#'
#' @description
#' \code{assignCatchments} Assigns a catchment to a site with a set of lat/lon points.
#'
#' @param sites dataframe containing 3 columns in specified order: 1)SiteID 2)Longitude 3)Latitude
#' @param catchmentShapefile SpatialPolygonsDataframe of the catchment shapefile over which the points will be matched
#' @param catchmentID a character vector of column name describing the catchment identifier in the shapefile (e.g. "FEATUREID" for NHDplusV2)
#' @param projectionString a CRS string of the spatial data projection of the shapefile and site coordinates
#' @return Returns the input dataframe with the catchment ID appended (column 4)
#' @details
#' This function uses a spatial overlay to assign catchment IDs to sites with associated lat/lon points. If a point does not match any catchment, NA is returned.
#' @export
assignCatchments <- function(sites, catchmentShapefile, catchmentID, projectionString){
  
  require(maptools)
  require(sp)
  
  # Save original names
  siteCol <- names(sites)[1]
  LonCol  <- names(sites)[2]
  LatCol  <- names(sites)[3]
  
  # Assign names for processing
  names(sites) <- c('siteID', 'Longitude', 'Latitude')
  
  # Convert lat/lon to format for overlay(s.p.d.f.)
  points <- SpatialPointsDataFrame(as.matrix(sites[,c('Longitude', 'Latitude')]), sites, proj4string = projectionString)
  
  # Overlay the points and return catchments to dataframe
  sitesOut <- data.frame(points@data, over(points,catchmentShapefile)[,c(catchmentID)])
  
  names(sitesOut) <- c(siteCol, LonCol, LatCol, catchmentID)

  # Return the covariates list
  return(sitesOut)
}

#' @title indexCovariateData
#'
#' @description
#' \code{indexCovariateData} Indexes values from the master list of covariates given either a list of catchments or a list of latitude and longitude values
#'
#' @param method a character vector describing the indexing method. Either "catchments" or "coordinates". Default is "catchments".
#' @param sites dataframe with the indexing information. If method = "catchments" then sites = a vector of catchment IDs as described by the "catchmentID" parameter. If method = "coordinates" then a dataframe containing 3 columns in specified order: 1)SiteID 2)Longitude 3)Latitude.
#' @param masterCovariates a dataframe of the covariates for the catchments IDs
#' @param fields a string of variables to pull from the covariates list. Deafault = ALL
#' @param catchmentID a character vector of column name describing the catchment identifier in the shapefile (e.g. "FEATUREID" for NHDplusV2)
#' @param catchmentShapefile required if method = "coordinates". SpatialPolygonsDataframe of the catchment shapefile over which the points will be matched
#' @param projectionString required if method = "coordinates". a CRS string of the spatial data projection of the shapefile and site coordinates
#' @return Returns a dataframe with the input data and the select covariate value fields
#' @details
#' This function indexes values from the master list of covariates for observed stream temperature sites.
#' If "Latitude" or "Longitude" exist within the "fields" and the "coordinates" method is used then the indexed values will be renamed to "catchmentLatitude" and "catchmentLongitude"
#' If fields are indicated that are not in the master covariate object, then they are returned as all NA.
#' @export
indexCovariateData <- function(method = "coordinates", sites, masterCovariates, fields = NULL,  catchmentID , catchmentShapefile = NULL, projectionString = NULL){
  
  # fields default
  if( is.null(fields) ) {fields <- names(masterCovariates) }
  
  # "coordinates" method assigns catchment IDs
  if( method == 'coordinates' ) { 
    
    if (is.null(catchmentShapefile)) {print("Error: missing catchmentShapefile")}
    if (is.null(projectionString))   {print("Error: missing projectionString")}
    
    sites <- assignCatchments(sites, catchmentShapefile, catchmentID, projectionString)
  
    featureIDs <- sites[,4]
  }
  
  # "catchments" method simply lists catchments
  if( method == 'catchments' ) { featureIDs <- sites }
  
  # Select fields
  selectUpstreamStats <- masterCovariates[masterCovariates[,names(masterCovariates)== catchmentID] %in% featureIDs, names(masterCovariates) %in% fields]  
  
  if(!all(fields %in% names(masterCovariates))){
    
    empty <- fields[!fields %in% names(masterCovariates)]
    emptyStats <- data.frame(matrix(data = NA, ncol = length(empty), nrow = nrow(selectUpstreamStats)))
    names(emptyStats) <- empty
  
    selectUpstreamStats <- data.frame(selectUpstreamStats, emptyStats)
    
    print(paste("Warning: Some fields do not exist in the master list. Returning NA in the following columns:", empty))
  }
  
  # If the "coordinates" method is used then rename the catchment centroid coordinates:
  if( method == 'coordinates' ) {
    if( 'Latitude'  %in% names(selectUpstreamStats) ) { names(selectUpstreamStats)[names(selectUpstreamStats) == 'Latitude' ] <- 'catchmentLatitude'  }
    if( 'Longitude' %in% names(selectUpstreamStats) ) { names(selectUpstreamStats)[names(selectUpstreamStats) == 'Longitude'] <- 'catchmentLongitude' }    
  
    newCovs <- merge(sites, selectUpstreamStats, by = catchmentID, all.x = T, all.y = F, sort = F)
  }
    
  # If "catchments" method is used then just retrun the indexed covariates.
  if( method == 'catchments') {newCovs <- selectUpstreamStats}
  
  # Return the covariates list
  return(newCovs)
}



#' @title indexDaymetTileByLatLon
#'
#' @description
#' \code{indexDaymetTileByLatLon} This function indexes Daymet tiles by the latitude and longitude of a site.
#'
#' @param SiteLat Site latitude
#' @param SiteLon Site longitude
#' 
#' @return Returns the Daymet tile that holds data for the coordinates given.
#' @details
#' This function indexes Daymet tiles by the latitude and longitude of a site.
#' @export
indexDaymetTileByLatLon <- function(SiteLat, SiteLon){

  Tile <- ifelse( SiteLat > 38 & SiteLat < 40 & SiteLon > -80 & SiteLon < -78, 11571,
          ifelse( SiteLat > 36 & SiteLat < 38 & SiteLon > -80 & SiteLon < -78, 11391,
          ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -74 & SiteLon < -72, 11754,
          ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -72 & SiteLon < -70, 11755,
          ifelse( SiteLat > 40 & SiteLat < 42 & SiteLon > -70 & SiteLon < -68, 11756,
          ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -74 & SiteLon < -72, 11934,
          ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -72 & SiteLon < -70, 11935,
          ifelse( SiteLat > 42 & SiteLat < 44 & SiteLon > -70 & SiteLon < -68, 11936,
          ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -74 & SiteLon < -72, 12114,  
          ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -72 & SiteLon < -70, 12115,
          ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -70 & SiteLon < -68, 12116,
          ifelse( SiteLat > 44 & SiteLat < 46 & SiteLon > -68 & SiteLon < -66, 12117,
          ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -72 & SiteLon < -70, 12295,     
          ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -70 & SiteLon < -68, 12296,     
          ifelse( SiteLat > 46 & SiteLat < 48 & SiteLon > -68 & SiteLon < -66, 12297,   
          "Tile Error")))))))))))))))

  return(Tile)
}


#' @title findNearestDaymetPoint
#'
#' @description
#' \code{findNearestDaymetPoint} This function takes the latitude and longitude of a site and calculates the nearest Daymet point based on the latitude/longitude arrays read in from Daymet. 
#'
#' @param siteLat A single latitude value of the site.
#' @param siteLon A single longitude value of the site.
#' @param dayLat An array of latitude points as read in from Daymet NetCDF file.
#' @param dayLon An array of longitude points as read in from Daymet NetCDF file.
#' @param currentTile The Daymet tile that is currently being indexed.
#' 
#' @return Returns array position (row, col) of the nearest Daymet point to the site.
#' @details
#' This function takes the latitude and longitude of a site and calculates the nearest Daymet point based on the latitude/longitude arrays read in from Daymet. It returns the array position (row, col) of the nearest Daymet point to the site. This is used to index the climate records from the variable array as read in from Daymet.
#' 
#' @examples
#' 
#' \dontrun{
#' findNearestDaymetPoint(siteLat = 44.0 , siteLon = -72.5, dayLat = lat, dayLon = lon, currentTile = 11934)
#' }
#' @export
findNearestDaymetPoint <- function(siteLat, siteLon, dayLat, dayLon, currentTile){
  
  require(sp)
  
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
}


#' @title indexLocalDaymetVariables
#'
#' @description
#' \code{indexLocalDaymetVariablesForObservedSites} takes a period of record for a location and returns the corresponding time series of climate data from Daymet. This is done by indexing the nearest Daymet point to the site location.
#'
#' @param record dataframe indicating the period of record ( limited to 1980 - 2013) over which to pull climate data. Dataframe columns must be in the following order: 1)Site ID 2)Year 3)Day of Year 4)Longitude 5)Latitude
#' @param variables character string listing the climate variables (as named by Daymet) to be indexed
#' @param daymetDirectory character string indicating the folder where Daymet files reside
#' 
#' @return Returns the record with added columns for specified Daymet variable records.
#' @details
#' This function takes a period of record for a location and returns the corresponding time series of climate data from Daymet. This is done by indexing the nearest Daymet point to the site location.
#' Input names for "record" can vary as long as order is correct. Original names will be included in output.
#' 
#' @examples
#' 
#' \dontrun{
#' indexLocalDaymetVariables(record = masterData , variables = c('prcp', 'tmin', 'tmax'), daymetDirectory = '//IGSAGBEBWS-MJO7/projects/dataIn/environmental/climate/daymet/unzipped/Daily')
#' }
#' @export
indexLocalDaymetVariables <- function(record, variables, daymetDirectory){
  
  require(ncdf)
  require(sp)
  
  # Save original names
  siteCol <- names(record)[1]
  yearCol <- names(record)[2]
  dOYCol  <- names(record)[3]
  LonCol  <- names(record)[4]
  LatCol  <- names(record)[5]
  
  names(record) <- c('site', 'year', 'dOY', 'Longitude', 'Latitude')
  
  # Add tile assignments
  record$tile <- indexDaymetTileByLatLon(record$Latitude, record$Longitude)
  
  # Generate tile list
  tiles <- unique(record$tile)
  
  # In case this object exists
  if( exists("daymet") ) { rm(daymet) }
  
  # Tile loop
  for (tile in tiles){
    
    # Index records of sites in tile
    tileRecord <- record[record$tile == tile,]
    
    # Read in a sample NetCDF of the tile to determine site positions within the Daymet array
    # ---------------------------------------------------------------------------------------
    # Open NetCDF    
    locNCDF <- open.ncdf(file.path(daymetDirectory, paste0(tile, '_1980/prcp.nc')))    #netcdf
    
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
        NCDF <- open.ncdf( file.path(daymetDirectory, paste0(tile, '_', year), paste0(variable, '.nc') ) )   #netcdf
                
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
    
  }# End tile loop
  
  outRecord <- merge(record, daymet, by = c('site', 'year', 'dOY'), all.x = T, all.y = F, sort = F)
  
  # Rename output
  names(outRecord)[names(outRecord) == 'site']      <- siteCol
  names(outRecord)[names(outRecord) == 'year']      <- yearCol
  names(outRecord)[names(outRecord) == 'dOY']       <- dOYCol
  names(outRecord)[names(outRecord) == 'Longitude'] <- LonCol
  names(outRecord)[names(outRecord) == 'Latitude']  <- LatCol
  
  # Remove "tile" column
  outRecord <- outRecord[,-which(names(outRecord) == "tile")]
  
  return(outRecord)
}


#' @title indexUpstreamDaymetVariablesByCatchment
#'
#' @description
#' \code{indexUpstreamDaymetVariablesByCatchment} calculates the spatial average of the Daymet variables within a watershed.
#'
#' @param record dataframe indicating the period of record (limited to 1980 - 2013) over which to pull climate data. Columns must be in the following order: 1)Catchment ID 2)Year 3)Day of Year 4)Site ID (if applicable)
#' @param variables A string of variables to pull from Daymet. Options: "prcp", "tmin", "tmax", "srad", "vp", "swe", & "dayl"
#' @param catchmentShapefile SpatialPolygonsDataframe of the catchment shapefile covering the area over which the sites will be delineated
#' @param catchmentID a character vector of column name describing the catchment identifier in the shapefile (e.g. "FEATUREID" for NHDplusV2)
#' @param delineatedCatchmentsList A list of catchment delineations for the region
#' @param daymetDirectory Directory
#' 
#' @return Returns the original dataframe with new columns for the Daymet variables.
#' @details
#' This function calculates the spatial average of the Daymet variables within a watershed.
#' @export
indexUpstreamDaymetVariablesByCatchment <- function(record, variables, catchmentShapefile, catchmentID, delineatedCatchmentsList, daymetDirectory){
  
  # Libraries
  require(ncdf)
  require(rgeos)
  require(sp)
  require(reshape2)
  
  # Save original names
  yearCol  <- names(record)[2]
  dOYCol   <- names(record)[3]
  if(ncol(record) > 3) {siteCol  <- names(record)[4]}
  
  # Temporarily rename record
  names(record) <- c(catchmentID, 'year', 'dOY', 'site')

  # Get list of tiles and years over which data exists
  sourceInfo <- colsplit(list.files(daymetDirectory), "_", c('tiles', 'years'))
  
  tiles <- unique(sourceInfo$tiles)
  minDayYr <- min(sourceInfo$years)
  maxDayYr <- max(sourceInfo$years)
  
  #------------------------------------------------------------------------------
  # Create a master list of the daymet coords by looping through all of the tiles
  #------------------------------------------------------------------------------
  for (i in seq_along(tiles)){
    
    # Open the NetCDF     
    NCDF <- open.ncdf(file.path(daymetDirectory, paste0(tiles[i], '_1980'), 'prcp.nc'))    #netcdf
    
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
  # Loop through catchmentIDs and NetCDFs, getting data.
  #----------------------------------------------  
  
  fids <- unique(record[,c(catchmentID)])
  
  masterLength <- length(delineatedCatchmentsList)
  
  for ( i in seq_along(fids) ){
    
    print(paste0(round(i/length(fids), digits = 3)*100, '% done.'))
    
    # Define the catchment polygon:
    #------------------------------
    featureID <- fids[i] 
    features <- delineatedCatchmentsList[[which(sapply(c(1:masterLength),FUN = function(x){delineatedCatchmentsList[[x]][1]==featureID})==TRUE)]]
    
    catchmentShape <- catchmentShapefile[catchmentShapefile$FEATUREID %in% features,]
    basinShape     <- gUnaryUnion(catchmentShape) #dissolve individual catchments
    a <- SpatialPoints(masterCoords, proj4string=CRS(proj4.NHD))
    
    inside <- as.data.frame(a[!is.na(over(a, basinShape)),])
    
    #If no point falls within the catchment, find the nearest one:
    #-------------------------------------------------------------
    if(nrow(inside) == 0 ){
      
      tempLat <- coordinates(basinShape)[,2]
      tempLon <- coordinates(basinShape)[,1]
      
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
    
    # Dataframe of coordinates and which tiles they fall into
    spatialLocs <- spatialLocs[ order(spatialLocs$Tile), ]
    subTiles <- unique(spatialLocs$Tile)
    
    
    #Link record with catchment area:
    #--------------------------------
    
    # Pull records for current catchmentID in range of years with data
    curRecord <- record[which(record[,c(catchmentID)] == fids[i]),c(catchmentID, 'year', 'dOY')]
    
    # Years in the record
    recYears <- unique(curRecord$year)
    
    # Don't loop over years with no Daymet data
    loopYears <- recYears[which(recYears >= minDayYr & recYears <= maxDayYr)]
    
    # Order
    loopYears <- loopYears[order(loopYears)]
    
    # First year
    begYr <- min(loopYears)
    
    # Loop through daymet variables
    for (j in 1:length(variables)){
      
      # Loop through record years
      for ( year in loopYears ){
        
        # Loop through tiles
        for ( t in 1:length(subTiles)){
                    
          # Open NetCDF
          NCDF <- open.ncdf(file.path(daymetDirectory, paste0(subTiles[t], '_1980'), paste0(variables[j], '.nc') ) )    #netcdf
          
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
          
          # Which coordinates are in the current tile
          xx <- spatialLocs[which(spatialLocs$Tile == subTiles[t]),]
          
          # Index data for points within the watershed
          for ( m in 1:nrow(xx) ){
            position <- which(lon == xx$Longitude[m] & lat == xx$Latitude[m], arr.in = TRUE)
            
            if ( t == 1 & m == 1) {tempVar <- data.frame(year, dOY, var[position[1], position[2], 1:365])} else(tempVar <- cbind(tempVar, var[position[1], position[2], 1:365]))
          }
          
        } # end tile loop
        
        # Take the average across all points
        ifelse( ncol(tempVar) > 3, R <- rowMeans(tempVar[,-c(1,2)], na.rm = TRUE, dims = 1),  R <- tempVar[,-c(1,2)] )
        
        # Replace with means
        tempVar <- data.frame(tempVar[,c(1,2)], R)
        
        # Name columns
        names(tempVar) <- c("year", "dOY", paste0(variables[j]))
        
        # Store spatial averages for current variable
        if ( year == begYr ) ( mainVar <- tempVar)
        if ( year >  begYr ) ( mainVar <- rbind(mainVar, tempVar))
        
        rm(tempVar, R)
      } # end year loop
      
      # Store spatial average for all variables
      if (j == 1) {allVars <- mainVar} else {allVars <- merge(allVars, mainVar, by = c('year','dOY'), all.x = T)}  
    } # end variable loop
    
    #Add data into main dataframe
    allVars[ ,paste(catchmentID)] <- fids[i]
    
    if (i == 1) {tempRecord <- allVars} else {tempRecord <- rbind(tempRecord, allVars)}
  }
  
  # Merge into original dataframe by catchmentID
  outRecord <- merge(record, tempRecord, by = c(catchmentID, "year", "dOY"), all.x = T, all.y = F, sort = F)
  
  # Rename output
  names(outRecord)[names(outRecord) == 'year']      <- yearCol
  names(outRecord)[names(outRecord) == 'dOY']       <- dOYCol
  if( exists('siteCol') ) { names(outRecord)[names(outRecord) == 'site'] <- siteCol}
  
  # Return record with Daymet variables
  return(outRecord)
}
