#=========================================================================================================
# This function creates plots of raw stream temperature data.
# It takes the following:
#    1) The stream temperature record (unique site ID, date, stream temperature, and air temperature columns)
#    2) The directory to save the plots to
#
# It returns a PNG file with plots of air and water plotted against each other and over time.
#=========================================================================================================

plotRawTemperatureData <- function(masterData, plotDirectory){
  
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  
  # Create a list of sites
  sites <- unique(masterData$site)
  
  # Loop through the sites creating plots for each
  for (j in 1:length(sites) ){
    
    print(paste0(round(j/length(sites), digits = 3)*100, '% done.'))
    
    # Index the current site
    curSite <- masterData[masterData$site == sites[j], c('date', 'temp', 'airTemp')]
    
    # Plot of air vs water temperature
    # --------------------------------
    gAirWater <- ggplot( curSite, aes(x = airTemp, y = temp)) + 
      geom_point() +
      ylim(0,30) +
      xlim(-10,30) +
      ylab(expression(paste("Water Temperature ( "[]^{o}, "C )")))+
      xlab(expression(paste("Air Temperature ( "[]^{o}, "C )")))+
      ggtitle('Air/Water Temperature Relationship') +
      theme(title = element_text(face="bold", size=18),
            axis.title.x = element_text(face="bold", size=16),
            axis.text.x = element_text(size=16),
            axis.title.y = element_text(face="bold", size=16),
            axis.text.y = element_text(size=16))
    
    # Plot of water and air temperature over time
    # -------------------------------------------
    # Melt dataframe for use in ggplot
    curSiteMelt <- melt(curSite[,c('date', 'airTemp', 'temp')],id.vars=c('date'))
    
    # Add for ggplotting (can't do this within melt)
    curSiteMelt$year <- as.numeric(strftime(curSiteMelt$date, '%Y'))
    curSiteMelt$dOY <- as.numeric(strftime(curSiteMelt$date, '%j'))
    
    # Plotting
    gTimeSeries <- ggplot( curSiteMelt, aes(x = dOY, y = value, colour = variable)) + 
      scale_color_manual(values=c('red', 'blue')) +
      geom_point() +
      ylim(-10, 30) +
      xlim(0, 366) +
      ylab(expression(paste("Temperature ( "[]^{o}, "C )")))+
      xlab('Day of Year') +
      facet_wrap(~year) + 
      ggtitle('Observed Record') +
      theme(title = element_text(face="bold", size=18),
            axis.title.x = element_text(face="bold", size=16),
            axis.text.x = element_text(size=16),
            axis.title.y = element_text(face="bold", size=16),
            axis.text.y = element_text(size=16))
    
    # Join plots into one object
    gOut <- arrangeGrob( gAirWater, gTimeSeries, ncol=1, main = paste0( 'Site: ', sites[j]))
    
    # Save the plots
    setwd(plotDirectory)
    
    ggsave(plot=gOut, file=paste0(sites[j],'.png'),dpi=300,width=6,height=8, units='in', scale=2)
  }
  
}




#=========================================================================================================
# This function creates a shapefile of the stream temperature data locations.
# It takes the following:
#    1) The stream temperature record (unique site ID, Latitude, Longitude, Agency)
#    2) The projection CRS of the output shapefile.
#
# It returns a points shapefile of the site locations.
#=========================================================================================================
createSiteLocationsShapefile <- function(record, projectionString){
  
  library(sp)
  
  siteLocs <- unique(record[,c('site', 'Latitude', 'Longitude', 'agency')])
  
  sitesShapefile <- SpatialPointsDataFrame(data.frame(siteLocs$Longitude, siteLocs$Latitude), siteLocs, coords.nrs = numeric(0), proj4string = CRS(proj4.NHD), match.ID = TRUE, bbox = NULL)
  
  return(sitesShapefile)
}



