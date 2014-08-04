#######################################################################################################
#                     Functions used in processing raw stream temperaure data
#######################################################################################################


#=========================================================================================================
# This function goes through a timeseries record and removes extra entries(NAs) on the beginning and end of the
#  record for each site in each year. This is mainly a space-saver.
#    1) The stream temperature record (unique site ID, latitude, and longitude columns )
#    2) The name of the column with the data you want to use to trim the dataframe
#
# It returns the same dataframe with extra entries .
#=========================================================================================================

trimNAsFromRecord <- function(record, columnToCheck){
  
  sites <- unique(record$site)
  
  for ( i in 1:length(sites)){
    
    print(i/length(sites))
    
    curSite <- record[record$site == sites[i],]
    
    vals <- which(!is.na(curSite[names(curSite) == columnToCheck]))
    
    trimCur <- curSite[min(vals):max(vals),]
    
    curYears <- unique(trimCur$year)
    
    for( j in 1:length(curYears) ){
      
      curSiteYear <- trimCur[trimCur$year == curYears[j],]
      
      vals <- which(!is.na(curSiteYear[names(curSiteYear) == columnToCheck]))
      
      newRec <- curSiteYear[min(vals): max(vals),]
      
      if( i == 1 & j == 1 ){outRec <- newRec} else( outRec <- rbind(outRec, newRec))
      
      rm(newRec)
    }# End year loop
  }# End site loop
  return(outRec)
}
