#' @title deriveMetricsYear
#'
#' @description
#' \code{deriveMetricsYear} Wrapper to calculate derived metrics in parallel using foreach %dopar%
#'
#' @param data Dataframe created with the prepData and predictTemp functions
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' blah, blah, blah, something, something
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
deriveMetricsYear <- function(data) {
  
  library(dplyr)
  library(tidyr)
  
  byfeatureid <- group_by(data, featureid)
  byfeatureidYear <- group_by(byfeatureid, year, add = TRUE)
  #(maxTempfeatureid <- dplyr::dplyr::summarise(byfeatureid, max(tempPredicted, na.rm = T)))
  
  derivedfeatureidMetrics <- byfeatureidYear %>%
    distinct
  
  # Mean maximum daily mean temperature by featureid (over years) - this must be calculated before totalObs to work with NA and left join if all data missing
  meanMaxTemp <- byfeatureidYear %>%
    dplyr::summarise(maxTemp = max(tempPredicted, na.rm = T)) 
  
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanMaxTemp, by = c("featureid", "year"))
  rm(meanMaxTemp)
  
  # total observations (days with data) per featureid
  totalObs <- byfeatureidYear %>%
    dplyr::filter(!is.na(temp)) %>%
    dplyr::summarise(totalObs = n()) 
  
  
  derivedfeatureidMetrics <- dplyr::left_join(derivedfeatureidMetrics, totalObs, by = c("featureid", "year")) %>%
    dplyr::mutate(totalObs = as.numeric(ifelse(is.na(totalObs), 0, as.numeric(totalObs))))
  derivedfeatureidMetrics <- dplyr::select(derivedfeatureidMetrics, featureid, year, totalObs, maxTemp)
  
  # Mean July Air Temp
  meanJulyAir <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month == 7) %>%
    dplyr::summarise(JulyAirTemp = mean(airTemp, na.rm = T))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanJulyAir, by = c("featureid", "year"))
  rm(meanJulyAir)
  
  # Mean July Water Temp
  meanJulyTemp <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month == 7) %>%
    dplyr::summarise(JulyTemp = mean(tempPredicted, na.rm = T))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanJulyTemp, by = c("featureid", "year"))
  rm(meanJulyTemp)
  
  # Mean July Obs Temp
  # meanJulyObs <- byfeatureidYear %>%
  #   dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
  #   dplyr::filter(month == 7)
  # 
  # nJulyObs <- meanJulyObs %>%
  #   dplyr::filter(!is.na(temp)) %>%
  #   dplyr::summarise(nJulyObs = n())
  # 
  # nJulyObsYears <- nJulyObs %>%
  #   dplyr::filter(nJulyObs >= 25) %>% # at least 25 days with observed data
  #   dplyr::summarise(nYearsJulyObs = n())
  # 
  # meanJulyObs <- nJulyObs %>%
  #   left_join(meanJulyObs) %>%
  #   dplyr::filter(nJulyObs >= 25) %>%
  #   dplyr::summarise(meanJulyObs = mean(temp, na.rm = T)) 
  # 
  # derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
  #   dplyr::left_join(meanJulyObs) %>%
  #   dplyr::left_join(nJulyObsYears)
  # rm(meanJulyObs)
  
  # Mean Aug Temp
  meanAugTemp <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month == 8) %>%
    dplyr::summarise(AugTemp = mean(tempPredicted, na.rm = T))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanAugTemp, by = c("featureid", "year"))
  rm(meanAugTemp)
  
  # Mean Summer Temp
  meanSummerTemp <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month >= 6 & month <= 8) %>%
    dplyr::summarise(SummerTemp = mean(tempPredicted, na.rm = T))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanSummerTemp, by = c("featureid", "year"))
  rm(meanSummerTemp)
  
  # Annual mean 30-day maximum mean daily temperature (#------- only works if 30 days of data otherwise throws error------------)
  # mean30Day <- byfeatureidYear %>%
  #   dplyr::arrange(featureid, year, dOY) %>%
  #   dplyr::mutate(mm30Day = rollapply(data = tempPredicted, 
  #                                     width = 30, 
  #                                     FUN = mean, 
  #                                     align = "right", 
  #                                     fill = NA, 
  #                                     na.rm = T)) %>%
  #   dplyr::summarise(max30DayYear = max(mm30Day, na.rm = T))
  # derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, mean30Day, by = c("featureid", "year"))
  # rm(mean30Day)
  
  # Number of days with stream temp > threshold
  derivedfeatureidMetrics <- calcThresholdDaysAnnual(byfeatureidYear, derivedfeatureidMetrics, 18) %>%
    dplyr::rename(days.18 = days) %>%
    dplyr::mutate(days.18 = as.numeric(days.18)) %>%
    dplyr::mutate(days.18 = as.numeric(ifelse(!is.na(maxTemp) & is.na(days.18), 0, days.18)))
  derivedfeatureidMetrics <- calcThresholdDaysAnnual(byfeatureidYear, derivedfeatureidMetrics, 20) %>%
    dplyr::rename(days.20 = days) %>%
    dplyr::mutate(days.20 = as.numeric(days.20)) %>%
    dplyr::mutate(days.20 = as.numeric(ifelse(!is.na(maxTemp) & is.na(days.20), 0, days.20)))
  derivedfeatureidMetrics <- calcThresholdDaysAnnual(byfeatureidYear, derivedfeatureidMetrics, 22) %>%
    dplyr::rename(days.22 = days) %>%
    dplyr::mutate(days.22 = as.numeric(days.22)) %>%
    dplyr::mutate(days.22 = as.numeric(ifelse(!is.na(maxTemp) & is.na(days.22), 0, days.22)))
  #if(class(threshold) == "numeric") derivedfeatureidMetrics <- calcThresholdDaysAnnual(byfeatureidYear, derivedfeatureidMetrics, threshold)
  
  # CT DEEP Thresholds
  #derivedfeatureidMetrics <- calcYearsCold(byfeatureidYear, derivedfeatureidMetrics, states = c("CT"))
  
  # Resistance to peak air temperature
  ## This probably makes the most sense during minimum flow periods but we don't have a sufficient flow model
  ## 60 or 90 days max air temp?
  # error if use standardized values rather than original scale
  # dOY.max.warm <- byfeatureidYear %>%
  #    mutate(warm.90 = rollsum(x = airTemp, 90, align = "right", fill = NA))
  #  dOY.max.warm <- dOY.max.warm %>%
  #   group_by(featureid, year, add = TRUE) %>%
  #  filter(warm.90 == max(warm.90, na.rm = T)) %>%
  # select(dOY)
  meanResist <- byfeatureidYear %>%
    #filter(dOY > dOY.max.warm$dOY - 90 & dOY <= dOY.max.warm$dOY) %>%
    dplyr::filter(dOY >= 152 & dOY < 244) %>% # clip to summer
    dplyr::mutate(absResid = abs(airTemp - tempPredicted)) %>%
    dplyr::summarise(resistance = sum(absResid))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanResist, by = c("featureid", "year"))
  rm(meanResist)
  
  # User broom and dplyr to get TS for each feature id but make sure it handles NA for entire featureid or entire datasets
  # Orange %>% group_by(Tree) %>% do(tidy(lm(age ~ circumference, data=.)))
  
  # Thermal Sensitivity ------------need to add year into it------------
  # byfeatureid_complete <- byfeatureidYear %>%
  #   dplyr::filter(!is.na(tempPredicted) & !is.na(airTemp))
  # 
  # group_count <- byfeatureid_complete %>%
  #   dplyr::summarise(count = n())
  # 
  # if(nrow(byfeatureid_complete) >= 5) { # actually want to check for 5 obs within any group
  #   byfeatureid_complete <- left_join(byfeatureid_complete, group_count) %>%
  #     dplyr::filter(count >= 5)
  #   if(nrow(byfeatureid_complete) >= 5 & length(unique(byfeatureid_complete$featureid)) > 1) {
  #     TS <- byfeatureid_complete %>% do(broom::tidy(lm(tempPredicted ~ airTemp, data = .))) %>%
  #       dplyr::filter(term == "airTemp") %>%
  #       dplyr::select(featureid, TS = estimate)
  #     derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, TS, by = c("featureid", "year"))
  #     rm(TS)
  #   } else {
  #     derivedfeatureidMetrics$TS <- NA_real_
  #   }
  # } else {
  #   derivedfeatureidMetrics$TS <- NA_real_
  # }
  
  
  # RMSE for each featureid (flag highest)
  bar <- byfeatureidYear %>%
    dplyr::filter(!(is.na(temp) & !is.na(tempPredicted))) %>%
    dplyr::mutate(error = temp - tempPredicted) %>%
    dplyr::select(featureid, year, temp, tempPredicted, error)
  
  if(dim(bar)[1] > 0) {
    error_metrics <- bar %>%
      dplyr::summarise(RMSE = rmse(error),
                       MAE = mae(error),
                       NSE = nse(temp, tempPredicted))
    derivedfeatureidMetrics <- dplyr::left_join(derivedfeatureidMetrics, error_metrics, by = c("featureid", "year"))
  } else {
    derivedfeatureidMetrics$RMSE <- NA_real_
    derivedfeatureidMetrics$MAE <- NA_real_
    derivedfeatureidMetrics$NSE <- NA_real_
    #     derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    #       dplyr::mutate(meanRMSE = NA) %>%
    #       dplyr::mutate(meanRMSE = as.numeric(meanRMSE))
  }
  rm(bar)
  
  ############# ADD RMSE by Trend in addition to predicted with AR1 #########
  
  # RMSE for each featureid (flag highest)
  bar <- byfeatureidYear %>%
    filter(!(is.na(temp) & !is.na(trend))) %>%
    mutate(error = temp - trend) %>%
    dplyr::select(featureid, year, temp, tempPredicted, error)
  
  if(dim(bar)[1] > 0) {
    error_metrics <- bar %>%
      dplyr::summarise(RMSE_trend = rmse(error),
                       MAE_trend = mae(error),
                       NSE_trend = nse(temp, tempPredicted))
    derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, error_metrics, by = c("featureid", "year"))
  } else {
    derivedfeatureidMetrics$RMSE_trend <- NA_real_
    derivedfeatureidMetrics$MAE_trend <- NA_real_
    derivedfeatureidMetrics$NSE_trend <- NA_real_
    #     derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    #       dplyr::mutate(meanRMSE = NA) %>%
    #       dplyr::mutate(meanRMSE = as.numeric(meanRMSE))
  }
  rm(bar)
  
  ########################
  
  #derived.site.metrics <- derivedfeatureidMetrics
  #rm(data)
  #rm(derivedfeatureidMetrics)
  #rm(foo)
  #gc()
  ################################
  
  # derived.site.metrics <- deriveMetrics(data)
  
  #derived_metrics
  #   if(exists("metrics") == FALSE) {
  #     metrics <- derivedfeatureidMetrics
  #     #print("no")
  #   } else {
  #     metrics <- rbind(metrics, derivedfeatureidMetrics)
  #   }
  
  #saveRDS(metrics.lat.lon, file = paste0(data_dir, "/derived_site_metrics.RData"))
  #write.table(metrics.lat.lon, file = paste0(data_dir, "/derived_site_metrics.csv"), sep = ',', row.names = F, append = TRUE)
  
  metrics <- as.data.frame(derivedfeatureidMetrics)
  
  return(metrics)
}


#' @title calcThresholdDaysAnnual
#'
#' @description
#' \code{calcThresholdDaysAnnual} Calculates the median number of days per year that a stream is predicted to exceed a threshold temperature
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param temp.threshold Optional numeric temperature threshold value in degrees C
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' blah, blah, blah, something, something
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calcThresholdDaysAnnual <- function(grouped.df, derived.df, temp.threshold, summer = FALSE) {
  var.name <- paste0("days.", temp.threshold)
  if(summer) {
    meanDays <- grouped.df %>%
      dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
      dplyr::filter(month >= 6 & month <= 8) %>%
      dplyr::filter(tempPredicted > temp.threshold) %>%
      dplyr::summarise(days = n()) %>%
      #dplyr::select(-month) %>%
      dplyr::mutate(days = ifelse(is.na(days), 0, days))
    #setNames(c("featureid", var.name))
  } else { 
    meanDays <- grouped.df %>%
      dplyr::filter(tempPredicted > temp.threshold) %>%
      dplyr::summarise(days = n())
    if(dim(meanDays)[1] > 0) {
      meanDays <- meanDays %>%
        dplyr::mutate(days = ifelse(is.na(days), 0, days))
      # setNames(c("featureid", var.name))
    } else {
      meanDays <- meanDays %>%
        dplyr::mutate(days = NA)
    }
  }
  derived.df <- dplyr::left_join(derived.df, dplyr::select(meanDays, featureid, days), by = "featureid")
  #derived.df[ , var.name][is.na(derived.df[ , var.name])] <- 0
  rm(meanDays)
  
  return(derived.df)
}



#' @title deriveMetrics
#'
#' @description
#' \code{deriveMetrics} Wrapper to calculate derived metrics in parallel using foreach %dopar%
#'
#' @param data Dataframe created with the prepData and predictTemp functions
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' blah, blah, blah, something, something
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
deriveMetrics <- function(data) {
  library(dplyr)
  
  byfeatureid <- group_by(data, featureid)
  byfeatureidYear <- group_by(byfeatureid, year, add = TRUE)
  #(maxTempfeatureid <- dplyr::dplyr::summarise(byfeatureid, max(tempPredicted, na.rm = T)))
  
  derivedfeatureidMetrics <- dplyr::select(byfeatureid, featureid) %>%
    distinct
  
  # Mean maximum daily mean temperature by featureid (over years) - this must be calculated before totalObs to work with NA and left join if all data missing
  meanMaxTemp <- byfeatureidYear %>%
    dplyr::summarise(maxTempPredicted = max(tempPredicted, na.rm = T)) %>%
    dplyr::summarise(meanMaxTemp = mean(maxTempPredicted))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanMaxTemp, by = "featureid")
  rm(meanMaxTemp)
  
  # total observations (days with data) per featureid
  totalObs <- byfeatureidYear %>%
    dplyr::filter(!is.na(temp)) %>%
    dplyr::summarise(Obs = n()) %>%
    dplyr::summarise(totalObs = sum(Obs)) %>%
    dplyr::select(featureid, totalObs)
  derivedfeatureidMetrics <- dplyr::left_join(derivedfeatureidMetrics, totalObs, by = "featureid") %>%
    dplyr::mutate(totalObs = as.numeric(ifelse(is.na(totalObs), 0, as.numeric(totalObs))))
  derivedfeatureidMetrics <- dplyr::select(derivedfeatureidMetrics, featureid, totalObs, meanMaxTemp)
  
  # Maximum max daily mean temperature
  maxMaxTemp <- byfeatureidYear %>%
    dplyr::summarise(maxTemp = max(tempPredicted, na.rm = T)) %>%
    dplyr::summarise(maxMaxTemp = max(maxTemp))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, maxMaxTemp, by = "featureid")
  rm(maxMaxTemp)
  
  # Mean July Air Temp
  meanJulyAir <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month == 7) %>%
    dplyr::summarise(JulyTemp = mean(airTemp, na.rm = T)) %>%
    dplyr::summarise(meanJulyAir = mean(JulyTemp))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanJulyAir, by = "featureid")
  rm(meanJulyAir)
  
  # Mean July Summer Temp
  meanJulyTemp <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month == 7) %>%
    dplyr::summarise(JulyTemp = mean(tempPredicted, na.rm = T)) %>%
    dplyr::summarise(meanJulyTemp = mean(JulyTemp))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanJulyTemp, by = "featureid")
  rm(meanJulyTemp)
  
  # Mean July Obs Temp
  meanJulyObs <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month == 7)
  
  nJulyObs <- meanJulyObs %>%
    dplyr::filter(!is.na(temp)) %>%
    dplyr::summarise(nJulyObs = n())
  
  nJulyObsYears <- nJulyObs %>%
    dplyr::filter(nJulyObs >= 25) %>% # at least 25 days with observed data
    dplyr::summarise(nYearsJulyObs = n())
  
  meanJulyObs <- nJulyObs %>%
    left_join(meanJulyObs) %>%
    dplyr::filter(nJulyObs >= 25) %>%
    dplyr::summarise(meanJulyObs = mean(temp, na.rm = T)) 
  
  derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    dplyr::left_join(meanJulyObs) %>%
    dplyr::left_join(nJulyObsYears)
  rm(meanJulyObs)
  
  # Mean Aug Temp
  meanAugTemp <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month == 8) %>%
    dplyr::summarise(AugTemp = mean(tempPredicted, na.rm = T)) %>%
    dplyr::summarise(meanAugTemp = mean(AugTemp))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanAugTemp, by = "featureid")
  rm(meanAugTemp)
  
  # Mean Summer Temp
  meanSummerTemp <- byfeatureidYear %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    dplyr::filter(month >= 6 & month <= 8) %>%
    dplyr::summarise(SummerTemp = mean(tempPredicted, na.rm = T)) %>%
    dplyr::summarise(meanSummerTemp = mean(SummerTemp))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanSummerTemp, by = "featureid")
  rm(meanSummerTemp)
  
  # Annual mean 30-day maximum mean daily temperature
  mean30Day <- byfeatureidYear %>%
    dplyr::arrange(featureid, year, dOY) %>%
    dplyr::mutate(mm30Day = rollapply(data = tempPredicted, 
                                      width = 30, 
                                      FUN = mean, 
                                      align = "right", 
                                      fill = NA, 
                                      na.rm = T)) %>%
    dplyr::summarise(max30DayYear = max(mm30Day, na.rm = T)) %>%
    dplyr::summarise(mean30DayMax = mean(max30DayYear, na.rm = T))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, mean30Day, by = "featureid")
  rm(mean30Day)
  
  # Number of days with stream temp > threshold
  derivedfeatureidMetrics <- calcThresholdDays(byfeatureidYear, derivedfeatureidMetrics, 18) %>%
    dplyr::rename(meanDays.18 = meanDays) %>%
    dplyr::mutate(meanDays.18 = as.numeric(meanDays.18)) %>%
    dplyr::mutate(meanDays.18 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(meanDays.18), 0, meanDays.18)))
  derivedfeatureidMetrics <- calcThresholdDays(byfeatureidYear, derivedfeatureidMetrics, 20) %>%
    dplyr::rename(meanDays.20 = meanDays) %>%
    dplyr::mutate(meanDays.20 = as.numeric(meanDays.20)) %>%
    dplyr::mutate(meanDays.20 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(meanDays.20), 0, meanDays.20)))
  derivedfeatureidMetrics <- calcThresholdDays(byfeatureidYear, derivedfeatureidMetrics, 22) %>%
    dplyr::rename(meanDays.22 = meanDays) %>%
    dplyr::mutate(meanDays.22 = as.numeric(meanDays.22)) %>%
    dplyr::mutate(meanDays.22 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(meanDays.22), 0, meanDays.22)))
  #if(class(threshold) == "numeric") derivedfeatureidMetrics <- calcThresholdDays(byfeatureidYear, derivedfeatureidMetrics, threshold)
  
  # CT DEEP Thresholds
  #derivedfeatureidMetrics <- calcYearsCold(byfeatureidYear, derivedfeatureidMetrics, states = c("CT"))
  # Number and frequency of years with mean max over threshold
  derivedfeatureidMetrics <- calcYearsMaxTemp(grouped.df = byfeatureidYear, derived.df = derivedfeatureidMetrics, temp.threshold = 18) %>%
    dplyr::mutate(yearsMaxTemp = as.numeric(ifelse(is.na(yearsMaxTemp), 0, as.numeric(yearsMaxTemp)))) %>%
    dplyr::rename(yearsMaxTemp.18 = yearsMaxTemp) %>%
    dplyr::mutate(freqMaxTemp.18 = yearsMaxTemp.18 / length(unique(byfeatureidYear$year)))
  derivedfeatureidMetrics[which(is.na(derivedfeatureidMetrics$meanMaxTemp)), "yearsMaxTemp.18"] <- NA
  derivedfeatureidMetrics[which(is.na(derivedfeatureidMetrics$meanMaxTemp)), "freqMaxTemp.18"] <- NA
  derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    dplyr::mutate(yearsMaxTemp.18 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(yearsMaxTemp.18), 0, yearsMaxTemp.18)),
                  yearsMaxTemp.18 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(yearsMaxTemp.18), 0, yearsMaxTemp.18)))
  
  derivedfeatureidMetrics <- calcYearsMaxTemp(grouped.df = byfeatureidYear, derived.df = derivedfeatureidMetrics, temp.threshold = 20) %>%
    dplyr::mutate(yearsMaxTemp = as.numeric(ifelse(is.na(yearsMaxTemp), 0, as.numeric(yearsMaxTemp)))) %>%
    dplyr::rename(yearsMaxTemp.20 = yearsMaxTemp) %>%
    dplyr::mutate(freqMaxTemp.20 = yearsMaxTemp.20 / length(unique(byfeatureidYear$year)))
  derivedfeatureidMetrics[which(is.na(derivedfeatureidMetrics$meanMaxTemp)), "yearsMaxTemp.20"] <- NA
  derivedfeatureidMetrics[which(is.na(derivedfeatureidMetrics$meanMaxTemp)), "freqMaxTemp.20"] <- NA
  derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    dplyr::mutate(yearsMaxTemp.20 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(yearsMaxTemp.20), 0, yearsMaxTemp.20)),
                  freqMaxTemp.20 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(freqMaxTemp.20), 0, freqMaxTemp.20)))
  
  derivedfeatureidMetrics <- calcYearsMaxTemp(grouped.df = byfeatureidYear, derived.df = derivedfeatureidMetrics, temp.threshold = 22) %>%
    dplyr::mutate(yearsMaxTemp = as.numeric(ifelse(is.na(yearsMaxTemp), 0, as.numeric(yearsMaxTemp)))) %>%
    dplyr::rename(yearsMaxTemp.22 = yearsMaxTemp) %>%
    dplyr::mutate(freqMaxTemp.22 = yearsMaxTemp.22 / length(unique(byfeatureidYear$year)))
  
  derivedfeatureidMetrics[which(is.na(derivedfeatureidMetrics$meanMaxTemp)), "yearsMaxTemp.22"] <- NA
  derivedfeatureidMetrics[which(is.na(derivedfeatureidMetrics$meanMaxTemp)), "freqMaxTemp.22"] <- NA
  derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    dplyr::mutate(yearsMaxTemp.22 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(yearsMaxTemp.22), 0, yearsMaxTemp.22)),
                  freqMaxTemp.22 = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(freqMaxTemp.22), 0, freqMaxTemp.22)))
  
  derivedfeatureidMetrics <- calcYearsMaxTemp(grouped.df = byfeatureidYear, derived.df = derivedfeatureidMetrics, temp.threshold = 23.5) %>%
    dplyr::mutate(yearsMaxTemp = as.numeric(ifelse(is.na(yearsMaxTemp), 0, as.numeric(yearsMaxTemp)))) %>%
    dplyr::rename(yearsMaxTemp.23.5 = yearsMaxTemp) %>%
    dplyr::mutate(freqMaxTemp.23.5 = yearsMaxTemp.23.5 / length(unique(byfeatureidYear$year)))
  derivedfeatureidMetrics[which(is.na(derivedfeatureidMetrics$meanMaxTemp)), "yearsMaxTemp.23.5"] <- NA
  derivedfeatureidMetrics[which(is.na(derivedfeatureidMetrics$meanMaxTemp)), "freqMaxTemp.23.5"] <- NA
  derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    dplyr::mutate(yearsAcute.MADEP = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(yearsMaxTemp.23.5), 0, yearsMaxTemp.23.5)),
                  freqAcute.MADEP = as.numeric(ifelse(!is.na(meanMaxTemp) & is.na(freqMaxTemp.23.5), 0, freqMaxTemp.23.5)))
  
  # Resistance to peak air temperature
  ## This probably makes the most sense during minimum flow periods but we don't have a sufficient flow model
  ## 60 or 90 days max air temp?
  # error if use standardized values rather than original scale
  # dOY.max.warm <- byfeatureidYear %>%
  #    mutate(warm.90 = rollsum(x = airTemp, 90, align = "right", fill = NA))
  #  dOY.max.warm <- dOY.max.warm %>%
  #   group_by(featureid, year, add = TRUE) %>%
  #  filter(warm.90 == max(warm.90, na.rm = T)) %>%
  # select(dOY)
  meanResist <- byfeatureidYear %>%
    #filter(dOY > dOY.max.warm$dOY - 90 & dOY <= dOY.max.warm$dOY) %>%
    dplyr::filter(dOY >= 152 & dOY < 244) %>% # clip to summer
    dplyr::mutate(absResid = abs(airTemp - tempPredicted)) %>%
    dplyr::summarise(resistance = sum(absResid)) %>%
    dplyr::summarise(meanResist = mean(resistance))
  derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, meanResist, by = "featureid")
  rm(meanResist)
  
  # User broom and dplyr to get TS for each feature id but make sure it handles NA for entire featureid or entire datasets
  # Orange %>% group_by(Tree) %>% do(tidy(lm(age ~ circumference, data=.)))
  
  # Thermal Sensitivity
  byfeatureid_complete <- byfeatureid %>%
    dplyr::filter(!is.na(tempPredicted) & !is.na(airTemp))
  
  group_count <- byfeatureid_complete %>%
    dplyr::summarise(count = n())
  
  if(nrow(byfeatureid_complete) >= 5) { # actually want to check for 5 obs within any group
    byfeatureid_complete <- left_join(byfeatureid_complete, group_count) %>%
      dplyr::filter(count >= 5)
    if(nrow(byfeatureid_complete) >= 5 & length(unique(byfeatureid_complete$featureid)) > 1) {
      TS <- byfeatureid_complete %>% do(broom::tidy(lm(tempPredicted ~ airTemp, data = .))) %>%
        dplyr::filter(term == "airTemp") %>%
        dplyr::select(featureid, TS = estimate)
      derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, TS, by = "featureid")
      rm(TS)
    } else {
      derivedfeatureidMetrics$TS <- NA_real_
    }
  } else {
    derivedfeatureidMetrics$TS <- NA_real_
  }
  
  
  # RMSE for each featureid (flag highest)
  bar <- byfeatureidYear %>%
    dplyr::filter(!(is.na(temp) & !is.na(tempPredicted))) %>%
    dplyr::mutate(error = temp - tempPredicted)
  
  if(dim(bar)[1] > 0) {
    error_metrics <- bar %>%
      dplyr::summarise(RMSE = rmse(error),
                       MAE = mae(error),
                       NSE = nse(temp, tempPredicted)) %>%
      dplyr::summarise(meanRMSE = mean(RMSE, na.rm = T),
                       meanMAE = mean(MAE, na.rm = T),
                       meanNSE = mean(NSE, na.rm = T))
    derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, dplyr::select(error_metrics, featureid, meanRMSE, meanMAE, meanNSE), by = "featureid")
  } else {
    derivedfeatureidMetrics$meanRMSE <- NA_real_
    derivedfeatureidMetrics$meanMAE <- NA_real_
    derivedfeatureidMetrics$meanNSE <- NA_real_
    #     derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    #       dplyr::mutate(meanRMSE = NA) %>%
    #       dplyr::mutate(meanRMSE = as.numeric(meanRMSE))
  }
  rm(bar)
  
  ############# ADD RMSE by Trend in addition to predicted with AR1 #########
  
  # RMSE for each featureid (flag highest)
  bar <- byfeatureidYear %>%
    filter(!(is.na(temp) & !is.na(trend))) %>%
    mutate(error = temp - trend)
  
  if(dim(bar)[1] > 0) {
    error_metrics <- bar %>%
      dplyr::summarise(RMSE_trend = rmse(error),
                       MAE_trend = mae(error),
                       NSE_trend = nse(temp, tempPredicted)) %>%
      dplyr::summarise(meanRMSE_trend = mean(RMSE_trend, na.rm = T),
                       meanMAE_trend = mean(MAE_trend, na.rm = T),
                       meanNSE_trend = mean(NSE_trend, na.rm = T))
    derivedfeatureidMetrics <- left_join(derivedfeatureidMetrics, dplyr::select(error_metrics, featureid, meanRMSE_trend, meanMAE_trend, meanNSE_trend), by = "featureid")
  } else {
    derivedfeatureidMetrics$meanRMSE_trend <- NA_real_
    derivedfeatureidMetrics$meanMAE_trend <- NA_real_
    derivedfeatureidMetrics$meanNSE_trend <- NA_real_
    #     derivedfeatureidMetrics <- derivedfeatureidMetrics %>%
    #       dplyr::mutate(meanRMSE = NA) %>%
    #       dplyr::mutate(meanRMSE = as.numeric(meanRMSE))
  }
  rm(bar)
  
  ########################
  
  #derived.site.metrics <- derivedfeatureidMetrics
  #rm(data)
  #rm(derivedfeatureidMetrics)
  #rm(foo)
  #gc()
  ################################
  
  # derived.site.metrics <- deriveMetrics(data)
  
  #derived_metrics
  #   if(exists("metrics") == FALSE) {
  #     metrics <- derivedfeatureidMetrics
  #     #print("no")
  #   } else {
  #     metrics <- rbind(metrics, derivedfeatureidMetrics)
  #   }
  
  #saveRDS(metrics.lat.lon, file = paste0(data_dir, "/derived_site_metrics.RData"))
  #write.table(metrics.lat.lon, file = paste0(data_dir, "/derived_site_metrics.csv"), sep = ',', row.names = F, append = TRUE)
  
  metrics <- as.data.frame(derivedfeatureidMetrics)
  
  return(metrics)
}


#' @title derive_metrics_par
#'
#' @description
#' \code{derive_metrics_par} Wrapper to pull data from database, prepare data, and calculate derived metrics in parallel using foreach %dopar%
#'
#' @param i Iteration passed from foreach loop
#' @param chunk.size Number of catchments to process at 1 time, default = 50 based on testing speed with pulling from the database and memory issues
#' @param catchmentid
#' @param springFallBPs
#' @param tempDataSync
#' @param df_covariates_upstream
#' @param featureid_huc8
#' @param featureid_site=featureid_site
#' @param featureid_lat_lon
#' @param coef.list
#' @param cov.list
#' @param var.names
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' blah, blah, blah, something, something
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
derive_metrics_par <- function(i, chunk.size = 50, catchmentid, springFallBPs, tempDataSync, df_covariates_upstream, featureid_huc8, featureid_site, featureid_lat_lon, coef.list, cov.list, var.names, con) {
  
  n.catches <- length(catchmentid)
  k <- i*chunk.size
  if(k <= n.catches) {
    catches <- catchmentid[(1+(i-1)*chunk.size):k]
  } else {
    catches <- catchmentid[(1+(i-1)*chunk.size):n.catches]
  }
  catches_string <- paste(catches, collapse = ', ')
  
  # reconnect to database if lost
  #   if(isPostgresqlIdCurrent(con) == FALSE) {
  #     drv <- dbDriver("PostgreSQL")
  #     con <- dbConnect(drv, dbname='sheds', host='ecosheds.org', user=options('SHEDS_USERNAME'), password=options('SHEDS_PASSWORD'))
  #   }
  
  # pull daymet records
  qry_daymet <- paste0("SELECT featureid, date, tmax, tmin, prcp, dayl, srad, vp, swe, (tmax + tmin) / 2.0 AS airTemp FROM daymet WHERE featureid IN (", catches_string, ") ;")
  
  rs <- dbSendQuery(con, statement = qry_daymet)
  climateData <- fetch(rs, n=-1)
  
  
  mean.spring.bp <- mean(dplyr::filter(springFallBPs, finalSpringBP != "Inf")$finalSpringBP, na.rm = T)
  mean.fall.bp <- mean(dplyr::filter(springFallBPs, finalFallBP != "Inf")$finalFallBP, na.rm = T)
  
  foo <- springFallBPs %>%
    dplyr::mutate(site = as.character(site),
                  featureid = as.integer(site),
                  finalSpringBP = ifelse(finalSpringBP == "Inf" | is.na(finalSpringBP), mean.spring.bp, finalSpringBP),
                  finalFallBP = ifelse(finalFallBP == "Inf" | is.na(finalFallBP), mean.fall.bp, finalFallBP))
  
  fullDataSync <- climateData %>%
    left_join(df_covariates_upstream, by=c('featureid')) %>%
    left_join(dplyr::select(tempDataSync, date, featureid, site, temp), by = c('date', 'featureid')) %>%
    left_join(featureid_huc8, by = c('featureid')) %>%
    left_join(featureid_lat_lon, by = c('featureid')) %>%
    dplyr::mutate(year = as.numeric(format(date, "%Y")),
                  airTemp = (tmax + tmin)/2) %>%
    left_join(dplyr::select(foo, -site), by = c('featureid', 'year')) %>%
    dplyr::mutate(dOY = yday(date)) %>%
    #dplyr::mutate(huc = huc8) %>%
    dplyr::filter(dOY >= finalSpringBP & dOY <= finalFallBP | is.na(finalSpringBP) | is.na(finalFallBP & finalSpringBP != "Inf" & finalFallBP != "Inf")) %>%
    dplyr::filter(dOY >= mean.spring.bp & dOY <= mean.fall.bp) %>%
    dplyr::filter(AreaSqKM < 400)
  
  # Order by group and date
  fullDataSync <- fullDataSync[order(fullDataSync$featureid, fullDataSync$year, fullDataSync$dOY),]
  
  # For checking the order of fullDataSync
  fullDataSync$count <- 1:length(fullDataSync$year)
  
  fullDataSync <- fullDataSync[order(fullDataSync$count),] # just to make sure fullDataSync is ordered for the slide function
  
  # moving means instead of lagged terms in the future
  fullDataSync <- fullDataSync %>%
    group_by(featureid, year) %>%
    arrange(featureid, year, dOY) %>%
    mutate(airTempLagged1 = lag(airTemp, n = 1, fill = NA),
           temp5p = rollapply(data = airTempLagged1, 
                              width = 5, 
                              FUN = mean, 
                              align = "right", 
                              fill = NA, 
                              na.rm = T),
           temp7p = rollapply(data = airTempLagged1, 
                              width = 7, 
                              FUN = mean, 
                              align = "right", 
                              fill = NA, 
                              na.rm = T),
           prcp2 = rollsum(x = prcp, 2, align = "right", fill = NA),
           prcp7 = rollsum(x = prcp, 7, align = "right", fill = NA),
           prcp30 = rollsum(x = prcp, 30, align = "right", fill = NA))
  
  var.names <- c("airTemp", 
                 "temp7p",
                 "prcp", 
                 "prcp2",
                 "prcp7",
                 "prcp30",
                 "dOY", 
                 "forest", 
                 "herbaceous", 
                 "agriculture", 
                 "devel_hi", 
                 "developed",
                 "AreaSqKM",  
                 "allonnet",
                 "alloffnet",
                 "surfcoarse", 
                 "srad", 
                 "dayl", 
                 "swe")
  
  fullDataSync <- fullDataSync %>%
    mutate(HUC8 = as.character(HUC8),
           huc8 = as.character(HUC8),
           site = as.numeric(as.factor(featureid))) 
  
  fullDataSyncS <- stdCovs(x = fullDataSync, y = tempDataSync, var.names = var.names)
  
  fullDataSyncS <- addInteractions(fullDataSyncS)
  
  fullDataSyncS <- indexDeployments(fullDataSyncS, regional = TRUE)
  
  fullDataSyncS <- predictTemp(data = fullDataSyncS, coef.list = coef.list, cov.list = cov.list, featureid_site = featureid_site)
  
  fullDataSync <- left_join(fullDataSync, select(fullDataSyncS, featureid, date, tempPredicted), by = c("featureid", "date"))
  
  mean.pred <- mean(fullDataSync$tempPredicted, na.rm = T)
  
  if(mean.pred == "NaN") {
    cat(paste0(i, " of ", n.loops, " loops has no predicted temperatures"))
  } 
  derived.site.metrics <- deriveMetrics(fullDataSync)
  
  return(derived.site.metrics)
}




#' @title calcThresholdDays
#'
#' @description
#' \code{calcThresholdDays} Calculates the median number of days per year that a stream is predicted to exceed a threshold temperature
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param temp.threshold Optional numeric temperature threshold value in degrees C
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' blah, blah, blah, something, something
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calcThresholdDays <- function(grouped.df, derived.df, temp.threshold, summer = FALSE) {
  var.name <- paste0("meanDays.", temp.threshold)
  if(summer) {
    meanDays <- grouped.df %>%
      dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
      dplyr::filter(month >= 6 & month <= 8) %>%
      dplyr::filter(tempPredicted > temp.threshold) %>%
      dplyr::summarise(days = n()) %>%
      dplyr::summarise(meanDays = round(median(days, na.rm = T))) %>%
      dplyr::select(-month) %>%
      dplyr::mutate(meanDays = ifelse(is.na(meanDays), 0, meanDays))
    #setNames(c("featureid", var.name))
  } else { 
    meanDays <- grouped.df %>%
      dplyr::filter(tempPredicted > temp.threshold) %>%
      dplyr::summarise(days = n())
    if(dim(meanDays)[1] > 0) {
      meanDays <- meanDays %>%
        dplyr::summarise(meanDays = round(median(days, na.rm = T))) %>%
        dplyr::mutate(meanDays = ifelse(is.na(meanDays), 0, meanDays))
      # setNames(c("featureid", var.name))
    } else {
      meanDays <- meanDays %>%
        dplyr::mutate(meanDays = NA)
    }
  }
  derived.df <- dplyr::left_join(derived.df, select(meanDays, featureid, meanDays), by = "featureid")
  #derived.df[ , var.name][is.na(derived.df[ , var.name])] <- 0
  rm(meanDays)
  
  return(derived.df)
}

#' @title calcYearsMaxTemp
#'
#' @description
#' \code{calcYearsMaxTemp} Calculates the number of years and frequency of years that a stream ever (even for 1 day) exceeds a threshold temperature
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param temp.threshold Optional numeric temperature threshold value in degrees C
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' blah, blah, blah, something, something
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calcYearsMaxTemp <- function(grouped.df, derived.df, temp.threshold, summer = FALSE) {
  library(lazyeval)
  var.name1 <- paste0("yearsMaxTemp.", temp.threshold)
  var.name2 <- paste0("freqMax.", temp.threshold)
  dots <- list(~n)
  if(summer) {
    yearsMaxTemp <- grouped.df %>%
      dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
      dplyr::filter(month >= 6 & month <= 8) %>%
      dplyr::summarise(maxTemp = max(tempPredicted, na.rm = T)) %>%
      dplyr::filter(maxTemp > temp.threshold) %>%
      dplyr::summarise(yearsMaxTemp = n()) %>%
      dplyr::mutate(yearsMaxTemp = ifelse(is.na(yearsMaxTemp), 0, yearsMaxTemp)) %>%
      dplyr::mutate(freqMax = yearsMaxTemp / length(unique(grouped.df$year))) %>%
      dplyr::mutate(freqMax = ifelse(is.na(freqMax), 0, freqMax)) %>%
      dplyr::select(-month)
    #setNames(c("featureid", var.name1, var.name2))
  } else {
    yearsMaxTemp <- grouped.df %>%
      dplyr::summarise(maxTemp = max(tempPredicted, na.rm = T)) %>%
      dplyr::filter(maxTemp > temp.threshold) %>%
      dplyr::summarise(yearsMaxTemp = n()) 
    if(dim(yearsMaxTemp)[1] > 0) {
      yearsMaxTemp <- yearsMaxTemp %>%
        dplyr::mutate(yearsMaxTemp = ifelse(is.na(yearsMaxTemp), 0, yearsMaxTemp)) #%>%
      #setNames(c("featureid", var.name1, var.name2))
    } else {
      yearsMaxTemp <- yearsMaxTemp %>%
        dplyr::mutate(yearsMaxTemp = NA) #%>%
      #setNames(c("featureid", var.name1, var.name2))
    }
  }
  derived.df <- dplyr::left_join(derived.df, yearsMaxTemp, by = "featureid")
  #derived.df[ , var.name1][is.na(derived.df[ , var.name1])] <- 0
  #derived.df[ , var.name2][is.na(derived.df[ , var.name2])] <- 0
  rm(yearsMaxTemp)
  
  return(derived.df)
}


#' @title calcConsecExceed
#'
#' @description
#' \code{calcConsecExceed} Calculate the mean and max (and quantiles?) of consecutive days per year where temperatures exceed some temperature threshold
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param threshold  Temperature threshold above which causes physiological problems for species of interest
#' @param summer Logical calculate metric just for June-August. Default FALSE
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' Maine is interested in adding another set of derived metrics that represented runs of consecutive days where the stream temperature exceeded some threshold. Maybe something like mean and max number of consecutive days per year where temp > 22 deg C, since that's closely tied to fish stress. 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calcConsecExceed <- function(y, threshold = 22) {
  if(all(is.na(y$tempPredicted))) {
    warning("all values are NA")
    maxConsec <- NA
    meanConsec <- NA
  } else if(all(y$tempPredicted <= threshold)) {
    maxConsec <- 0
    meanConsec <- 0
  } else {
    exceed <- y$tempPredicted > threshold
    consec <- rle(exceed)
    consecExceed <- consec$lengths[consec$values]
    maxConsec = max(consecExceed, na.rm = T)
    meanConsec = mean(consecExceed, na.rm = T)
  }
  d <- data.frame(maxConsec, meanConsec)
  names(d) <- paste(names(d), '_', threshold, sep='')
  return(d)
}

# calcConsecExceed <- function(grouped.df, derived,df, threshold, summer = FALSE) {
#   
#   if(summer) {
#     consecExceed <- grouped.df %>%
#       dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
#       dplyr::filter(month >= 6 & month <= 8) %>%
#       dplyr::mutate(exceed = ifelse(tempPredicted > temp.threshold, 1, 0),
#                     lag_exceed = lag(exceed, n = 1),
#                     #sum_exceed = exceed + lag_exceed,
#                     event_start = ifelse(exceed == 1 & lag_exceed != 1, 1, 0)) %>%#,
#                     #event_day = if(lag_exceed > 0, event_day[i-1]+1, exceed)
#                    # event = if(exceed == 1 & lag_exceed == 1, "name", NA)  # don't know how to name this without a for loop through 1.8 billion records
#       # dplyr::group_by(event)
#       
#       # easy to calculate the number of events but not the duration without resorting to a for loop or using LOTS of lag columns
#       dplyr::summarise(nExceedEvents = sum(event_start)) %>%
#       dplyr::summarise(nExceedEventsMean = mean(nExceed))
#   }
#   
#   derived.df <- dplyr::left_join(derived.df, consecExceed, by = "featureid")
#   
#   rm(consecExceed)
#   
#   return(derived.df)
# }


#' @title calc7DADM
#'
#' @description
#' \code{calc7DADM} Calculates the 7-day average of the daily (max) temperatures in any 7-contiguous day period from June - September
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param threshold  Temperature threshold above which causes physiological problems for species of interest
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' This is part of MassDEP CALM criteria for evaluating Cold Water Streams for 305b reporting. MassDEP uses 20C as their threshold ("endpoint").
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calc7DADM <- function(grouped.df, derived.df, threshold) {

  pct7DADM <- grouped.df %>%
    dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
    
    dplyr::arrange(featureid, year, date) %>%
    dplyr::mutate(dadm = rollmean(x = tempPredicted, 7, align = "right", fill = NA))
    dplyr::filter(month >= 6 & month <= 9) %>% # filter has to come last so don't have NA at the beginning throwing off % of days
    dplyr::ungroup() %>%
    dplyr::group(featureid) %>%
    dplyr::mutate(overThreshold = ifelse(dadm > threshold, 1, 0)) %>%
    dplyr::summarise(pct7DADM = sum(overThreshold)/n())
  
  derived.df <- dplyr::left_join(derived.df, pct7DADM, by = "featureid")
  
  rm(pct7DADM)
  
}




#' @title calcYearsCold
#'
#' @description
#' \code{calcYearsCold} Calculate derived metrics from predicted stream temperatures
#'
#' @param grouped.df Dataframe grouped by featureid then year
#' @param derived.df Dataframe of derived metrics
#' @param state Character string for specific State 
#' 
#' @return Returns Dataframe of derived metrics (e.g. max temperature predicted over daymet record) for each featureid across all years
#' @details
#' Different states have different classifications of cold, cool, and warm streams for regulatory purposes. This function takes the state abreviation and returns the number and frequency of years that any stream reach is expected to fall within each category.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
calcYearsCold <- function(grouped.df, derived.df, temp.threshold, states) {
  library(lazyeval)
  if(class(states) != "character") stop("calcYearsCold function expect state to be a character string of 2 letter abreviations")
  avail.states <- c("CT")
  if(!all(states %in% avail.states)) stop(paste("In calcYearsCold function, water temperature classifications only available for ", avail.states))
  
  if(any(states == "CT")) {
    # Cold
    var.name1 <- paste0("years.cold.", CT)
    var.name2 <- paste0("freq.cold.", CT)
    dots <- list(~n)
    yearsCold <- grouped.df %>%
      dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
      dplyr::filter(month >= 6 & month <= 8) %>%
      dplyr::summarise(meanTemp = mean(tempPredicted, na.rm = T)) %>%
      dplyr::filter(meanTemp < 18.29) %>%
      dplyr::summarise(yearsCold = n()) %>%
      dplyr::mutate(freqCold = yearsCold / length(unique(grouped.df$year))) %>%
      dplyr::mutate(yearsCold = ifelse(is.na(yearsCold), 0, yearsCold)) %>%
      dplyr::mutate(freqCold = ifelse(is.na(freqCold), 0, freqCold)) %>%
      dplyr::select(-month)
    setNames(c("featureid", var.name1, var.name2))
    derived.df <- merge(derived.df, yearsCold, by = "featureid")
    #derived.df[ , var.name1][is.na(derived.df[ , var.name1])] <- 0
    #derived.df[ , var.name2][is.na(derived.df[ , var.name2])] <- 0
    rm(yearsCold)
    # Cool
    var.name1 <- paste0("years.cool.", CT)
    var.name2 <- paste0("freq.cool.", CT)
    dots <- list(~n)
    yearsCool <- grouped.df %>%
      dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
      dplyr::filter(month >= 6 & month <= 8) %>%
      dplyr::summarise(meanTemp = mean(tempPredicted, na.rm = T)) %>%
      dplyr::filter(meanTemp >= 18.29 & meanTemp <= 21.70) %>%
      dplyr::summarise(yearsCool = n()) %>%
      dplyr::mutate(freqCold = yearsCool / length(unique(grouped.df$year))) %>%
      dplyr::mutate(yearsCold = ifelse(is.na(yearsCold), 0, yearsCold)) %>%
      dplyr::mutate(freqCold = ifelse(is.na(freqCold), 0, freqCold)) %>%
      dplyr::select(-month)
    setNames(c("featureid", var.name1, var.name2))
    derived.df <- left_join(derived.df, yearsCool, by = "featureid")
    #derived.df[ , var.name1][is.na(derived.df[ , var.name1])] <- 0
    #derived.df[ , var.name2][is.na(derived.df[ , var.name2])] <- 0
    rm(yearsCool)
    # Warm
    var.name1 <- paste0("years.warm.", CT)
    var.name2 <- paste0("freq.warm.", CT)
    dots <- list(~n)
    yearsCool <- grouped.df %>%
      dplyr::mutate(month = as.numeric(format(date, "%m"))) %>%
      dplyr::filter(month >= 6 & month <= 8) %>%
      dplyr::summarise(meanTemp = mean(tempPredicted, na.rm = T)) %>%
      dplyr::filter(meanTemp > 21.70) %>%
      dplyr::summarise(yearsCool = n()) %>%
      dplyr::mutate(freqCold = yearsCool / length(unique(grouped.df$year))) %>%
      dplyr::mutate(yearsCold = ifelse(is.na(yearsCold), 0, yearsCold)) %>%
      dplyr::mutate(freqCold = ifelse(is.na(freqCold), 0, freqCold)) %>%
      dplyr::select(-month)
    setNames(c("featureid", var.name1, var.name2))
    derived.df <- left_join(derived.df, yearsCool, by = "featureid")
    #derived.df[ , var.name1][is.na(derived.df[ , var.name1])] <- 0
    #derived.df[ , var.name2][is.na(derived.df[ , var.name2])] <- 0
    rm(yearsCool)
  }
  
  return(derived.df)
}