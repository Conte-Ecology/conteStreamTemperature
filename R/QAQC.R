#' @title obs_freq
#'
#' @description
#' \code{obs_freq} Median observations per day for each time series
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @return Median observations per day for each time series
#' @details
#' blah, blah, blah, something, something
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
obs_freq <- function(data) {
  obs_per_day <- data %>%
    dplyr::group_by(series_id, date) %>%
    dplyr::summarise(obs_per_day = n())
  
  median_obs <- obs_per_day %>%
    #dplyr::mutate(series_id_alt = paste0(series_id, "-", "a")) %>%
    dplyr::group_by(series_id) %>%
    #dplyr::mutate(obs_per_day = ifelse(is.na(obs_per_day), NA_integer_, obs_per_day)) %>%
    dplyr::summarise(median_freq = median(obs_per_day, na.rm = T)*1.0, min_n90 = median_freq*0.9) 
    #dplyr::summarise(median_freq = mean(obs_per_day, na.rm = T), min_n90 = median_freq*0.9) # make median when dplyr fixed
  
  data <- data %>%
    left_join(obs_per_day, by = c("series_id", "date")) %>%
    left_join(median_obs, by = c("series_id"))
  
  return(data)
}


#' @title flag_interval
#'
#' @description
#' \code{flag_interval} Uniformity of records within a day
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @return flag inconsistent time intervals within a day at a given location or series (logger deployment)
#' @details
#' check to make sure datetime was done correctly when guessing 12-hr vs 24-hr formats
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_interval <- function(data) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- ungroup(data)
  data <- group_by(data, series_id, date)
  data$row <- 1:nrow(data)
  data$time_prev <- c(NA_real_, data[2:nrow(data)-1, ]$datetime)
  data$series_prev <- c(NA_character_, data[1:nrow(data)-1, ]$series_id)
  data_series <- data %>%
    group_by(series_id, date) %>%
    mutate(series_start = ifelse(is.na(series_prev), 1, ifelse(series_id == series_prev, NA_real_, 1)),
           time_prev = ifelse(is.na(series_start), time_prev, NA_real_),
           d_time = ifelse(median_freq > 1, datetime - time_prev, NA_real_)) %>%
    summarise(flag_interval = ifelse(max(d_time, na.rm = T) != min(d_time, na.rm = T), TRUE, FALSE))
  
  data <- left_join(data, data_series, by = c("series_id", "date"))
  
  return(data)
}


#' @title flag_incomplete
#'
#' @description
#' \code{flag_incomplete} Flag incomplete days
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @return Flag days with less than 90% of median number of observations for that series
#' @details
#' blah, blah, blah, something, something
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_incomplete <- function(data) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- data %>%
    dplyr::mutate(flag_incomplete = ifelse(obs_per_day <= min_n90, TRUE, FALSE))
  return(data)
}


#' @title flag_hourly_rise
#'
#' @description
#' \code{flag_hourly_rise} Check for out of water (out of water vs. dry stream?)
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @param deg Threshold change in degrees for flag
#' @return Flag days with Rate of change in temperature per hour above deg
#' @details
#' Rate of change in temperature per hour
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_hourly_rise <- function(data, deg = 5) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- ungroup(data)
  data$row <- 1:nrow(data)
  data$temp_prev <- c(NA_real_, data[2:nrow(data)-1, ]$temp)
  data$series_prev <- c(NA_character_, data[1:nrow(data)-1, ]$series_id)
  data <- data %>%
    #group_by(series_id) %>%
    mutate(series_start = ifelse(is.na(series_prev), 1, ifelse(series_id == series_prev, NA_real_, 1)),
           temp_prev = ifelse(is.na(series_start), temp_prev, NA_real_),
          # d_temp = ifelse(median_freq == 24, temp - temp_prev, NA_real_), # restore this line once dplyr median fixed
          d_temp = ifelse(median_freq > 1, temp - temp_prev, NA_real_),
           flag_hourly_rise = ifelse(abs(d_temp) < deg | is.na(d_temp), FALSE, TRUE))
  
  return(data)
}


#' @title flag_daily_rise
#'
#' @description
#' \code{flag_daily_rise} Check for out of water (out of water vs. dry stream?)
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @param deg Threshold change in degrees for flag
#' @return Flag days with Rate of change in temperature per day above deg
#' @details
#' Rate of change in temperature per day
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_daily_rise <- function(data, deg = 5) {
  data <- ungroup(data)
  data$row <- 1:nrow(data)
  data$temp_prev <- c(NA_real_, data[2:nrow(data)-1, ]$temp)
  data$series_prev <- c(NA_character_, data[1:nrow(data)-1, ]$series_id)
  data <- data %>%
    #group_by(series_id) %>%
    mutate(series_start = ifelse(is.na(series_prev), 1, ifelse(series_id == series_prev, NA_real_, 1)),
           temp_prev = ifelse(is.na(series_start), temp_prev, NA_real_),
           d_temp = temp - temp_prev,
           flag_daily_rise = ifelse(abs(d_temp) < deg | is.na(d_temp), FALSE, TRUE))
  
  return(data)
}

#' @title convert_neg
#'
#' @description
#' \code{convert_neg} Check for excessively cold events
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @return temp Convert -1 - 0 C readings to 0
#' @details
#' Convert -1 - 0 C readings to 0
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
convert_neg <- function(data) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- data %>%
    mutate(temp = ifelse(temp < 0 & temp >= -1 & median_freq > 1, 0, temp))
  
  return(data)
}


#' @title flag_cold_obs
#'
#' @description
#' \code{flag_cold_obs} Flagging observations < -1 C
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @return Flagging observations < -1 C
#' @details
#' Flagging observations < -1 C
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_cold_obs <- function(data) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- data %>%
    mutate(flag_cold_obs = ifelse(temp < -1 & median_freq > 1, TRUE, FALSE))
  
  return(data)
}

#' @title flag_cold_days
#'
#' @description
#' \code{flag_cold_days} Flagging days < -1 C
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @return Flagging days < -1 C
#' @details
#' Flagging days < -1 C
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_cold_days <- function(data) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- data %>%
    mutate(flag_cold_days = ifelse(temp < -1 & median_freq == 1, TRUE, FALSE))
  
  return(data)
}


#' @title flag_hot_obs
#'
#' @description
#' \code{flag_hot_obs} Flagging hot observations
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @param threshold Temperature threshold to flag. Default = 35
#' @return logical Flag for hot observations
#' @details
#' Flagging observations > 30 C might not work for the database, as it isn't all small trout streams.
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_hot_obs <- function(data, threshold = 35) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- data %>%
    mutate(flag_hot_obs = ifelse(temp > threshold & median_freq > 1, TRUE, FALSE))
  
  return(data)
}

#' @title flag_hot_days
#'
#' @description
#' \code{flag_hot_days} Flagging hot days
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @param threshold Temperature threshold to flag. Default = 35
#' @return logical Flag for hot days
#' @details
#' Flagging observations > 30 C might not work for the database, as it isn't all small trout streams.
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_hot_days <- function(data, threshold = 35) {
  if(!("median_freq" %in% colnames(data, threshold))) {
    data <- obs_freq(data)
  }
  data <- data %>%
    mutate(flag_hot_days = ifelse(temp > threshold & median_freq == 1, TRUE, FALSE))
  
  return(data)
}


#' @title flag_extreme_obs
#'
#' @description
#' \code{flag_extreme_obs} Flagging extreme observations
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @param qlo Lower quantile to flag. Default = 0.001
#' @param qhi Upper quantile to flag. Default = 0.999
#' @return logical Flag for extreme observations
#' @details
#' Flagging upper and lower 5th percentiles probably wont work for the database, as those are best reverified by whoever collected/uploaded the data.
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_extreme_obs <- function(data, qlo = 0.001, qhi = 0.999) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- data %>%
    mutate(flag_extremes = ifelse(temp > quantile(temp, c(qhi), na.rm = T) | temp < quantile(temp, c(qlo) & median_freq > 1, na.rm = T), FALSE, TRUE))
  
  return(data)
}


#' @title flag_extreme_days
#'
#' @description
#' \code{flag_extreme_days} Flagging extreme days
#'
#' @param data time series data pull from the postgres database (df_values table)
#' @param qlo Lower quantile to flag. Default = 0.001
#' @param qhi Upper quantile to flag. Default = 0.999
#' @return logical Flag for extreme days
#' @details
#' Flagging upper and lower 5th percentiles probably wont work for the database, as those are best reverified by whoever collected/uploaded the data.
#' @examples
#' 
#' \dontrun{
#' 
#' }
#' @export
flag_extreme_days <- function(data, qlo = 0.001, qhi = 0.999) {
  if(!("median_freq" %in% colnames(data))) {
    data <- obs_freq(data)
  }
  data <- data %>%
    mutate(flag_extremes = ifelse(temp > quantile(temp, c(qhi), na.rm = T) | temp < quantile(temp, c(qlo) & median_freq == 1, na.rm = T), FALSE, TRUE))
  
  return(data)
}

