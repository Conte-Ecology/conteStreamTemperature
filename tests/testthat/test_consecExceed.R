context("Consecutive Exceedance Days")

library(dplyr)
library(lubridate)

featureid <- 1:4
year <- 2000:2003
month <- 6
day <- 1:5
df <- expand.grid(featureid, year, month, day, stringsAsFactors = FALSE)
names(df) <- c("featureid", "year", "month", "day")

df$tempPredicted <- c(rnorm(n = length(year)*length(month)*length(day), 22, 2),
                      rnorm(n = length(year)*length(month)*length(day), 22, 2),
                      rnorm(n = length(year)*length(month)*length(day), 22, 2),
                      rnorm(n = length(year)*length(month)*length(day), 22, 2))
df <- df %>%
  dplyr::mutate(date = ymd(paste0(year, "-", month, "-", day))) %>%
  dplyr::arrange(., featureid, date)

#derivedfeatureidMetrics <- dplyr::select(df, featureid) %>%
  #distinct

# test result size
test_that("the number of rows produced is equal to # of groups", {  
  df_consec <- group_by(df, featureid, year) %>%
    do(calcConsecExceed(., 22))
  
  expect_equal(dim(df_consec)[1], length(unique(df$featureid)) * length(unique(df$year)))
})

# test that individual NA ignored
i <- runif(n = 3, 1, length(df$tempPredicted))
df$tempPredicted[i] <- NA

test_that("individual NA are ignored by summary metrics", {
  df_consec <- group_by(df, featureid, year) %>%
    do(calcConsecExceed(., 22))
  
    expect_equal(dim(df_consec)[1], length(unique(df$featureid)) * length(unique(df$year)))
})

# test the class

# test the type of objects in the group_by

# test that returns a 0 if all values <= threshold
test_that("a zero is returned if all values below threshold", {
  df <- data.frame(tempPredicted = runif(10, 2, 20))
  consecExceed <- calcConsecExceed(df, 22)
  expect_equal(consecExceed$maxConsec_22, 0)
  expect_equal(consecExceed$meanConsec_22, 0)
})

# test if all values > threshold
test_that("a zero is returned if all values below threshold", {
  df <- data.frame(tempPredicted = runif(10, 23, 30))
  consecExceed <- calcConsecExceed(df, 22)
  expect_equal(consecExceed$maxConsec_22, length(df$tempPredicted))
  expect_equal(consecExceed$meanConsec_22, length(df$tempPredicted))
})

# test if all values are NA
test_that("NA is returned when all values are NA", {
  df <- data.frame(tempPredicted = rep(NA, times = 10))
  consecExceed <- calcConsecExceed(df, 22)
  expect_true(is.na(consecExceed$maxConsec_22))
  expect_true(is.na(consecExceed$meanConsec_22))
  expect_warning(calcConsecExceed(df, 22))
})

# test for correct mean and max
test_that("returns correct mean and max consecutive lengths", {
  df <- data.frame(tempPredicted = c(10, 10, 25, 25, 25, 10, 10, 25, 25, 10, 25))
  consecExceed <- calcConsecExceed(df, 22)
  expect_equal(consecExceed$maxConsec_22, 3)
  expect_equal(consecExceed$meanConsec_22, 2)
})


# get max/mean run of each featureid, year

