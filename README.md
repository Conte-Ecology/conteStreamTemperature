This package is used to clean, prepared, and analyze stream temperature data collected at the daily or sub-daily level. Data and model predictions can be visualized using `ggplot2` and our helper functions.

An example of how to use this package can be found in the vignettes...

to run the vignettes...

A full example where we use this package can be found on the Conte-Ecology GitHub package in the `conteStreamTemperature_northeast` repository.

To install and load this package you will need to have `RTools` installed and the R package `remotes` (or `devtools` which calls to `remotes`) installed. This package provides the `install_github` function so you can install our package directly from GitHub as shown below:

```
remotes::install_github("Conte-Ecology/conteStreamTemperature")
library(conteStreamTemperature)
```

This is the project folder for the stream temperature work underway at the USGS S.O. Conte Anadromous Fish Research Center in Turners Falls, MA.

The stream temperature model estimates effects of landscape variables (% forest cover, % agriculture, elevation, etc.) and time varying variables (solar radiation, snow-water equivalent, air temperature, etc.) on daily stream water temperature. For each site/year combination, the estimates are limited to the times of the year where air temperature and water temperature are synchronized.
