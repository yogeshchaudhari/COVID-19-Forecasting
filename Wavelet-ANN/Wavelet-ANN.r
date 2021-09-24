
#install.packages("TSA")
library(forecast)
library(TSA)
library(zoo)
library(Metrics)
library(stats)
#install.packages("WaveletANN")
library(WaveletANN)
#library(WaveletArima)


zero_reference = 18320 #Difference between 1970-01-01 to 2020-02-29
skip = 38
ts_size = 544     # Number of Sample/Timesteps in entire time series
tr_size = 222 #334 # Number of Samples/Timesteps in Training set 236
te_size = ts_size - tr_size # Number of Samples/Timesteps in Training set 140
horizons = 14  # Number of days to forecast ahead in the future
kt = 1       # Retraining Frequency (No. of Samples)
rt = 0         # Total Number of Retrainings Occured

# Pull in the dataset from given URL
covid_data <- read.csv(file = "2021-08-26-21-42-06-OWID.csv")

dates = covid_data$date
#Convert reponse column into "zoo" object
response_variable = covid_data$new_deaths
response_variable = ifelse(covid_data$new_deaths ==0, 1, covid_data$new_deaths)

dates[skip + tr_size + 1]

covid_ts <- zoo(response_variable, seq(as.Date("2020-01-22"), as.Date("2021-08-25"), by = "days"), frequency = 1)

smape_scores = matrix(data = NA, nrow = horizons, ncol = 22)

forecast_matrix <- matrix(data = c(-1.2), nrow = horizons, ncol = te_size)

for (i in 1:te_size) {
  print(i + " of " + te_size)
  start = i + zero_reference
  model_data <- window(log(covid_ts), start = as.Date(start), end = as.Date(start + tr_size - 1), frequency = 1)
  
  WaveletForecast<-WaveletFittingann(ts=ts(model_data),Waveletlevels=floor(log(length(model_data))),
                                     boundary='periodic',FastFlag=TRUE,nonseaslag,seaslag,NForecast=14)
  
  forecasts <- WaveletForecast$Finalforecast
  for (k in 1:horizons) {
    if (k + i - 1 <= te_size) { # Condition to prevent overflow at the end of the test set;
      forecast_matrix[k, k+i-1] = forecasts[k]
    }
  }
}

for(i in 2:horizons){
  for(j in 1:i-1){
    forecast_matrix[i,j] = forecast_matrix[i-1,j]
  }
}
forecast_matrix2 <- t(forecast_matrix)
forecast_matrix2 <- exp(forecast_matrix2)
test_data <- window(covid_ts, start = as.Date((tr_size+1) + zero_reference), end = as.Date(ts_size + zero_reference), frequency = 1)

for (i in 1:horizons) {
  smape_scores[i, p] <- mean(abs((test_data - forecast_matrix2[, i]) / (test_data + forecast_matrix2[, i]))) * 200 #smape formula
  print(paste0("h=",i,",smape=",smape_scores[i,p]))
}
for (i in 1:horizons) {
  smape_scores[i, p] <- mean(abs(test_data - forecast_matrix2[, i]))
  print(paste0("h=",i,",mae=",smape_scores[i,p]))
}