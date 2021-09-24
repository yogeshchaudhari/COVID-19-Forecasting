
#install.packages("TSA")
library(forecast)
library(TSA)
library(zoo)
library(Metrics)
library(stats)
library(stringr)
#install.packages("WaveletANN")
library(WaveletANN)
#library(WaveletArima)


zero_reference = 18364 #Difference between 1970-01-01 to 2020-02-29
skip = 0
ts_size = 368     # Number of Sample/Timesteps in entire time series
tr_size = 228 #334 # Number of Samples/Timesteps in Training set 236
te_size = ts_size - tr_size # Number of Samples/Timesteps in Training set 140
horizons = 14  # Number of days to forecast ahead in the future
kt = 1       # Retraining Frequency (No. of Samples)
rt = 0         # Total Number of Retrainings Occured
p = 1          # Non-seasonal Auto-Regressive Order
d = 1          # Non-seasonal Differencing/
q = 0          # Non-seasonal Moving Average Order
PP = 3         # Seasonal Auto-Regressive Order
DD = 1         # Seasonal Differencing
QQ = 1         # Seasonal Moving Average Order
s = 7          # Seasonal Period 

# Pull in the dataset from given URL
covid_data <- read.csv(file = "state_daily_deaths.csv")
cases_data <- read.csv(file = "state_daily_cases.csv")
state_list <- read.csv(file = "state_X_time.csv")$X
state_num = 0

state = "Ohio"
response_variable = ifelse(covid_data[[state]] <=0, 1, covid_data[[state]])
exog_variable = ifelse(cases_data[[state]] <=0, 1, cases_data[[state]])
smape_scores = matrix(data = NA, nrow = horizons, ncol = 22)
covid_ts <- zoo(response_variable, seq(as.Date("2020-04-13"), as.Date("2021-04-15"), by = "days"), frequency = 1)

country_test_data <- window(covid_ts, start = as.Date((tr_size+1) + zero_reference), end = as.Date(ts_size + zero_reference), frequency = 1)
xyz<- window(covid_ts, start = as.Date((tr_size+1) + zero_reference), end = as.Date(ts_size + zero_reference), frequency = 1)
country_forecast_matrix <- matrix(data = c(0), nrow = te_size, ncol = horizons)
forecast_matrix <- matrix(data = c(0), nrow = 14, ncol = te_size)
for(state in state_list){
  print(state)
  state_num = state_num + 1
  if(state_num != 46) {
    response_variable = ifelse(covid_data[[state]] <=0, 1, covid_data[[state]])
    exog_variable = ifelse(cases_data[[state]] <=0, 1, cases_data[[state]])
    smape_scores = matrix(data = NA, nrow = horizons, ncol = 22)
    covid_ts <- zoo(response_variable, seq(as.Date("2020-04-13"), as.Date("2021-04-15"), by = "days"), frequency = 1)
    exog_ts <- zoo(exog_variable, seq(as.Date("2020-04-13"), as.Date("2021-04-15"), by = "days"), frequency = 1)
    for (i in 1:te_size) {
      start = i + zero_reference
      model_data <- window(log(covid_ts), start = as.Date(start), end = as.Date(start + tr_size - 1), frequency = 1)
      exog_data <- window(log(exog_ts), start = as.Date(start), end = as.Date(start + tr_size - 1), frequency = 1)
      if (i == 1) {
        arima_model <- Arima(y = model_data, order = c(p, d, q), seasonal = list(order = c(PP, DD, QQ), period = s), method = "CSS")
      } else { 
        if ((i %% kt) == 0) { #if iteration mod retrain frequency is 0, then it is time to retrain the model with new training data
          rt = rt + 1
          arima_model <- Arima(y = model_data, order = c(p, d, q), seasonal = list(order = c(PP, DD, QQ), period = s), method = "CSS")
        } else {
          arima_model <- Arima(y = model_data, model = arima_model, xreg=exog_data)
          
        }
      }
      #future <- forecast(arima_model, h = horizons,level = 0)
      newxreg = window(log(exog_ts), start = as.Date(start - horizons + tr_size), end = as.Date(start + tr_size - 1), frequency = 1)
      future <- forecast(arima_model, h = horizons, level = 0)
      
      forecasts <- future$mean
      for (k in 1:horizons) {
        if (k + i - 1 <= te_size) { # Condition to prevent overflow at the end of the test set;
          forecast_matrix[k, k+i-1] = forecasts[k]
        }
      }
    }
    test_data <- window(covid_ts, start = as.Date((tr_size+1) + zero_reference), end = as.Date(ts_size + zero_reference), frequency = 1)
    country_test_data = country_test_data + test_data
  }
  forecast_matrix2 <- t(forecast_matrix)
  forecast_matrix2 <- exp(forecast_matrix2)
  country_forecast_matrix = country_forecast_matrix + forecast_matrix2
}

country_test_data = country_test_data - xyz

country_forecast_matrix <- t(country_forecast_matrix)
for(i in 2:horizons){
  for(j in 1:i-1){
    country_forecast_matrix[i,j] = country_forecast_matrix[i-1,j]
  }
  
}
country_forecast_matrix <- t(country_forecast_matrix)
View(covid_data)

for (i in 1:horizons) {
  smape_scores[i, p] <- mean(abs((country_test_data - country_forecast_matrix[, i]) / (country_test_data + country_forecast_matrix[, i]))) * 200 #smape formula
  #smape_scores[i, p] <- mean(abs(test_data - forecast_matrix2[, i]))
  print(paste0("h=",i,",smape=",smape_scores[i,p]))
}
forecast_matrix2 <- t(forecast_matrix)
forecast_matrix2 <- exp(forecast_matrix2)

write.table(country_forecast_matrix,file="state_data_country_forecast_matrix.csv", sep = ",") # keeps the rownames
write.table(country_test_data,file="state_data_country_test_data.csv", sep = ",") # keeps the rownames

state_num = 0
forecast_matrix <- matrix(data = c(-1.2), nrow = ts_size, ncol = 52)
for(state in state_list){
  state_num = state_num + 1
  if(state_num != 46) {
    response_variable = ifelse(covid_data[[state]] <=0, 1, covid_data[[state]])
    covid_ts <- zoo(response_variable, seq(as.Date("2020-04-13"), as.Date("2021-04-15"), by = "days"), frequency = 1)
    in_sample_data <- window(log(covid_ts), start = as.Date(1 + zero_reference), end = as.Date(544 + zero_reference), frequency = 1)
    in_sample_model <- Arima(y=in_sample_data, order = c(p,d,q), seasonal = list(order=c(PP,DD,QQ),period=s), method="CSS")
    in_ts <- as.ts(in_sample_data)
    fitted_values <- in_ts - abs(in_sample_model$residuals)
    forecast_matrix[,state_num] = exp(fitted_values)
      
  }
  
}


forecast_matrix2 <- t(forecast_matrix)
forecast_matrix2 <- exp(forecast_matrix2)

write.table(forecast_matrix2,file="state_data_in_sample_forecast_matrix2.csv", sep = ",") # keeps the rownames

#In-sample Analysis
xreg_lag = 14
in_sample_data <- window(log(covid_ts), start = as.Date(1 + zero_reference), end = as.Date(544 + zero_reference), frequency = 1)
in_sample_exog <- window(log(exog_ts), start = as.Date(1 + zero_reference - xreg_lag), end = as.Date(544 + zero_reference - xreg_lag), frequency = 1)
in_ts <- as.ts(in_sample_data)
ex_ts <- as.ts(in_sample_exog)

in_sample_model <- Arima(y=in_sample_data, order = c(p,d,q), seasonal = list(order=c(PP,DD,QQ),period=s), method="CSS", xreg = in_sample_exog)
fitted_values <- in_ts - abs(in_sample_model$residuals)
in_sample_smape <- mean(abs(exp(in_ts) - exp(fitted_values)) / (abs(exp(in_ts)) + abs(exp(fitted_values)))) * 200
print(xreg_lag)
print(in_sample_smape)

plot(exp(fitted_values))
plot(exp(in_sample_data))
length(fitted_values)

write.table(exp(fitted_values) ,file="fitted_values_new_dataset.csv", sep = ",") # keeps the rownames
sum = 0.0
e <- exp(in_ts) - exp(fitted_values)
cnt = 0
for(i in 1:length(in_ts)){
  denom <- abs(exp(in_ts[i])) + abs(exp(fitted_values[i]))
  if(denom != 0){
    sum = sum + (abs(e[i]) / denom)
    cnt = cnt + 1
  }
}
print("In-Sample SMAPE:")
print(200 * sum/cnt)


#write.table(smape_scores,file="D:\\Courses\\7200-MSProject\\R\\test.txt") # keeps the rownames
#unlink("test.txt")
