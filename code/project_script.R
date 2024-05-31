###############################################
#             IMPORTING DATA                  #
###############################################

library(stats)
library(astsa)
library(tseries)
library(dplyr)
library(ggplot2)
library(rugarch)
library(rmgarch)
library(quantmod)
library(tidyquant)
library(dplyr)
library(forecast)

# data_microsoft <- read.csv("C:/Users/Miguel Silva/Desktop/universidade/4th year/2nd semester/TS/project/data_microsoft/Microsoft_Stock.csv")
data_microsoft <- read.csv("~/Desktop/Time-Series-project/data/Microsoft_Stock.csv")
data_microsoft$Date <- as.Date(data_microsoft$Date, format="%m/%d/%Y %H:%M:%S")
plot <- chartSeries(data_microsoft)

# Series is daily stock return
# series <- diff(data_microsoft$Close)/data_microsoft$Close[-length(data_microsoft$Close)]
# series_df <- c(NA, series)
# data_microsoft$stock_return <- series_df
# ts.plot(series)

###############################################
#  DATA PREPARATION AND ACF ANALYSIS         #
###############################################
log_returns <- tq_transmute(data_microsoft,
                            select = Close,
                            mutate_fun = periodReturn,
                            period = "daily",
                            type = "log")

par(mfrow = c(2, 1))
ts.plot(log_returns) # if trend is visible (it is!), model should be like sarima + garch
ts.plot(log_returns$daily.returns, main = "Stock Daily Returns over Time",
        ylab = "Daily Return")

acf(log_returns)
pacf(log_returns)

# d_series = diff(series)
# ts.plot(d_series)

###############################################
#.     INITIAL EXPLORATION AND ANALYSIS       #
###############################################

summary_stats <- c(
    mean_daily_return = mean(log_returns$daily.returns),
    median_vdaily_return = median(log_returns$daily.returns),
    variance_daily_return = var(log_returns$daily.returns),
    sd_daily_return = sd(log_returns$daily.returns),
    min_daily_return = min(log_returns$daily.returns),
    max_daily_return = max(log_returns$daily.returns)
  )
print(summary_stats)

# Plot the Volume over time
ggplot(log_returns, aes(x=Date, y=daily.returns)) +
  geom_line(color="blue") +
  labs(title="Daily return over Time",
       x="Date",
       y="Daily Return") +
  theme_minimal()

# histogram and density of Daily Stock Returns
p2 <- ggplot(log_returns) 

p2 + geom_histogram(aes(x=daily.returns, y=..density..), bins = 100, color="steelblue",
                    fill="grey", size=1) +
  stat_function(fun = dnorm, args = list(mean = mean(log_returns$daily.returns, na.rm = T),
                                         sd = sd(log_returns$daily.returns, na.rm = T)), size=1) +
  xlab("Daily Returns") +
  ggtitle("Histogram and density of Daily Stock Returns")
# looks like normal distribution


########### STATIONARITY ##########

# stationarity - series
adf.test(log_returns$daily.returns)
kpss.test(log_returns$daily.returns, null = c("Level", "Trend"), lshort = TRUE)

# stationarity - first diffs
# adf.test(d_series) #stationary 
# kpss.test(d_series, null = c("Level", "Trend"), lshort = TRUE)

###############################################
#               ACF and PACF                  #
###############################################

# acf_series = acf(series, lag.max=70, plot=FALSE)
# plot(acf_series, type="h", xlab="lag", col='red')
# 
# pacf_series = pacf(series, lag.max=70, plot=FALSE)
# plot(pacf_series, type="h", xlab="lag", col='red')

###############################################
#                DIFFERENTIAL                 #
###############################################

# ts.plot(d_series)
# 
# acf_d_series = acf(d_series, lag.max=70, plot=FALSE)
# plot(acf_d_series, type="h", xlab="lag", col='red')
# 
# pacf_d_series = pacf(d_series, lag.max=70, plot=FALSE)
# plot(pacf_d_series, type="h", xlab="lag", col='red')


###############################################
#                  MODELING                   #
###############################################

########### VOLATILITY #######################
rolling_vol <- rollapply(log_returns$daily.returns, width = 30,
                         FUN = sd, na.rm = TRUE, align='right')
vol <- data.frame(index(rolling_vol), rolling_vol)
colnames(vol) <- c("date", "volatility")
p3 <- ggplot(vol, aes(x=date, y=volatility))
p3 + 
  geom_line( color="steelblue") +
  labs(title="Daily return volatility over Time",
       x="Day",
       y="Volatility")

############# SARIMA #################
# STEP 1 - Find the most suitable sarima model
# sarima(log_returns$daily.returns, 0, 0, 2, 0, 0, 0, 0)

# autoarima
arima_fit <- auto.arima(log_returns$daily.returns)
residuals_arima <- residuals(arima_fit)

############# GARCH #################
# STEP 2 - Find the most suitable garch model in respect to choosen sarima model (basing on residuals from SARIMA model)
# Specify the GARCH model - here we use a basic GARCH(1,1) model
garch_spec <- ugarchspec(variance.model=list(model="sGARCH",
                                             garchOrder=c(1,1)),
                         mean.model=list(armaOrder=c(0,0)),
                         distribution.model = "norm")
# just for now - probably needed to find more efficient one

# Fit the GARCH model to the differenced volume series
# garch_fit <- ugarchfit(spec = garch_spec, data = log_returns$daily.returns)
garch_fit <- ugarchfit(spec = garch_spec, data = residuals_arima)

# Print the summary of the model fit
print(garch_fit)

############## RESIDUAL ANALYSIS #################
residuals <- residuals(garch_fit)
plot.ts(residuals)
squared_residuals <- residuals^2
plot.ts(squared_residuals, main = "Squared Residuals of GARCH model", ylab = "Squared Residuals")

# normality test for residuals
shapiro.test(as.numeric(residuals))

# ACF analysis for residuals
acf(squared_residuals, main = "ACF of squared residuals")
