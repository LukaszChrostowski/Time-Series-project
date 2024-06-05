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
library(trend)

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

closed_returns <- data_microsoft$Close
squared_returns <- closed_returns^2

par(mfrow=c(1,1))
plot.ts(closed_returns)
acf(squared_returns, main = "ACF of Squared Closing Price")
# highly correlated observations
# data have a trend

# let us consider daily stock returns with log transformation
log_returns <- tq_transmute(data_microsoft,
                            select = Close,
                            mutate_fun = periodReturn,
                            period = "daily",
                            type = "log")
squared_log_returns <- log_returns$daily.returns^2

par(mfrow = c(1, 1))
ts.plot(log_returns$daily.returns, main = "Stock Daily Returns over Time",
        ylab = "Daily Return")

par(mfrow = c(2, 1))
acf(log_returns$daily.returns, main = "ACF of Daily Returns")
pacf(log_returns$daily.returns, main = "PACF of Daily Returns")

par(mfrow = c(2, 1))
acf(squared_log_returns, main = "ACF of Squared Daily Returns")
pacf(squared_log_returns, main = "PACF of Squared Daily Returns")
# Conclusions: Observations are highly correlated

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
# adf.test(log_returns$daily.returns)
# kpss.test(log_returns$daily.returns, null = c("Level", "Trend"), lshort = TRUE)

# Trend
ts.plot(log_returns)
model <- lm(data = log_returns, daily.returns ~ Date)
# slowly increasing trend is visible but not significant from statistical point of view

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

############# ARIMA #################

# Since we know that we have (?) trend in data - we have to remove it using autoarima
arima_fit <- auto.arima(log_returns$daily.returns)
# ARIMA(2,0,1) with non-zero mean 
residuals_arima <- residuals(arima_fit)


############# GARCH #################
# function to fit the model
do_garch <- function(garch.model, garch.r, garch.s, data, dist.model, ar = 0, ma = 0) {
  # Specify the GARCH model - here we use a basic GARCH(1,1) model
  garch_spec <- ugarchspec(variance.model=list(model=garch.model,
                                               garchOrder=c(garch.r, garch.s)),
                           mean.model=list(armaOrder=c(ar, ma)),
                           distribution.model = dist.model)
  # just for now - probably needed to find more efficient one
  
  # Fit the GARCH model to the differenced volume series
  # garch_fit <- ugarchfit(spec = garch_spec, data = log_returns$daily.returns)
  garch_fit <- ugarchfit(spec = garch_spec,
                         data = data)
  
  ############## RESIDUAL ANALYSIS #################
  residuals <- residuals(garch_fit, standardize = TRUE)
  squared_residuals <- residuals^2
  
  list(garch_spec = garch_spec,
       fit = garch_fit,
       residuals = residuals,
       squared_residuals = squared_residuals)
}
############# Assuming TREND in DATA -> basing on residuals from ARIMA #################
# (1, 1)
garch <- do_garch(garch.model = "sGARCH",
                  garch.r = 1,
                  garch.s = 1,
                  data = residuals_arima,
                  dist.model = "norm")
print(garch$fit)

# Conclusion from fit
# 1 No serial correlation between garch residuals
# Akaike       -5.5615
# Bayes        -5.5474

par(mfrow = c(1, 1))
plot.ts(garch$residuals, main = "Residuals of GARCH model", ylab = "Residuals")
plot.ts(garch$squared_residuals, main = "Squared Residuals of GARCH model", ylab = "Squared Residuals")

# ACF analysis for residuals
# This will allow you to check the adequacy of the GARCH(1,1) model by examining whether the residuals
# and their squared values exhibit any significant autocorrelation. If the ACF and PACF plots show no
# significant autocorrelation, and have a N(0,1) white noise behaviour, the model is considered a good fit for the data.
par(mfrow = c(2, 2))
acf(garch$residuals, main = "ACF of Standardized Residuals")
pacf(garch$residuals, main = "PACF of Standardized Residuals")
acf(garch$squared_residuals, main = "ACF of Squared Standardized Residuals")
pacf(garch$squared_residuals, main = "PACF of Squared Residuals")
par(mfrow = c(1, 1))


############# Assuming NO TREND in DATA #################

# Define the grid for GARCH orders (excluding the (0,0) point)
garch_orders <- expand.grid(garch.r = 0:2, garch.s = 0:2)
garch_orders <- garch_orders[!(garch_orders$garch.r == 0 & garch_orders$garch.s == 0), ]

# Define the grid for ARMA orders (0 to 2)
arma_orders <- expand.grid(ar = 0:2, ma = 0:2)

# Initialize variables to store the results
results <- data.frame(garch.r = integer(),
                      garch.s = integer(),
                      ar = integer(),
                      ma = integer(),
                      AIC = numeric(),
                      BIC = numeric(),
                      SHIB = numeric(),
                      HAN = numeric(),
                      stringsAsFactors = FALSE)

# Nested loops to fit the GARCH models with different GARCH and ARMA orders
for (i in 1:nrow(garch_orders)) {
  garch.r <- garch_orders$garch.r[i]
  garch.s <- garch_orders$garch.s[i]
  
  for (j in 1:nrow(arma_orders)) {
    ar <- arma_orders$ar[j]
    ma <- arma_orders$ma[j]
    
    # Fit the GARCH model with ARMA orders
    try({
      model <- do_garch(garch.model = "sGARCH", 
                        garch.r = garch.r, 
                        garch.s = garch.s, 
                        data = log_returns$daily.returns, 
                        dist.model = "norm", 
                        ar = ar, 
                        ma = ma)
      
      # Extract the AIC and BIC
      aic <- infocriteria(model$fit)[1]
      bic <- infocriteria(model$fit)[2]
      shib <- infocriteria(model$fit)[3]
      han <- infocriteria(model$fit)[4]
      
      # Store the results
      results <- rbind(results, data.frame(garch.r = garch.r, garch.s = garch.s,
                                           ar = ar, ma = ma,
                                           AIC = aic, BIC = bic,
                                           SHIB = shib, HAN = han))
    }, silent = TRUE)
  }
}

# Find the model with the lowest AIC and BIC
best_aic_model <- results[which.min(results$AIC), ]
best_bic_model <- results[which.min(results$BIC), ]
best_shib_model <- results[which.min(results$SHIB), ]
best_han_model <- results[which.min(results$HAN), ]

# Print the best models
print("Best ARMA model by AIC:")
print(best_aic_model)
print("Best ARMA model by BIC:")
print(best_bic_model)
print("Best ARMA model by SHIB:")
print(best_shib_model)
print("Best ARMA model by HAN:")
print(best_han_model)

########## FINAL MODEL 
garch <- do_garch(garch.model = "sGARCH", 
                  garch.r = 1, 
                  garch.s = 1, 
                  data = log_returns$daily.returns, 
                  dist.model = "norm", 
                  ar = 1, 
                  ma = 2)

par(mfrow = c(1, 1))
plot.ts(garch$residuals, main = "Residuals of GARCH model", ylab = "Residuals")
plot.ts(garch$squared_residuals, main = "Squared Residuals of GARCH model", ylab = "Squared Residuals")

# ACF analysis for residuals
# This will allow you to check the adequacy of the GARCH(1,1) model by examining whether the residuals
# and their squared values exhibit any significant autocorrelation. If the ACF and PACF plots show no
# significant autocorrelation, and have a N(0,1) white noise behaviour, the model is considered a good fit for the data.
par(mfrow = c(2, 2))
acf(garch$residuals, main = "ACF of Standardized Residuals")
pacf(garch$residuals, main = "PACF of Standardized Residuals")
acf(garch$squared_residuals, main = "ACF of Squared Standardized Residuals")
pacf(garch$squared_residuals, main = "PACF of Squared Standardized Residuals")
par(mfrow = c(1, 1))

########## FORECASTING #################################

