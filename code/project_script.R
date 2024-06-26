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
# stationarity - original time series
adf.test(data_microsoft$Close)
kpss.test(data_microsoft$Close, null = c("Level", "Trend"), lshort = TRUE)
# stationarity - log returns
adf.test(log_returns$daily.returns)
kpss.test(log_returns$daily.returns, null = c("Level", "Trend"), lshort = TRUE)
# stationarity - squared_log_returns
adf.test(squared_log_returns)
kpss.test(squared_log_returns, null = c("Level", "Trend"), lshort = TRUE)

# Trend
ts.plot(log_returns)
model <- lm(data = log_returns, daily.returns ~ Date)

# plot
ggplot(log_returns, aes(x = Date, y = daily.returns)) +
  geom_line(color = "blue") +  # Plot the time series data
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Add the linear trend line
  labs(title = "Trend Plot for Daily Returns",
       x = "Date",
       y = "Value") +
  theme_minimal()
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

############# GARCH #################
# function to fit the model
do_garch <- function(garch.model,
                     garch.r, # parameter for GARCH model
                     garch.s, # parameter for GARCH model
                     data,
                     dist.model,
                     ar = 0, # parameter for AR part
                     ma = 0) # parameter for MA part
  {
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

# shapiro.test(as.vector(garch$residuals))

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
# Since we found the most suitable model for our time series.
# now we can spllit data and do forecast to evaluate our model on test data
# Train-test split (80-20)
set.seed(123)
n <- length(log_returns$daily.returns)
train_size <- floor(0.8 * n)
train_data <- log_returns$daily.returns[1:train_size]
test_data <- log_returns$daily.returns[(train_size + 1):n]
test_dates <- log_returns$Date[(train_size + 1):n]

# Fit GARCH model on training data
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(1, 2)),
                   distribution.model = "norm")
fit <- ugarchfit(spec, train_data, solver = "hybrid")

# Forecast based on test data
forecast_horizon <- length(test_data)
garch_forecast <- ugarchforecast(fit, n.ahead = forecast_horizon)


# FORECAST OF DAILY RETURNS ###################

# Extract forecasted values and confidence intervals
forecasted_values <- fitted(garch_forecast)
forecasted_sigma <- sigma(garch_forecast)
forecasted_upper <- forecasted_values + 1.96 * forecasted_sigma
forecasted_lower <- forecasted_values - 1.96 * forecasted_sigma

# Create data frame for plotting
forecast_df <- data.frame(Date = test_dates,
                          Forecasted_Returns = forecasted_values,
                          Upper_CI = forecasted_upper,
                          Lower_CI = forecasted_lower)
colnames(forecast_df) <- c("Date", "Forecasted_Returns", "Upper_CI", "Lower_CI")

# Plot forecasted returns with confidence intervals
ggplot() +
  geom_line(data = log_returns, aes(x = Date, y = daily.returns), color = "blue") +
  geom_line(data = forecast_df, aes(x = Date, y = Forecasted_Returns), color = "red") +
  geom_ribbon(data = forecast_df, aes(x = Date, ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2, fill = "red") +
  labs(title = "GARCH Model Forecast with Confidence Intervals",
       x = "Date",
       y = "Daily Returns") +
  theme_minimal()

# FORECAST OF VOLATILITY (VARIABILITY) ######################

# Initialize variables
test_start <- train_size + 1
rolling_forecast_sigma <- numeric(length(test_data))  # Adjust size according to rolling window

# Rolling forecast on test data
for (i in test_start:n) {
  train_data <- log_returns$daily.returns[1:i]

  # Fit the GARCH model
  spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(1, 2)),
                     distribution.model = "norm")

  fit <- ugarchfit(spec, train_data, solver = "hybrid")

  # Forecast one step ahead
  forecast <- ugarchforecast(fit, n.ahead = 1)

  # Store the forecasted sigma
  rolling_forecast_sigma[i - test_start + 1] <- sigma(forecast)
}

# Create data frame for plotting
rolling_forecast_dates <- log_returns$Date[test_start:n]
rolling_forecast_df <- data.frame(Date = rolling_forecast_dates, Predicted_Volatility = rolling_forecast_sigma)

# Extract actual values for comparison
actual_volatility <- log_returns$daily.returns[test_start:n]

# Plot predicted volatility with actual values
ggplot() +
  geom_line(data = rolling_forecast_df, aes(x = Date, y = Predicted_Volatility), color = "red") +
  geom_line(aes(x = rolling_forecast_dates, y = actual_volatility), color = "blue") +
  labs(title = "One-Step Ahead Rolling Forecast of Volatility with Actual Values",
       x = "Date",
       y = "Volatility") +
  theme_minimal()
