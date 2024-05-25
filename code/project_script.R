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

# data_microsoft <- read.csv("C:/Users/Miguel Silva/Desktop/universidade/4th year/2nd semester/TS/project/data_microsoft/Microsoft_Stock.csv")
data_microsoft <- read.csv("~/Desktop/Time-Series-project/data/Microsoft_Stock.csv")
data_microsoft$Date <- as.Date(data_microsoft$Date, format="%m/%d/%Y %H:%M:%S")
plot <- chartSeries(data_microsoft)

# Series is daily stock return
series <- diff(data_microsoft$Close)/data_microsoft$Close[-length(data_microsoft$Close)]
series_df <- c(NA, series)
data_microsoft$stock_return <- series_df
ts.plot(series)

d_series = diff(series)
ts.plot(d_series)

###############################################
#.     INITIAL EXPLORATION AND ANALYSIS       #
###############################################

summary_stats <- c(
    mean_volume = mean(series),
    median_volume = median(series),
    variance_volume = var(series),
    sd_volume = sd(series),
    min_volume = min(series),
    max_volume = max(series)
  )
print(summary_stats)

# Plot the Volume over time
ggplot(data_microsoft, aes(x=Date, y=stock_return)) +
  geom_line(color="blue") +
  labs(title="Volume over Time",
       x="Date",
       y="Volume") +
  theme_minimal()

# histogram and density of Volume
p2 <- ggplot(data_microsoft) 

p2 + geom_histogram(aes(x=stock_return, y=..density..), bins = 100, color="steelblue",
                    fill="grey", size=1) +
  stat_function(fun = dnorm, args = list(mean = mean(data_microsoft$stock_return, na.rm = T),
                                         sd = sd(data_microsoft$stock_return, na.rm = T)), size=1)


########### STATIONARITY ##########

# stationarity - series
adf.test(series)
kpss.test(series, null = c("Level", "Trend"), lshort = TRUE)

# stationarity - first diffs
adf.test(d_series)#stationary 
kpss.test(d_series, null = c("Level", "Trend"), lshort = TRUE)

###############################################
#               ACF and PACF                  #
###############################################

acf_series = acf(series, lag.max=70, plot=FALSE)
plot(acf_series, type="h", xlab="lag", col='red')

pacf_series = pacf(series, lag.max=70, plot=FALSE)
plot(pacf_series, type="h", xlab="lag", col='red')

###############################################
#                DIFFERENTIAL                 #
###############################################

ts.plot(d_series)

acf_d_series = acf(d_series, lag.max=70, plot=FALSE)
plot(acf_d_series, type="h", xlab="lag", col='red')

pacf_d_series = pacf(d_series, lag.max=70, plot=FALSE)
plot(pacf_d_series, type="h", xlab="lag", col='red')


###############################################
#                  MODELING                   #
###############################################

########### VOLATILITY #######################
rolling_vol <- rollapply(data_microsoft$stock_return, width = 30,
                         FUN = sd, na.rm = TRUE, align='right')
vol <- data.frame(index(rolling_vol), rolling_vol)
colnames(vol) <- c("date", "volatility")
p3 <- ggplot(vol, aes(x=date, y=volatility))
p3 +
  geom_line( color="steelblue")

# sarima
sarima(series, 0, 0, 2, 0, 0, 0, 0)

############# GARCH #################

# Specify the GARCH model - here we use a basic GARCH(1,1) model
garch_spec <- ugarchspec(variance.model=list(model="sGARCH",
                                             garchOrder=c(1,1)),
                         mean.model=list(armaOrder=c(0,0)),
                         distribution.model = "norm")

# Fit the GARCH model to the differenced volume series
garch_fit <- ugarchfit(spec = garch_spec, data = series)

# Print the summary of the model fit
print(garch_fit)
