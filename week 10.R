#install.packages("TTR")
#install.packages("forecast")


library("TTR")
library("forecast")

###################### Reading Time Series Data ########################

kings <- scan("http://robjhyndman.com/tsdldata/misc/kings.dat",skip=3)
kingstimeseries <- ts(kings)
kingstimeseries
plot.ts(kingstimeseries)



souvenir <- scan("http://robjhyndman.com/tsdldata/data/fancy.dat")
souvenirtimeseries <- ts(souvenir, frequency=12, start=c(1987,1))
souvenirtimeseries
plot.ts(souvenirtimeseries)

logsouvenirtimeseries <- log(souvenirtimeseries)
plot.ts(logsouvenirtimeseries)



births <- scan("http://robjhyndman.com/tsdldata/data/nybirths.dat")
birthstimeseries <- ts(births, frequency=12, start=c(1946,1))
birthstimeseries
plot.ts(birthstimeseries)

###################### Decomposition Time Series ########################

#####  Decomposing Non-Seasonal Data ####

kingstimeseriesSMA3 <- SMA(kingstimeseries,n=3)
plot.ts(kingstimeseriesSMA3)

kingstimeseriesSMA8 <- SMA(kingstimeseries,n=8)
plot.ts(kingstimeseriesSMA8)

####### Decomposing Seasonal Data ######

birthstimeseriescomponents <- decompose(birthstimeseries)

birthstimeseriescomponents$seasonal 

plot(birthstimeseriescomponents)

birthstimeseriescomponents <- decompose(birthstimeseries)
birthstimeseriesseasonallyadjusted <- birthstimeseries - birthstimeseriescomponents$seasonal

plot(birthstimeseriesseasonallyadjusted)

##################### Forecasts using Smoothing ##########################

####### Simple Exponential Smoothing #######

rain <- scan("http://robjhyndman.com/tsdldata/hurst/precip1.dat",skip=1) 
rainseries <- ts(rain,start=c(1813))
plot.ts(rainseries)

rainseriesforecasts <- HoltWinters(rainseries, beta=FALSE, gamma=FALSE)
rainseriesforecasts 

rainseriesforecasts$fitted
rainseriesforecasts$x

plot(rainseriesforecasts)

rainseriesforecasts$SSE

HoltWinters(rainseries, beta=FALSE, gamma=FALSE, l.start=23.56)

rainseriesforecasts2 <- forecast(rainseriesforecasts,h=8)

plot(rainseriesforecasts2)

acf(rainseriesforecasts2$residuals, lag.max=20 , na.action = na.pass)

Box.test(rainseriesforecasts2$residuals, lag=20, type="Ljung-Box")

plot.ts(rainseriesforecasts2$residuals)

plotForecastErrors <- function(forecasterrors)
{
  # make a histogram of the forecast errors:
  mybinsize <- IQR(forecasterrors)/4 
  mysd	<- sd(forecasterrors)
  mymin <- min(forecasterrors) - mysd*5 
  mymax <- max(forecasterrors) + mysd*3
  # generate normally distributed data with mean 0 and standard deviation mysd
  mynorm <- rnorm(10000, mean=0, sd=mysd)
  mymin2 <- min(mynorm)
  mymax2 <- max(mynorm)
  if (mymin2 < mymin) { mymin <- mymin2 }
  if (mymax2 > mymax) { mymax <- mymax2 }
  # make a red histogram of the forecast errors, with the normally distributed data overlaid:
  mybins <- seq(mymin, mymax, mybinsize)
  hist(forecasterrors, col="red", freq=FALSE, breaks=mybins)
  # freq=FALSE ensures the area under the histogram = 1
  # generate normally distributed data with mean 0 and standard deviation mysd
  myhist <- hist(mynorm, plot=FALSE, breaks=mybins)
  # plot the normal curve as a blue line on top of the histogram of forecast errors:
  points(myhist$mids, myhist$density, type="l", col="blue", lwd=2)
}



rainseriesforecasts2$residuals <-rainseriesforecasts2$residuals[!is.na(rainseriesforecasts2$residuals)]
plotForecastErrors(rainseriesforecasts2$residuals)



##### Holt-Winters Exponential Smoothing (additive model with increasing or decreasing trend and no seasonality) #####

skirts <- scan("http://robjhyndman.com/tsdldata/roberts/skirts.dat",skip=5) 
skirtsseries <- ts(skirts,start=c(1866))
plot.ts(skirtsseries)

skirtsseriesforecasts <- HoltWinters(skirtsseries, gamma=FALSE)
skirtsseriesforecasts 

plot(skirtsseriesforecasts)

HoltWinters(skirtsseries, gamma=FALSE, l.start=608, b.start=9)

skirtsseriesforecasts2 <- forecast(skirtsseriesforecasts, h=19)
plot(skirtsseriesforecasts2)

acf(skirtsseriesforecasts2$residuals, lag.max=20, na.action = na.pass)
Box.test(skirtsseriesforecasts2$residuals, lag=20, type="Ljung-Box")

plot.ts(skirtsseriesforecasts2$residuals) # make time series plot

skirtsseriesforecasts2$residuals <-skirtsseriesforecasts2$residuals[!is.na(skirtsseriesforecasts2$residuals)]

plotForecastErrors(skirtsseriesforecasts2$residuals) # make a histogram


##################### Holt-Winters Exponential Smoothing (additive model with increasing or decreasing trend and seasonality) #########

logsouvenirtimeseries <- log(souvenirtimeseries)
souvenirtimeseriesforecasts <- HoltWinters(logsouvenirtimeseries)
souvenirtimeseriesforecasts

plot(souvenirtimeseriesforecasts)

souvenirtimeseriesforecasts2 <- forecast(souvenirtimeseriesforecasts,h=48)
plot(souvenirtimeseriesforecasts2)

acf(souvenirtimeseriesforecasts2$residuals, lag.max=20 , na.action = na.pass)
Box.test(souvenirtimeseriesforecasts2$residuals, lag=20, type="Ljung-Box")


plot.ts(souvenirtimeseriesforecasts2$residuals) # make time series plot

souvenirtimeseriesforecasts2$residuals <- souvenirtimeseriesforecasts2$residuals [!is.na(souvenirtimeseriesforecasts2$residuals)]

plotForecastErrors(skirtsseriesforecasts2$residuals) # make a histogram


 ###################Arima #########################

#### Differencing a Time Series #####

skirtsseriesdiff1 <- diff(skirtsseries, differences=1)
plot.ts(skirtsseriesdiff1)


skirtsseriesdiff2 <- diff(skirtsseries, differences=2)
plot.ts(skirtsseriesdiff2)

####### Selecting a Candidate ARIMA Model ######

kingtimeseriesdiff1 <- diff(kingstimeseries, differences=1)
plot.ts(kingtimeseriesdiff1)

acf(kingtimeseriesdiff1, lag.max=20)
acf(kingtimeseriesdiff1, lag.max=20, plot=FALSE)

pacf(kingtimeseriesdiff1, lag.max=20)
pacf(kingtimeseriesdiff1, lag.max=20, plot=FALSE)
auto.arima(kings)
auto.arima(kingtimeseriesdiff1)


volcanodust <- scan("http://robjhyndman.com/tsdldata/annual/dvi.dat", skip=1) 
volcanodustseries <- ts(volcanodust,start=c(1500))
plot.ts(volcanodustseries)

acf(volcanodustseries, lag.max=20)
# plot a correlogram
acf(volcanodustseries, lag.max=20, plot=FALSE) # get the values

pacf(volcanodustseries, lag.max=20)
pacf(volcanodustseries, lag.max=20, plot=FALSE)

auto.arima(volcanodust)

auto.arima(volcanodust,ic='bic')


####### Forecasting Using an ARIMA Model #######

kingstimeseriesarima <- arima(kingstimeseries, order=c(0,1,1)) # fit an ARIMA(0,1,1) model
kingstimeseriesarima 

kingstimeseriesforecasts <- forecast(kingstimeseriesarima, h=5)
kingstimeseriesforecasts

plot(kingstimeseriesforecasts)

acf(kingstimeseriesforecasts$residuals, lag.max=20)
Box.test(kingstimeseriesforecasts$residuals, lag=20, type="Ljung-Box")

plot.ts(kingstimeseriesforecasts$residuals)
plotForecastErrors(kingstimeseriesforecasts$residuals) # make a histogram

volcanodustseriesarima <- arima(volcanodustseries, order=c(2,0,0))
volcanodustseriesarima 

volcanodustseriesforecasts <- forecast(volcanodustseriesarima, h=31)
volcanodustseriesforecasts

plot(volcanodustseriesforecasts)

acf(volcanodustseriesforecasts$residuals, lag.max=20)
Box.test(volcanodustseriesforecasts$residuals, lag=20, type="Ljung-Box")

plot.ts(volcanodustseriesforecasts$residuals)
plotForecastErrors(volcanodustseriesforecasts$residuals) 

