# TIMES SERIES ANALYSIS ASSiGNMENT
# GERARD HOLIAN
# 16170571

install.packages("TSA")
install.packages("tseries")
install.packages("forecast")

library(TSA)
library(tseries)
library(forecast)



kts1<-ts(kayak08.18, freq=12, start=c(2008,1))


# test for random walk or white noise

dev.new(width=14, height= 10)

plot(kts1, main="Kayaking searches worldwide 01/2008-12/2017", ylab="Number of searches", type="l")


acf1<-acf(kts1)


# shows there is autcorrelations and can rule out white noise

acf2<-acf(diff(kts1))


# shows there is autocorrelations and can rule out random walk

# Data is suitable for analysis




# Need to remove 10% of observations 
length(kts1)

# Length = 120, Therefore need to remove 12 observations

newkts1<- head(kts1, -12)

length(newkts1)


# 108 observations now

ktstrain<-ts(newkts1, freq=12, start=c(2008,1))

# plot of new time series

plot(ktstrain, main="Kayaking searches worldwide 01/2008-12/2016", ylab= "Number of Searches", xlab= "Time", type="l")


# looks like overall and seasonal trend



# test if series is stationary

# Shapiro wilk test used to test for normal distribution

shapiro.test(ktstrain)

# data not normally distributed therefore do boxcox test

hist(ktstrain)
# boxcox test determines how to transform

bck1<-BoxCox.ar(ktstrain, lambda=seq(-2, 2, 0.1))

bck1$mle

bck1$ci

# no transformation suggested



adf.test(ktstrain)


# Decompose series

decompk1<-decompose(ktstrain)



plot(decompk1)
# again looks like seasonal and overall trend

diff(range(ktstrain))

diff(range(decompk1$seasonal, na.rm=T))

# seasonality trend of 56.11 units

diff(range(decompk1$random, na.rm=T))

# noise of 16.97 units 

diff(range(decompk1$trend, na.rm=T))

# overall trend of 12.38 units




# appears to require differencing and seasonal differencing
ktstrain2 <- diff(diff(ktstrain,lag=1),lag=12)

plot(ktstrain2)

# looks stationary

hist(ktstrain2)

adf.test(ktstrain2)
# no more differencing required

shapiro.test(ktstrain2) #normally distributed






# decompose again

decompk2<-decompose(ktstrain2)

plot(decompk2)

diff(range(ktstrain2))

diff(range(decompk2$seasonal, na.rm=T))

# seasonality trend reduced to 3.9

diff(range(decompk2$random, na.rm=T))

# noise of 22.33

diff(range(decompk2$trend, na.rm=T))

# overall trend reduced to 2.29









# models

acf(ktstrain2)

#This suggests that not only is there seasonal autocorrelation (correlation between 
#values 12 months apart) but also short term autocorrelation (correlation between neighbouring values, i.e., one month apart)


pacf(ktstrain2)

eacf(ktstrain2)

ktarma<-armasubsets(ktstrain2, nar=14, nma = 14)
plot(ktarma)



# models

model1 <- arima(ktstrain, order = c(1,1,4), seasonal=list(order=c(1,1,0), period=12))
summary(model1)
# aic 499.11
model2 <- arima(ktstrain, order = c(1,1,10), seasonal=list(order=c(1,1,0), period=12))
summary(model2)
# aic = 504.29
model3<- arima(ktstrain, order = c(1,1,6), seasonal=list(order=c(1,1,0), period=12))
summary(model3)
# aic =501.81
model4<- arima(ktstrain, order = c(0,1,3), seasonal=list(order=c(1,1,0), period=12))
summary(model4)
# aic=498.42
model5<- arima(ktstrain, order = c(0,1,0), seasonal=list(order=c(1,1,0), period=12))
summary(model5)
# aic=514.58


#pick model4 as lowest AIC and simplest model


# overfitting models ARIMA (0,1,3) X (1,1,0)
summary(model4)
# aic=498.42

model4.a<- arima(ktstrain, order = c(1,1,3), seasonal=list(order=c(1,1,0), period=12))
summary(model4.a)
# aic = 501.6
model4.b<- arima(ktstrain, order = c(2,1,3), seasonal=list(order=c(1,1,0), period=12))
summary(model4.b)
# aic = 500.47
# Leave out AR term


model4.c<- arima(ktstrain, order = c(0,1,2), seasonal=list(order=c(1,1,0), period=12))
summary(model4.c)
# aic = 498.59
model4.e<- arima(ktstrain, order = c(0,1,4), seasonal=list(order=c(1,1,0), period=12))
summary(model4.e)
# aic = 497.26
model4.d<- arima(ktstrain, order = c(0,1,1), seasonal=list(order=c(1,1,0), period=12))
summary(model4.d)
# aic = 496.8
# use MA(1) term


model4.f<- arima(ktstrain, order = c(0,1,1), seasonal=list(order=c(1,1,1), period=12))
summary(model4.f)
# aic = 498.64
model4.g<- arima(ktstrain, order = c(0,1,1), seasonal=list(order=c(1,1,2), period=12))
summary(model4.g)
# aic = 500.41
# leave out seasonal MA term



model4.h<- arima(ktstrain, order = c(0,1,1), seasonal=list(order=c(0,1,0), period=12))
summary(model4.h)
# aic = 510.46
model4.i<- arima(ktstrain, order = c(0,1,1), seasonal=list(order=c(2,1,0), period=12))
summary(model4.i)
# aic = 498.61
# leave seasonal AR term as 1



# final models is model4 and model4.d




#residuals

res4d<-rstandard(model4.d)
res4e<-rstandard(model4.e)

lbres3<-tsdiag(model4.d)
lbres4<-tsdiag(model4.e) #4e has better lung box p-values


plot(res4e)



qqnorm(res4d); qqline(res4d)
qqnorm(res4e); qqline(res4e)







hist(res4d)
hist(res4e)


shapiro.test(res4e) # not normal
shapiro.test(res4d) #normal



plot(res4d)
plot(res4e)




LB.test(model4.e, lag=12) # white noise 
LB.test(model4.d, lag=12) # not white noise
 # white noise

# choose 4e

pred<-predict(model4.e,n.ahead=12)





plot(kts1, type='o',col='blue')
lines(pred$pred, col='red')
lines(pred$pred-pred$se, col='green')
lines(pred$pred+pred$se, col='green')

#predictions and confidence intervals
pred$pred
pred$pred-pred$se
pred$pred+pred$se
