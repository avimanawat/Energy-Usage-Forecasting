setwd("C:/Users/avman/OneDrive/Desktop/Time Series/TSF excel files")
electric = read.csv("FINAL PROJECT.csv", header = TRUE)

install.packages("caret")

library("caret")
library(forecast)
library(nnet)
library(fGarch)
library(tseries)
library(fUnitRoots)
library(mda)


#dimention of data
head(electric)
tail(electric)
dim(electric)

#reading the hour wise data
electric.consumption = ts(electric[,10],start=c(2013,01,01), frequency=24*365.25)
plot(electric.consumption, col="red")

#building correlation matrix
cor(electric[,-1])

#######################################################################################

#Since the data is too big we will take small part and apply to whole series
setwd("C:/Users/avman/OneDrive/Desktop/Time Series/TSF excel files")
electric = read.csv("elec.csv", header = TRUE)

electric = ts(electric[,10],start=c(1), frequency=24)
plot(electric, col="blue")

#checking autocorrelation
acf(electric)
pacf(electric)

#convert to log form
log.data=log(electric)
acf(electric)
pacf(electric)

#Since there is trend taking difference of data
diff.data=diff(electric)
plot(diff.data)
acf(diff.data)
pacf(diff.data)

#dividing the data into test and training
length(electric)

stepah=500
training = length(electric)-stepah

train.ts = window(electric, end = c(1,training))
valid.ts = window(electric, start= c(1,training+1), end = c(1,training+stepah))

log.train.ts=window(log.data, end = c(1,training))

##Now we decompose the ts data to see trend and season in the time series
decsale1=decompose(electric, type="additive")
plot(decsale1)

decsale1=decompose(electric, type="multiplicative")
plot(decsale1)

#We can see some trend and proper seasonlity in data

?stl
?ts

sales.stl=stl(train.ts,t.window=13,s.window="periodic",robust=TRUE)
autoplot(sales.stl)

sales.stl2=stl(train.ts,t.window=13,s.window=24,robust=TRUE)
autoplot(sales.stl2)

sales.St = seasonal(sales.stl2)
plot(sales.St)

pred = forecast(sales.stl,h=500)
pred2 = forecast(sales.stl2,h=500)

accuracy(pred)
accuracy(pred2)

electric.pred.naive = stlf(sales.St, method='naive',h=12)
accuracy(electric.pred.naive)
?stlf


electric.pred.ets = stlf(electric, method="ets",h=12)
accuracy(electric.pred.ets)

electric.pred.rwdrift = stlf(electric,method="rwdrift",h=12)
accuracy(electric.pred.rwdrift)




sales.stl.pred <- forecast(sales.stl,h=12)

accuracy(sales.stl.pred) #RMSE 84.26



### we will use exponential smoothing method



plot(electric.pred.naive)

autoplot(electric.pred.naive)




################################### Linear Model ######################################

plot(train.ts)

tslm.model = tslm(train.ts ~ trend +season) #fitting quadradic model
summary(tslm.model)

pred10=forecast(tslm.model, h=500)

accuracy(pred10, valid.ts)


################################## Exponential smoothing ###################################

es.model <- ets(train.ts)
es.fcast <- forecast(es.model,h=500)
plot(es.fcast)
accuracy(es.fcast,valid.ts)

################################## Naive Method #######################################
naive.mtd=naive(train.ts,h=500)
accuracy(naive.mtd,valid.ts)

  ################################ Smoothing Method #########################################

smoothing.mtd=ses(train.ts,h=500)
accuracy(smoothing.mtd,valid.ts)


################################### Arima Models #######################################

arima.model=auto.arima(log.data)
summary(arima.model)

tsdisplay(arima.model$residuals)
Box.test(arima.model$residuals) #no autocorrelation

acf(arima.model$residuals^2)
Box.test(arima.model$residuals^2) #auto-correlation


#there is GARCH effect in the model

garch.model=garchFit(~arma(3,0)+garch(1,1),data=train.ts, trace = FALSE)
summary(garch.model)

sresi=residuals(garch.model,standardize=T) # Obtain standardized residuals
sigma.t=volatility(garch.model)  # obtain the fitted volatility sigma_t.

Box.test(sresi,10,type='Ljung')
Box.test(sresi^2,10,type='Ljung')

qqnorm(sresi)
qqline(sresi)

predict(sigma.t,500)


############################# Try Other Arima Model ######################

RMSE.matrix <- matrix(data = NA, nrow=500, ncol=7)

count=0
for(a in 0:2){
  for(b in 0:1){
    for(c in 0:2){   
      for(d in 0:2){
        for(e in 0:1){
          for(f in 0:2){
            tryCatch({
              fit <- arima(log.train.ts, order=c(a,e,b),seasonal= list(order=c(c,f,d),period=24))
              acc<- accuracy(fit)
              RMSE.matrix[count+1,1] <- acc[[2]]#provide RMSE
              RMSE.matrix[count+1,2]<- a
              RMSE.matrix[count+1,3]<- e
              RMSE.matrix[count+1,4]<- b
              RMSE.matrix[count+1,5]<- c   
              RMSE.matrix[count+1,6]<- f  
              RMSE.matrix[count+1,7]<- d  
              count= count+1
            },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
          }
        }
      }
    }
  }
}

colnames(RMSE.matrix)= c("RMSE","AR","MA","Season AR","Season AR","Season AR","Season MA")

RMSE.matrix 
RMSE.matrix[which.min(RMSE.matrix[,1]),] #find lowest RMSE
df.RMSE <- as.data.frame(RMSE.matrix) #make data frame for sorting
df.RMSE


###to find both aic and rmse value and compare together
x=list()
z=list()
no=list()


for (p in 0:2){
  for (d in 0:1){
    for (q in 0:2){
      for (P in 0:2){
        for (D in 0:2){
          for (Q in 0:2){
            tryCatch({
              m=arima(log.train.ts,order=c(p,d,q),seasonal=list(order=c(P,D,Q),period=24))
              aicvalue=m$aic
              x=append(x,aicvalue)
              a=accuracy(m)
              rmsevalue=a[2]
              z=append(z,rmsevalue)
              b=list(c(p,d,q,P,D,Q))
              no=append(no,b)
            },error=function(e){})
          }}}}}}



x
z
no

#best model arima(1,1,0)(0,1,2)[24]

l=rbind(x,z,no)
arima.model.1=arima(log.train.ts, order = c(1,1,0),seasonal = list(order=c(0,1,2), period=24))
summary(arima.model.1)

pred.arima.log=forecast(arima.model.1,h=500)
plot(pred.arima.log)
pred.arima.exp=exp(pred.arima.log$mean) #data back to normal form


####################### Dynamic Regression ######################################

#try to fit other city data into one model

data1=electric
data.csv=read.csv("elec.csv", header = TRUE)
data.csv.ts=ts(data.csv, start = c(1), frequency = 24)
data2=data.csv.ts[,2]
data3=data.csv.ts[,3]
data4=data.csv.ts[,4]
data5=data.csv.ts[,5]
data6=data.csv.ts[,6]
data7=data.csv.ts[,7]
data8=data.csv.ts[,8]
data9=data.csv.ts[,9]

#try to fit various graphs
plot(data1,ylim=c(4000,20000), col="red")
lines(data6)
lines(data8,col="green")

#see correlationship
plot(data1,data6,lty=15,pch=20)

plot(diff(data1), diff(data6), lty=15, pch=20)


dynamic.model.1=lm(data1~data6)
summary(dynamic.model.1)

acf(dynamic.model.1$residuals) #violation
Box.test(dynamic.model.1$residuals,lag=10,type="Ljung")

ar(data1)
adfTest(data1,lag=30)#test

#various dynamic model
dynamic.model.1 = auto.arima(data1, xreg = data3)#model looking good but AIC score is very high
summary(dynamic.model.1)
tsdiag(dynamic.model.1)

#forecasting the model
dynamic.forecast.1 <- forecast(dynamic.model.1,xreg=rep(mean(data3),500), h=500)

autoplot((dynamic.forecast.1),xlab="Time" ,ylab="consumption")
?autoplot


############################## Random Forest ##############################################

electric.binary = read.csv("elec.csv", header = TRUE)

mean(electric.binary$PJMW)#takin mean to see the demand of future will increase of decrease

electric.binary$Date <- as.Date(electric.binary$Datetime, format="%m/%d/%Y")
tail(electric.binary$Date)

electric.binary$increase.demand <- ifelse(electric.binary$PJMW > 5848.158, 1, 0)

nPeriods <- length(electric.binary$increase.demand)
electric.binary$Lag1 <- c(NA,electric.binary$PJMW[1:(nPeriods-1)])
electric.binary$t <- seq(1, nPeriods, 1)

electric.binary$Seasonal_sine = sin(2 * pi * electric.binary$t / 365.25)
electric.binary$Seasonal_cosine = cos(2 * pi * electric.binary$t / 365.25)

train.binary <- electric.binary[electric.binary$Date <= as.Date("07/13/2018", format="%m/%d/%Y"), ]
train.binary1 <- train.binary[-1,]

valid.binary <- electric.binary[electric.binary$Date > as.Date("07/13/2018", format="%m/%d/%Y"), ]
valid.binary2 <- valid.binary[, c(13,15)]

binary.linear.model <- glm(increase.demand ~ Lag1+Seasonal_sine, data = train.binary, family = "binomial")
summary(binary.linear.model)

binary.linear.model.pred <- predict(binary.linear.model, valid.binary2, type = "response") 

accuracy(binary.linear.model)

confusion(ifelse(binary.linear.model$fitted > 0.5, 1, 0), train.binary1$increase.demand, positive="1")

confusion(ifelse(binary.linear.model.pred > 0.5, 1, 0), valid.binary$increase.demand, positive="1")


##################################### try some neural network and more advance model ##############################################

?nnetar

nn.model <- nnetar(train.ts, p=23, P = 2, size =7)
nn.fcast <- forecast(nn.model,h=500)
plot(nn.fcast)
accuracy(nn.fcast,valid.ts)
accuracy(nn.model)



################################## TBAT model ###############################################################

tbat.model=tbats(train.ts)
summary(tbat.model)
tbats.fcasr=forecast(tbat.model,h=500)
plot(tbats.fcasr)
accuracy(tbats.fcasr,valid.ts)
