  rm(list=objects())
  
  #2545
library(tidyverse)
library(lubridate)
library(ranger)

#2009-06-01
#2017-09-30
#2545 correspond à 2016-06-01
#2666 correspond à 2016-09-30
##on choisit de considérer 7 années complètes consécutives pour l'entrainement, puis 1 année consécutive complète
###le train va donc du 2009-06-01 au  2016-06-01 (7*365 jours) et le test va du 2016-06-02	au 2017-06-03	 (365 jours)

#on va creer artificiellement un test de 1 an et on ignore donc les derniers jours des données (119 jours ignorés)

library(xts)
library(dygraphs)
library(forecast)
train <- read.csv(file="montsouris.csv", sep=",", dec='.')
Date = as.POSIXct(strptime(train$DATE, "%Y-%m-%d ",tz="GMT"))

n0=2545 # kes 7 années d'entrainement durent 2545 jours
n1=365  #l'annee test dure 365 jours

Date_train=Date[1:2545]
n0=length(Date_train)
Date_test=Date[2546:2910]
Date_total=Date[1:(n0+n1)]
n1=length(Date_test)
length(Date_total)-n0-n1



which(is.na(Date_train))

TAVG_train=train$TAVG[1:n0]
TAVG_test=train$TAVG[(n0+1):(n0+n1)]
TAVG_total=train$TAVG[1:(n0+n1)]
length(TAVG_test)
TAVG_xts=xts(TAVG_train,order.by = Date_train)
NROW(TAVG_train)
TAVG_ts=ts(TAVG_train,frequency=365)

tab_x=c(1:length(Date_train))


plot(TAVG_train)
plot(TAVG_test)

TAVG_test_xts<-xts(TAVG_test,order.by=Date_test)
TAVG_train_xts<-xts(TAVG_train,order.by=Date_train)

TAVG_total_xts<-xts(TAVG_total,order.by=Date_total)

graphe <- cbind( TAVG_test_xts,TAVG_train_xts)
#dygraph(graphe,main=)
dygraph(graphe,main="Relevés journaliers de la température moyenne ,Parc Montsouris",xlab="Date",ylab="Temperaure Moyenne (en °F)")

###spectres non exploitables.
spectrum(TAVG_test)
spectrum(TAVG_train)


###utilisation d'un filtrage par moyenne mobile.On choisit d'utiliser une fenêtre de largeur 365, car 365 est la periode de notre série.
##Ainsi,on filtre 
TAVG_detrend=stats::filter(TAVG_train, filter=array(1/365,dim=365), method = c("convolution"),
           sides = 2, circular = FALSE)



plot(Date_train, TAVG_train, type = "l")
lines(Date_train,TAVG_detrend,col="red")   


TAVG_detrend<-xts(TAVG_detrend,order.by=Date_train)
plot(TAVG_train)
lines(TAVG_detrend)
which(is.na(TAVG_detrend))
dygraph(TAVG_detrend)

pacf(TAVG_train)
dygraph(TAVG_detrend, 
              main = "Tendance obtenue par filtrage par moyenne mobile ", 
              xlab = "Année", 
              ylab = "Température (°F )")


par(mfrow = c(2, 2))

#10,25,100,1000
####################regression à noyau

noyau <- ksmooth(tab_x, TAVG_train, kernel = c("normal"),range.x=range(tab_x),n.points=max(length(tab_x),length(TAVG_train)),
                 bandwidth =500)
                   
plot(Date_train, TAVG_train, type = "l", xlab = "",
     ylab = "Température observée", col = "black")


#####################

legend=c("Temperature observee","Lissage par noyau gaussien",col=c("red","blue"))







  lines(Date_train, noyau$y, col = "red", lwd = 2)
  
  
  legend(1, 95, legend=c("Line 1", "Line 2"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
plot(Date_train, TAVG_train - noyau$y, type = "l",
     xlab = "", ylab = "Residus", col = "black")




acf(TAVG_train - noyau$y)


pacf(TAVG_train - noyau$y)



par(mfrow = c(1, 1))



trend=TAVG_train-TAVG_detrend
plot(trend)




plot(tab_x, TAVG_train, type="b", pch=19, col="black", xlab="x", ylab="y")
# Ajouter une ligne
lines(tab_x, noyau$y, pch=18, col="red", type="b", lty=2)
# Ajouter une légende
legend(1, 20, legend=c( "Regression sur noyau gaussien"),
       col=c("red", "blue"), lty=2:2, cex=0.9)


              


prev=arima.sim(list(fit_arma$coef),n=100)
plot(prev)

predict(auto.arima)
mean(trend,na.rm=TRUE)
plot(trend)
plot(TAVG_train,type='l')
spectrum(TAVG_train)
lines(TAVG_detrend,col='red')
plot(TAVG_train)
###juste pour vérifier nos resultats).
TAVG_dec=decompose(TAVG_ts,type="additive")
plot(TAVG_dec)



X=TAVG_train
t=c(1:length(Date_train))
###################################################estimation de la tendance
#################regression
#reg<-lm(X~t+I(t^2))
reg<-lm(X~t)
summary(reg)
ychap.lm <- reg$fitted
# ychap.lm<-xts(as.numeric(reg$fitted),order.by=Date)
# ychap.lm<-xts(reg$fitted,order.by=Date)
plot(X,type='l')
lines(ychap.lm,col='red')
plot(TAVG_dec$trend)

reg0<-lm(TAVG_detrend~Date_train)

summary(reg0)




#################polynomes locaux
lo<-loess(X~t, degree=1, span=0.9)
ychap.lo<-xts(lo$fitted,order.by=Date_train)
plot(X,type='l')
lines(ychap.lo,col='red')
?loess

plot(ychap.lo)



which(is.na())
TAVG_detrend[1:4]=59.9
TAVG_detrend[2662:2666]=62.7

TAV
spectrum(TAVG_detrend*TAVG_detrend)
spectrum(TAVG_train)

which(is.na(TAVG_dec$trend))


###################################################estimation de la partie saisonniere
#################regression

t=c(1:length(Date_train))
w=2*pi/365
fourier<-cbind(cos(w*t), sin(w*t))
K<-5
for(i in c(2:K))
{
  fourier<-cbind(fourier,cos(i*w*t), sin(i*w*t))
}
matplot(fourier,type='l')
dim(fourier)

#eg<-lm(TAVG_train~fourier[,1:2])
#reg<-lm(TAVG_train~fourier[,1:2])
##50
#reg<-lm(TAVG_train~fourier[,1:2])

###pour aller plus loin
reg<-lm(noyau$y~fourier[,1:2])

###
reg_th_fourier=lm(TAVG_dec$seasonal~fourier[,1:2])

summary(reg_th_fourier)
reg_th_fourier$coefficients

plot(reg_th_fourier$fitted.values)




########################
par(mfrow = c(1, 1))

TAVG_s=TAVG_train-reg$fitted.values
TAVG_eps=TAVG_train-reg_th_fourier$fitted.values
plot(TAVG_eps)
mean(TAVG_eps)


TAVG_essai=TAVG_train-noyau$y
acf(TAVG_essai)

plot(TAVG_eps)
mean(TAVG_s)


par(mfrow=c(1,1))
a=0.80
acf(TAVG_eps,lag.max=30)
lines(c(0:30),a^c(0:30),col='red')



pacf(TAVG_eps)



################################################

pacf(TAVG_s_th)
auto.arima(TAVG_s,trace=TRUE,ic="aic")


modele_arma=auto.arima(TAVG_eps,trace=TRUE,ic="aic")
modele_arma$coef
plot(modele_arma$fitted)
par(mfrow=c(1,1))


####simulation d'un ARMA
####simulation=arima.sim(list(ar=modele_arma$coef["ar1"],ma=modele_arma$coef["ma1"]),sd=modele_arma$sigma,2545) +reg_th_fourier$fitted.values
####lines(Date_train,TAVG_train,col="red")
#####

sim=arima.sim(list(ar=modele_arma$coef["ar1"],ma=modele_arma$coef["ma1"]),sd=modele_arma$sigma,2545)
plot(sim)


length(Date_total)
n0+n1



plot(modele_arma$residuals)
acf(modele_arma$residuals)

######

library(forecast)


par(mfrow=c(1,1))



coefs=fourier(TAVG_dec$seasonal,K=1)

coefs_train=fourier(TAVG_dec$seasonal,K=1)
coefs_test=fourier(TAVG_dec$seasonal,K=1,h=365)
fit <- Arima(TAVG_train, order=c(1,0,1), xreg=coefs_train)

prevision_arma=forecast(fit,xreg=coefs_test,level=c(80,95,99))

plot(forecast(fit, h=2*m, xreg=coefs_test))
lines(Date_total[1:n0],TAVG_train,type='l')
length(TAVG_test)


n1
n0
n0+n1
length(prevision_arma$mean)
length(prevision_arma$mean)
plot(prevision_arma$mean)
plot(prevision_arma$upper)


#######################simulation

par(mfrow=c(1,1))
n0=length(TAVG_train)
length(Date_total)
length(train$TAVG)
plot(Date_total,TAVG_total,type='l')
length(Date_total[(n0+1):length(Date_total)])
lines(Date_test,TAVG_test,col='green')


lines(Date_total[1:n0],TAVG_train,col='black')
lines(Date_total[(n0+1):(n0+n1)],TAVG_test,col='green')
lines(Date_total[(n0+1):(n0+n1)],prevision_arma$mean ,col='red')
lines(Date_train,fit$fitted,col="red")
lines(Date_total[(n0+1):(n0+n1)],prevision_arma$upper ,col='red')
length(prevision_arma$upper)

length(prevision_arma$mean)
plot(prevision_arma$upper)
head(prevision_arma$mean)
head(prevision_arma$lower)
plot(prevision_arma$lower[368])


length(prevision_arma$upper[,1])
p1_xts=xts(prevision_arma$upper[,1],order.by=Date_test)
p2_xts=xts(prevision_arma$upper[,2],order.by=Date_test)
p3_xts=xts(prevision_arma$upper[,3],order.by=Date_test)
pm_xts=xts(prevision_arma$mean,order.by=Date_test)
p4_xts=xts(prevision_arma$lower[,1],order.by=Date_test)
p5_xts=xts(prevision_arma$lower[,2],order.by=Date_test)
p6_xts=xts(prevision_arma$lower[,3],order.by=Date_test)




graphe_prev <- cbind( p1_xts,p2_xts,p2_xts,p3_xts,p4_xts,p5_xts,p6_xts,pm_xts)


u=dygraph(graphe_prev,main="Prévisions pour la température moyenne journalière ,Parc Montsouris",ylab="Temperaure Moyenne (en °F)",xlab="Date")


###verification des residus
checkresiduals(modele_arma$residuals,plot=TRUE)
acf(modele_arma$residuals)




#######
#####modele obtenu=modele ARIMA(1,0,1)=ARMA(1,1) whith non zero mean


#####################???
Box.test(modele_arma$residuals, lag =5, type = "Box-Pierce", fitdf = 2)


Box.test(modele_arma$fitted, lag = 1, type = c("Box-Pierce"),fitdf=2)



####le test de Box-Pierce a été concluant

x <- rnorm (100)
Box.test (x, lag = 1)
Box.test (x, lag = 1, type = "Ljung")
Box.test (x, lag = 1, type = "Box-Pierce")





###utilisation des fonctions de lissages
source(file="smoothing.R")

TAVG.smooth<-expSmooth(TAVG_train,0.5)
plot(TAVG_train,type='l')
lines(TAVG.smooth,col='red')

plot(tail(TAVG_train,n-1),type='l')

lines(head(TAVG.smooth,n-1),col='red')

tail(TAVG_train,n-1)
mean((tail(TAVG_train,n-1)-head(TAVG.smooth,n-1))^2)

#=17.94


###########
alpha<-seq(0.05,0.95,length=100)
forecast<-lapply(alpha,expSmooth,x=TAVG_train)
erreur_simple<-unlist(
  lapply(forecast,
         function(x){mean((tail(TAVG_train,n-1)-head(x,n-1))^2)}))
plot(alpha,erreur_simple,type='l')

TAVG.smooth<-expSmooth(TAVG_train,alpha[which.min(erreur_simple)])
plot(TAVG_train,type='l')
lines(TAVG.smooth,col='red')
mean((tail(TAVG_train,n-1)-head(TAVG.smooth,n-1))^2)
#mean=14.6

####utilisation du lissage double
alpha<-seq(0.05,0.95,length=100)
forecast<-lapply(alpha,DoubleExpSmooth,x=TAVG_train)
erreur_double<-unlist(
  lapply(forecast,
         function(x){mean((tail(TAVG_train,n-1)-head(x$smooth,n-1))^2)}))
plot(alpha,erreur_simple,type='l')

TAVG_train.smooth<-DoubleExpSmooth(TAVG_train,alpha[which.min(erreur_double)])
plot(TAVG_train,type='l')
lines(TAVG_train.smooth$smooth,col='red')

plot(TAVG_train.smooth$l,type='l',ylim=range(TAVG_train.smooth$l,TAVG_train.smooth$b),col='blue')
lines(TAVG_train.smooth$b,col='purple')



n0<-length(Date_train)
n1<-length(Date_test)
n0
n1
length(TAVG_test)
length(Date_total[(n0+1):length(Date_total)])
n0
 

par(mfrow=c(1,1))

plot(Date_total,train$TAVG,type='l')
lines(Date_total[1:n0],TAVG_train,type='l')
length(Date_total[(n0+1):length(Date_total)])
lines(Date_total[(n0+1):length(Date_total)],TAVG_test,col='green')
#predict.expSmooth<-function(Xsmooth,inst,horizon,smooth.type="double")
TAVG_train.smooth.simple<-expSmooth(TAVG_train,alpha=alpha[which.min(erreur_simple)])


TAVG_train.smooth.double<-DoubleExpSmooth(TAVG_train,alpha[which.min(erreur_double)])

TAVG.smooth.simple.forecast<-predict.expSmooth(TAVG_train.smooth.simple,n0,n1,smooth.type="simple")
TAVG.smooth.double.forecast<-predict.expSmooth(TAVG_train.smooth.double,n0,n1,smooth.type="double")

#plot(TAVG.smooth.simple.forecast,col='blue')
#plot(TAVG.smooth.double.forecast,col='purple')
lines(Date_total[1:length(Date_total)],TAVG.smooth.simple.forecast,col='red')
lines(Date_total[1:length(Date_total)],TAVG.smooth.double.forecast,col='blue')




###approche holt-winters

ht=HoltWinters(TAVG_ts, alpha = NULL, beta = NULL, gamma = NULL,
            seasonal = c("additive", "multiplicative"),
            start.periods = 2, l.start = NULL, b.start = NULL,
            s.start = NULL,
            optim.start = c(alpha = 0.3, beta = 0.1, gamma = 0.1),
            optim.control = list())

length(ht$x)
#plot(ht$x)
#lines(Date_total[1:length(Date_train)],ht$x)

###prevision 
ht_f=predict(ht,n.ahead=n1
,prediction.interval = FALSE,level=0.95)
length(ht_f)
plot(ht_f)

rmse=function(obs,pre)
{
  return(sqrt((mean(obs-pre))^2))
}

rmse(TAVG_ts,ht_f)

#lines(Date_total[(n0+1):(n0+n1)],ht_f,col="red")


