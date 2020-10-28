
is=read.csv("D:/01.First Drive/Data/data/issac.csv",header=T)

timee=is[,1]
gg=is[,2]
tmtt=is[,3]
sstg=is[,4]
vvent=is[,5]
treatment1=is[,10]
length(tmtt)
length(vvent)
length(treatment1)
age=is[,8]

library(splines)
library(survival)



***testing which parametric model to use for survival time with no covariate

weib=survreg(Surv(timee,vvent)~1,dist="w")
expp=survreg(Surv(timee,vvent)~1,dist="exponential")
loglog=survreg(Surv(timee,vvent)~1,dist="logl")


****** AIC for no each model with no predictor*
extractAIC(weib)[2]
extractAIC(expp)[2]
extractAIC(loglog)[2]




**** Checking Weibull assumption for no predictor 
kp=survfit(Surv(timee,vvent)~1)
summ=summary(kp)
tim=summ$time
suv=summ$surv 
risk=summ$n.risk 
lgtym=log(tim) 
lgsv=log(-log(suv))
plot(lgtym,lgsv,xlab="log(t)",ylab=expression(log(-log(hat(S)*"(t)"))))
abline(lm(lgsv ~ lgtym))


************  weibull survivor curve for no predictor 
plot(tim,1-pweibull(tim,1/2.1,exp(11)),ylim=c(0.85,1),xlab="time",
ylab="survival probability",type="o")


********* weibull hazard curve for no predictor

plot(tim,pweibull(tim,1/2.1,exp(11)),ylim=c(0.02,0.3),xlab="time",
ylab="survivor probability",type="l")





**** ploting the weibul survivor curve and kaplan meier for no covariate 
kp=survfit(Surv(timee,vvent)~1)
plot(tim,1-pweibull(tim,1/2.1,exp(11)),ylim=c(0.86,1),xlab="time"
,ylab="survivor probability",type="l",col="red") 
lines(kp,conf.int=FALSE)
legend("bottomleft",legend=c("black = kaplan Meir","red=Weibull"),col=1:2)




********* Estimating weibull for no predictor 
w=survreg(Surv(timee,vvent) ~1,dist="weibull",scale=0)
summary(w)


**********   AIC selection for full model

wete=survreg(Surv(timee,vvent)~gg+as.factor(tmtt)+as.factor(sstg),dist="w")

expte = survreg(Surv(timee,vvent) ~ gg+as.factor(tmtt)+as.factor(sstg),
 dist="exponential") 

 logte=survreg(Surv(timee,vvent)~gg+as.factor(tmtt)+as.factor(sstg),dist="logl")
 extractAIC(wete)[2]
 extractAIC(expte)[2]
 extractAIC(logte)[2]


********* checking weibull assumption for cancer stage 
stage=survfit(Surv(timee,vvent) ~as.factor(sstg))
sum=summary(stage) 
suv3=survfit(Surv(timee[sstg==3],vvent[sstg==3])~ 1)
suv4=survfit(Surv(timee[sstg==4],vvent[sstg==4]) ~ 1)
sum3=summary(suv3) 
sum4=summary(suv4)
suvstg3=sum3$surv 
time3=sum3$time 
suvstg4=sum4$surv 
time4=sum4$time
ti=sum$time 
tim3log=log(time3)
suvstg3log=log(-log(suvstg3))
tim4log=log(time4) 
suvtg4log=log(-log(suvstg4)) 
par(mfrow=c(2,1)) 
plot(tim3log,suvstg3log,xlab="log(t)",ylab="log(-logS(t))")
abline(lm(suvstg3log ~ tim3log))
plot(tim4log,suvtg4log,xlab="log(t)",ylab="log(-logS(t))")
abline(lm(suvtg4log ~tim4log)) 





*******  ploting the weibull survival for cancer stages
par(mfrow=c(1,1))
plot(tim,1-pweibull(tim,1/2.04,exp(47.062)),xlab="time(days)",ylab="survivor probability",ylim=c(0.84,1),type="l")
lines(tim,1-pweibull(tim,1/2.04,exp(12.638)),col=2)
lines(tim,1-pweibull(tim,1/2.04,exp(11.448)),col=3)
lines(tim,1-pweibull(tim,1/2.04,exp(10.757)),col=4)
legend("bottomleft",legend=c("black = stage 1","red = stage 2",
 "green = stage 3", "blue = stage 4"),col=1:2)


*********  ploting the weibull hazard for cancer stages

plot(tim,pweibull(tim,1/2.04,exp(12.6)),ylim=c(0.01,0.2),xlab="time(days)",ylab="hazard probability",type="b") 
lines(tim,pweibull(tim,1/2.04,exp(11.)),type="l",col=2)
lines(tim,pweibull(tim,1/2.04,exp(10.4)),type="o",col=4)
legend("topleft",legend=c("black = stage 2","red = stage 3", "green =
 stage 3"),col=1:3)


*********fitting Weibull model************

weee=survreg(Surv(timee,vvent)~as.factor(sstg),dist="w") 
summary(weee)
