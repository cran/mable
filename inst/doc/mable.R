## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo = FALSE-------------------------------------------------------------
library(mable)

## -----------------------------------------------------------------------------
data(Vaal.Flow)
head(Vaal.Flow, 3)

## ----results = "hide", warning = FALSE, message=FALSE-------------------------
vaal<-mable(Vaal.Flow$Flow, M = c(2,100), interval = c(0, 3000), IC = "all",
       controls = mable.ctrl(sig.level = 1e-8, maxit = 2000, eps = 1.0e-9))

## ----fig.align='center', fig.cap="Vaal River Annual Flow Data\\label{fig:vaal-river-data-plot}", fig.width=7, fig.height=7----
op <- par(mfrow = c(1,2))
layout(rbind(c(1, 2), c(3, 3)))
plot(vaal, which = "likelihood", cex = .5)
plot(vaal, which = "change-point", lgd.x = "topright")
hist(Vaal.Flow$Flow, prob = TRUE, xlim = c(0,3000), ylim =c(0,.0022), breaks =100*(0:30), 
    main = "Histogram and Densities of the Annual Flow of Vaal River",
    border = "dark grey",lwd = 1, xlab = "Flow", ylab = "Density", col ="light grey")
lines(density(x<-Vaal.Flow$Flow, bw = "nrd0", adjust = 1), lty = 2, col = 2,lwd = 2)
lines(y<-seq(0, 3000, length=100), dlnorm(y, mean(log(x)), sqrt(var(log(x)))), 
    lty = 4, col = 4, lwd = 2)
plot(vaal, which = "density", add = TRUE, lwd = 2)
legend("topright", lty = c(1, 4, 2), col = c(1, 4, 2), bty = "n",lwd = 2, 
c(expression(paste("MABLE: ",hat(f)[B])), expression(paste("Log-Normal: ",hat(f)[P])),
    expression(paste("KDE: ",hat(f)[K]))))
par(op)

## ----fig.align='center', fig.cap="AIC and BIC based on Vaal River Data\\label{fig:Vaal-River-data-AIC-BIC-plot}", fig.width=6, fig.height=4----
M <- vaal$M[1]:vaal$M[2]
aic <- vaal$ic$AIC
bic <- vaal$ic$BIC
qhc <- vaal$ic$QHC
vaal.gcp <- optim.gcp(vaal) # choose m by gamma change-point model
lr <- vaal.gcp$lr
plot(M, aic, cex = 0.7, col = 1, xlab = "m", ylab = "", ylim = c(ymin<-min(aic,bic,qhc,lr),
      ymax<-max(aic,bic,qhc,lr)), main = "AIC, BIC, QHC, and LR")
points(M, bic, pch = 2, cex = 0.7, col = 2)
points(M, qhc, pch = 3, cex = 0.7, col = 3)
points(M[-1], lr, pch = 4, cex = 0.7, col = 4)
segments(which.max(aic)+M[1]-1->m1, ymin, m1, max(aic), lty = 2)
segments(which.max(bic)+M[1]-1->m2, ymin, m2, max(bic), lty = 2, col = 2)
segments(which.max(qhc)+M[1]-1->m3, ymin, m3, max(qhc), lty = 2, col = 3)
segments(which.max(lr)+M[1]->m4, ymin, m4, max(lr), lty = 2, col = 4)
axis(1, c(m1,m2, m3, m4),  as.character(c(m1,m2,m3,m4)), col.axis = 4)
legend("topright", pch=c(1,2,3,4), c("AIC", "BIC", "QHC", "LR"), bty="n", col=c(1,2,3,4))

## -----------------------------------------------------------------------------
summary(vaal)

## -----------------------------------------------------------------------------
p<-vaal$p
p

## -----------------------------------------------------------------------------
data(chicken.embryo)
head(chicken.embryo, 2)

## ----results = "hide", warning = FALSE----------------------------------------
a <- 0
b <- 21
day <- chicken.embryo$day
nT <- chicken.embryo$nT
embryo<-mable.group(x = nT, breaks = a:b, M=c(2,100), interval = c(a, b), IC = "aic",
    controls = mable.ctrl(sig.level = 1e-6,  maxit = 2000, eps = 1.0e-7))

## ----fig.align='center', fig.cap="Chicken Embryo Data\\label{fig:chicken-embryo-data-plot}", fig.width=6, fig.height=7----
Day <- rep(day,nT)
op <- par(mfrow = c(1,2), lwd = 2)
layout(rbind(c(1, 2), c(3, 3)))
plot(embryo, which = "likelihood")
plot(embryo, which = "change-point")
fk <- density(x = rep((0:20)+.5, nT), bw = "sj", n = 101, from = a, to = b)
hist(Day, breaks = seq(a,b,  length = 12), freq = FALSE, col = "grey", border = "white", 
     main = "Histogram and Density Estimates")
plot(embryo, which = "density", cex = 0.7, add = TRUE)
lines(fk, lty = 2, col = 2)
legend("top", lty = c(1:2), c("MABLE", "Kernel"), bty = "n", col = c(1:2))
par(op)

## ----fig.align='center', fig.cap="AIC and BIC based on Chicken Embryo Data\\label{fig:chicken-embryo-data-AIC-BIC-plot}", fig.width=6, fig.height=4----
M <- embryo$M[1]:embryo$M[2]
aic <- embryo$ic$AIC
bic <- embryo$ic$BIC
res.gcp <- optim.gcp(embryo) # choose m by gamma change-point model
lr <- res.gcp$lr
plot(M, aic, cex = 0.7, col = 1, xlab = "m", ylab = "", ylim = c(ymin<-min(aic,bic,lr),
    ymax<-max(aic,bic,lr)), main = "AIC, BIC, and LR")
points(M, bic, pch = 2, cex = 0.7, col = 2)
points(M[-1], lr, pch = 3, cex = 0.7, col = 4)
segments(which.max(aic)+M[1]-1->m1, ymin, m1, max(aic), lty = 2)
segments(which.max(bic)+M[1]-1->m2, ymin, m2, max(bic), lty = 2, col = 2)
segments(which.max(lr)+M[1]->m3, ymin, m3, max(lr), lty = 2, col = 4)
axis(1, c(m1,m2, m3),  as.character(c(m1,m2,m3)), col.axis = 4)
legend("right", pch = 1:3, c("AIC", "BIC", "LR"), bty = "n", col = c(1,2,4))

## -----------------------------------------------------------------------------
summary(embryo)

## ----results = "hide", warning = FALSE----------------------------------------
set.seed(123)
mu <- 1; sig <- 2; a <- mu - sig*5; b <- mu + sig*5;
gn <- function(x) dnorm(x, 0, 1)
n <- 50;
x <- rnorm(n, mu, sig); e <- rnorm(n); y <- x + e;
decn <- mable.decon(y, gn, interval = c(a,b), M = c(5, 50))

## ----fig.align='center', fig.cap="Simulated Normal Data\\label{fig:simulated-normal-data-plot}", fig.width=7, fig.height=7----
op <- par(mfrow = c(2,2), lwd = 2)
plot(decn, which = "likelihood")
plot(decn, which = "change-point", lgd.x = "right")
plot(xx<-seq(a, b, length = 100), yy<-dnorm(xx, mu, sig), type = "l", xlab = "x",
     ylab = "Density", ylim = c(0, max(yy)*1.1))
plot(decn, which = "density", add = TRUE, lty = 2, col = 2)
# kernel density based on pure data
lines(density(x), lty = 5, col = 4)
legend("topright", bty = "n", lty = c(1,2,5), col = c(1,2,4), c(expression(f),     
    expression(hat(f)),  expression(tilde(f)[K])))
plot(xx, yy<-pnorm(xx, mu, sig), type = "l", xlab = "x", ylab = "Distribution Function")
plot(decn, which = "cumulative", add = TRUE, lty = 2, col = 2)
legend("bottomright", bty = "n", lty = c(1,2), col = c(1,2), c(expression(F), 
    expression(hat(F))))
par(op)

## ----results = "hide", warning = FALSE----------------------------------------
library(coxinterval) 

## -----------------------------------------------------------------------------
bcos=cosmesis
head(bcos, 3)

## ----results = "hide"---------------------------------------------------------
bc.res0 <- mable.ic(bcos[bcos$treat == "RT",1:2], M = c(1,50), IC = "none")
bc.res1 <- mable.ic(bcos[bcos$treat == "RCT",1:2], M = c(1,50), IC = "none")

## ----message=FALSE, warning=FALSE, fig.align='center', fig.cap="Breast Cosmesis Data\\label{fig:Breast-Cosmesis-Data-plot}", fig.width=7, fig.height=7, warning = FALSE----
op <- par(mfrow = c(2,2),lwd = 2)
plot(bc.res0, which = "change-point", lgd.x = "right")
plot(bc.res1, which = "change-point", lgd.x = "right")
plot(bc.res0, which = "survival", xlab = "Months", ylim = c(0,1), main = "Radiation Only")
plot(bc.res1, which = "survival", xlab = "Months", main = "Radiation and Chemotherapy")
par(op)

## -----------------------------------------------------------------------------
data(faithful)
head(faithful, 3)

## ----message=FALSE, warning=FALSE, results = "hide"---------------------------
a <- c(0, 40); b <- c(7, 110)
#faith2 <- mable.mvar(faithful, M = c(60,30), interval = cbind(a,b))
faith2 <- mable.mvar(faithful, M = c(46,19), search =FALSE, interval = cbind(a,b))

## ----fig.align='center', fig.cap="Density Estimate based on the Old Faithful Data\\label{fig:Old-Faithful-Data-plot}", fig.width=6, fig.height=6, warning = FALSE----
plot(faith2, which = "density")

## -----------------------------------------------------------------------------
summary(faith2)

## -----------------------------------------------------------------------------
futime2 <- ovarian$futime
futime2[ovarian$fustat==0] <- Inf
ovarian2 <- data.frame(age = ovarian$age, futime1 = ovarian$futime, futime2 = futime2)
head(ovarian2, 3)

## ----results = "hide", warning = FALSE, message=FALSE-------------------------
ova<-mable.ph(cbind(futime1, futime2) ~ age, data = ovarian2, M = c(2,35), g = .16, x0=35)

## -----------------------------------------------------------------------------
summary(ova)

## ----fig.align='center', fig.cap="Ovarian Cancer Data\\label{fig:ovarian-Data-plot}", fig.width=7, fig.height=7, warning = FALSE, message=FALSE----
op <- par(mfrow = c(2,2))
plot(ova, which = "likelihood")
plot(ova, which = "change-point")
plot(ova, y=data.frame(age=60), which="survival", type="l", xlab="Days", main="Age = 60")
plot(ova, y=data.frame(age=65), which="survival", type="l", xlab="Days", main="Age = 65")
par(op)

## ----results = "hide", warning = FALSE----------------------------------------
ova1 <- mable.reg(cbind(futime1, futime2) ~ age, data = ovarian2, M = c(2,35))

## -----------------------------------------------------------------------------
bcos2 <- data.frame(bcos[,1:2], x = 1*(bcos$treat == "RCT"))

## ----results = "hide", warning=FALSE, message=FALSE---------------------------
g <- 0.41 # Hanson and Johnson 2004, JCGS
aft.res<-mable.aft(cbind(left, right) ~ x, data = bcos2, M =c(1, 30), g=.41,  tau =100, x0=1)

## -----------------------------------------------------------------------------
summary(aft.res)

## ----fig.align='center', fig.cap="AFT Model Fit for Breast Cosmesis Data\\label{fig:Breast-Cosmesis-Data-aft-plot}", fig.width=7, fig.height=4, warning = FALSE----
op <- par(mfrow = c(1,2), lwd = 1.5)
plot(x = aft.res, which = "likelihood")
plot(x = aft.res, y = data.frame(x = 0), which = "survival", type = "l", col = 1, 
    main = "Survival Function")
plot(x = aft.res, y = data.frame(x = 1), which = "survival", lty = 2, col = 1, add = TRUE)
legend("bottomleft", bty = "n", lty = 1:2, col = 1, c("Radiation Only", 
    "Radiation and Chemotherapy"), cex = .7)
par(op)

## ----results = "hide"---------------------------------------------------------
aft.res1 <- mable.reg(cbind(left, right) ~ x, data = bcos2, 'aft', M = c(1, 30), 
      g=.41, tau=100, x0=1)

## ----message=FALSE, warning = FALSE, results = "hide"-------------------------
# Hosmer and Lemeshow (1989):                                           
# ages and the status of coronary disease (CHD) of 100 subjects         
x<-c(20, 23, 24, 25, 26, 26, 28, 28, 29, 30, 30, 30, 30, 30, 32,        
32, 33, 33, 34, 34, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39,         
40, 41, 41, 42, 42, 42, 43, 43, 44, 44, 45, 46, 47, 47, 48, 49,         
49, 50, 51, 52, 55, 57, 57, 58, 60, 64)                                 
y<-c(25, 30, 34, 36, 37, 39, 40, 42, 43, 44, 44, 45, 46, 47, 48,        
48, 49, 50, 52, 53, 53, 54, 55, 55, 56, 56, 56, 57, 57, 57, 57,         
58, 58, 59, 59, 60, 61, 62, 62, 63, 64, 65, 69) 
a<-20; b<-70
regr<-function(x) cbind(1,x)                                            
chd.mable<-mable.dr(x, y, M=c(1, 15), regr, interval = c(a,b))   

## ----fig.align='center', fig.cap="DR Model Fit for Coronary Heart Disease  Data\\label{fig:CHD-data-plot}", fig.width=7, fig.height=4, warning = FALSE----
z<-seq(a,b,length=512)
f0hat<-dmixbeta(z, p=chd.mable$p, interval=c(a, b))
rf<-function(x) chd.mable$regr((x-a)/(b-a))
f1hat<-dtmixbeta(z, p=chd.mable$p, alpha=chd.mable$alpha, 
                 interval=c(a, b), regr=rf)
op<-par(mfrow=c(1,2),lwd=1.2, cex=.7, mar=c(5,4,1,1))
hist(x, freq=F, col = "light grey", border = "white", xlab="Age", 
  ylab="Density", xlim=c(a,b), ylim=c(0,.055), main="Control")
lines(z, f0hat, lty=1, col=1)
hist(y, freq=F, col = "light grey", border = "white", xlab="Age", 
  ylab="Density", xlim=c(a,b), ylim=c(0,.055), main="Case")
lines(z, f1hat, lty=1, col=1)
par(op)

## ----message=FALSE, warning = FALSE, results = "hide"-------------------------
data(pancreas)
head(pancreas,3)
x<-log(pancreas$ca199[pancreas$status==0])
y<-log(pancreas$ca199[pancreas$status==1])
a<-min(x,y); b<-max(x,y)
M<-c(1,29)
regr<-function(x) cbind(1,x,x^2)                                            
m=maple.dr(x, y, M, regr=regr, interval=c(a,b), controls=mable.ctrl(sig.level=.001))$m
pc.mable<-mable.dr(x, y, M=m, regr=regr, interval=c(a,b),
                  controls=mable.ctrl(sig.level=1/length(c(x,y))))
#pc.mable   

## ----fig.align='center', fig.cap="DR Model Fit for Pancreatic Cancer Biomarker Data\\label{fig:pcb-data-plot}", fig.width=7, fig.height=3, warning = FALSE----
z<-seq(a,b,length=512)
# baseline is "case"
f1hat<-dmixbeta(z, p=pc.mable$p, interval=c(a, b))
rf<-function(x) pc.mable$regr((x-a)/(b-a))
f0hat<-dtmixbeta(z, p=pc.mable$p, alpha=pc.mable$alpha, 
                 interval=c(a, b), regr=rf)
op<-par(mfrow=c(1,2),lwd=1.2, cex=.7, mar=c(5,4,1,1))
hist(x, freq=F, col = "light grey", border = "white", xlab="log(CA19-9)", 
  ylab="Density", xlim=c(a,b),  main="Control")
lines(z, f0hat, lty=1, col=1)
hist(y, freq=F, col = "light grey", border = "white", xlab="log(CA19-9)", 
  ylab="Density", xlim=c(a,b), ylim=c(0,.5), main="Case")
lines(z, f1hat, lty=1, col=1)
par(op)

