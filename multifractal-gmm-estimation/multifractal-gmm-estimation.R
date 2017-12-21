library(Matrix)
library(compiler)
library(ACDm)
library(numDeriv)
library(moments)

############################# Functions ##############################
## Some global parameters
init <- 1
j <- 1
k=7
invReweighted <- 1
S <- 1


## Parameter vector theta = (b,g,lam,mu)
TransP <- function(b,g,k,j){1-(1-g)^(b^(j-k))}


Lambda <- function(n,b,g,lam, mu,k){
  ts <- c(1:n)
  init <<- 1 + numeric(k)
  first <- lam*prod(init)
  ts2 <- sapply(ts[2:n],function(x){
    Ms <- c(1:k)
    Ms <- sapply(Ms,function(x){
      P <- TransP(b,g,k,x)
      Rand <- runif(1,0,1)
      if(P > Rand){
        rlnorm(1,-mu,sqrt(2*mu))
      } else {
        init[x]
      }
    })
    init <<- Ms
    lam*prod(Ms)
  })
  return(c(first,ts2))
}

cLambda <- cmpfun(Lambda)

Dur <- function(n,theta,k){
  b <- theta[1]
  g <- theta[2]
  lam <- theta[3]
  mu <- theta[4]
  exps <- rexp(n,1)
  lambda <- cLambda(n,b,g,lam,mu,k)
  return(exps/lambda)
}

NewAuto <- function(t,b,g,lam,mu,k){
  Vec <- numeric(k)
  Ed <- exp(2*mu)
  Ed2 <- exp(6*mu)
  for(j in 1:k){
    gam <- TransP(b,g,k,j)
    Vec[j] <- (1-(1-gam)^t) * Ed^2 + (1-gam)^t * Ed2}
  Vec <- prod(Vec)
  if(t == 0){e = 2}
  else{e = 1}
  return((e*Vec)/lam^2)
}

NewAuto22 <- function(t,b,g,lam,mu,k){
  Vec <- numeric(k)
  Ed <- exp(2*mu)
  Ed2 <- exp(6*mu)
  Ed4 <- exp(20*mu)
  for(j in 1:k){
    gam <- TransP(b,g,k,j)
    Vec[j] <- (1-(1-gam)^t) * Ed2^2 + (1-gam)^t * Ed4}
  Vec <- prod(Vec)
  if(t == 0){e = 2}
  else{e = 1}
  return(4*(e*Vec)/lam^4)
}

SampleAuto <- function(lag,data,n=1){
  len <- length(data)
  r <- sum(data[(1+lag):len]^n*data[1:(len-lag)]^n)/(len-lag)
  return(r)
}

DurMom <- function(n,lam,mu,k){
  factorial(n)*exp(n*(n+1)*mu*k)/lam^n
}

## One half of the g-function in GMM, giving the analytic moments
gfunc <- function(theta){
  b <- theta[1]
  g <- theta[2]
  lam <- theta[3]
  mu <- theta[4]
  Moments <- c(
    DurMom(1,lam,mu,k),
    DurMom(2,lam,mu,k),
    NewAuto(1,b,g,lam,mu,k),
    NewAuto(10,b,g,lam,mu,k),
    NewAuto22(5,b,g,lam,mu,k))
  return(Moments)
}

## Other half of the g-function in GMM, giving sample moments.
Sample <- function(data){
  len <- length(data)
  Moments <- c(
    (sum(data)/len),
    (sum(data^2)/len),
    SampleAuto(1,data),
    SampleAuto(10,data),
    SampleAuto(5,data,n=2))
  return(Moments)
}


## Optimized function to calculate W in 2-step GMM
hfunc <- function(theta,data){
  len <- length(data)
  maxlag <- 200
  Moments <- gfunc(theta)
  j <<- 1
  mat <- lapply(data[1:(len-200)], function (x) {
    j <<- j + 1
    SampleMom <- c(
      x,
      x^2,
      x*data[j],
      x*data[j+9],
      x^2*data[j+4]^2)
    Mom <- Moments - SampleMom
    return(tcrossprod(Mom))
  }
  )
  return(mat)
}

cfunc <- cmpfun(hfunc)

## Objective function for first step in 2-step GMM
f1 <- function(theta){
  M <- gfunc(theta)
  Q <- crossprod(M-S)
  return(Q)
}

## Objective function for second step in 2-step GMM
f2 <- function(theta){
  
  M <- gfunc(theta)
  Q <- t(M-S) %*% (invReweighted %*% (M-S))
  return(Q)
}

ReWeight <- function(Data,est){
  interm <- cfunc(est,Data)
  sumInterm <- Reduce('+',interm)
  Reweighted <- sumInterm/(length(Data)-20)
  invReweighted <<- solve(Reweighted)
  invReweighted
}

## Calculates 2-step GMM for Data, requires bounds on b.
GMMest <- function(Data){
  j <<- 1
  S <<- Sample(Data)
  firstest <- nlminb(c(3,0.8,2,0.01),f1,lower=c(1,0,0,0),upper=c(100,0.95,20,2),control=list(eval.max=100000,iter.max=100000))
  invReweighted <- ReWeight(Data,firstest$par)
  secondest <- nlminb(firstest$par,f2,lower=c(1,0,0,0),upper=c(100,0.95,20,2),control=list(eval.max=100000,iter.max=100000))
  return(secondest$par)
}
library(numDeriv)
VarCoVar <- function(par,mat){
  solve(t(jacobian(gfunc,par)) %*% mat %*% jacobian(gfunc,par))
}
#################### Monte Carlo Simulations #########################

set.seed(20458543)
bs <- numeric(1000)
gs <- numeric(1000)
lams <- numeric(1000)
mus <- numeric(1000)
for(i in 1:1000){
  data <- Dur(20000,c(2,0.5,3,0.05),7)
  est <- GMMest(data)
  bs[i] <- est[1]
  gs[i] <- est[2]
  lams[i] <- est[3]
  mus[i] <- est[4]
  print(i)
}

########## Importing and working with data #################

data <- read.csv("HFdurs.csv")
len <- length(data$times2)
bn <- ceiling(max(diff(data$durations)))

## Compute Nadarya-Watson estimates for relevant points
ndw <- ksmooth(data$times2[1:len],
               data$durations[1:len],
               kernel="normal",
               bandwidth = bn,
               n.points  = max(data$times2)+1)

## Relace 0 with suitably small value
durs <- sapply(data$durations,function(x){
  if(x==0){
    x = 0.0009
  } else {x = x}
})


## Compute weights to remove seasonality
weightz <- numeric(len)
for(i in 1:(length(data$times2[1:len]))){
  j = data$times2[i]
  weightz[i] <- ndw$y[j + 1]
}

# d*, the transformed data
newDur <- durs/weightz

# Estimates for duration
parest <- GMMest(newDur)

# Variance-covaraince matrix for estimates
VCV <- VarCoVar(parest,invReweighted)/length(newDur)

# ACD Fit

Fit <- acdFit(newDur)

#################### Forecasting ###########################

n=10000
a <- numeric(100000)
for(i in 1:100000){
  a[i] = NewAuto(i,parest[1],parest[2],parest[3],parest[4],7) - DurMom(1,parest[3],parest[4],7)
}
a <- c(DurMom(2,parest[3],parest[4],7),a)

newDur0 <- newDur - mean(newDur)

n=10000
vs <- numeric(n)
tempphis <- numeric(n)
phis <- numeric(n)
phisplus <- 0

vs[1] <- a[1]
tempphis[1] <- a[2]/a[1]

for(m in 1:n){
  if(vs[m] == 0){phisplus <- 0}
  else{phisplus <- (a[m+1] - (tempphis[1:m] %*% rev(a[2:(m+1)]))[1,1])/vs[m]}
  phis <- tempphis[1:m] - phisplus*rev(tempphis[1:m])
  vs[m+1] <- (1-phisplus^2)*vs[m]
  tempphis <- c(phis,phisplus)
}
n = n+1
vsh <- vs[n]
phish <- matrix(nrow = n,ncol=30)
phish[,1] <- tempphis

for(h in 2:30){
  phish[,h][n] <- (a[n+h] - phis %*% a[(h+1):(h+n-1)])/vs[n-1]
  phish[,h][1:(n-1)] <- phish[,(h-1)][2:(n+1)] + phish[,1][1]*phis - phish[,h][n]*rev(phis)
  vsh <- vsh + (phish[,(h-1)][1]^2 + phish[,h][n]^2)*vs[n-1]
}

# Partitions of duration, usually slices list [1:i], for i in training2
training <- seq(1,2147600,1)

# Reverse the zero mean process for easier list slicing
Actual <- rev(newDur0)

# MSE for Naive forecast for lag h. Use mean or desired statistic on output
Naive <- function(h){
  trainingh <- seq(1+h,2147600+h,1)
  return((Actual[trainingh]-Actual[training])^2)
}

# Calculate MSE for MSMD, returns a list. Use mean or desired statistic on output
forecastMSMD <- function(h){
  coef <- phish[,h]
  data <- Actual[(1+h):(2147600+h)]
  Forecasts <- numeric(213760)
  j = 0
  for(i in c(1:213760)){
    j = j + 1
    temp <- data[(1+j):(n+j)]
    Forecasts[j] <- temp %*% coef
  }
  return((Forecasts - Actual[1:213760])^2)
}

# Helper to recursively calculate forecasts
ForecastHelper <- function(h,w,a,b,init){
  if(h==1){w+(a+b)*init}
  else{w+(a+b)*ForecastHelper(h-1,w,a,b,init)}
}


# Static Window Forecasting for ACD

fits <- acdFit(newDur[1:1000000])$par
w = fits$Coef[1]
a = fits$Coef[2]
b = fits$Coef[3]

# Calculate MSE for ACD
ForecastACD <- function(h){
  fits <- acdFit(newDur[1:100000],output=FALSE)$par
  w = fits$Coef[1]
  a = fits$Coef[2]
  b = fits$Coef[3]
  act <- rev(newDur)
  data <- act[(1+h):(2147600+h)]
  Forecasts <- numeric(length(data))
  j = 0
  for(i in data){
    j = j + 1
    Forecasts[j] <- ForecastHelper(h,w,a,b,i)
  }
  return((act[1:2147600]-Forecasts)^2)
}
mean(ForecastACD(1))

mean(ForecastACD(5))

mean(ForecastACD(10))

## Results for Forecasting

## MSE Ratios

for(i in c(1,5,10,15,20)){
  print(mean(forecastMSMD(i))/mean(Naive(i)))
}

for(i in c(1,5,10,15,20)){
  print(mean(ForecastACD(i))/mean(Naive(i)))
}

## MAE Ratios

for(i in c(1,5,10,15,20)){
  print(mean(sqrt(forecastMSMD(i)))/mean(sqrt(Naive(i))))
}

for(i in c(1,5,10,15,20)){
  print(mean(sqrt(ForecastACD(i)))/mean(sqrt(Naive(i))))
}

## Sample Forecasts

fits <- acdFit(newDur[1:1000000])$par
w = fits$Coef[1]
a = fits$Coef[2]
b = fits$Coef[3]

Dat <- newDur[1:1230030]
Dat0 <- newDur0[1:1230000]
ACDSampleF <- c(1:20)
for(i in 1:20){
  ACDSampleF[i] <- ForecastHelper(i,w,a,b,newDur[1230000])
}

MSMDSampleF <- c(1:20)
for(i in 1:20){
  MSMDSampleF[i] <- phish[,i] %*% Dat0[1220000:1230000]
}
MSMDSampleF <- MSMDSampleF + mean(newDur)

NaiveSampleF <- rep(newDur[1230000],20)

## Summary Statistics
mean(newDur)
var(newDur)
skewness(newDur)
kurtosis(newDur)

## Plots
# All data
ts.plot(newDur)
# Subset of data
ts.plot(newDur[19900:20000],ylab="Duration")

# Sample Forecasts
ts.plot(Dat[1229950:1230020],ylab="Duration")
lines(MSMDSampleF,x=seq(51,70,1),col="blue")

lines(ACDSampleF,x=seq(51,70,1),col="red")

lines(NaiveSampleF,x=seq(51,70,1),col="green")

# Density plots
musd <- density(mus)
lamsd <- density(lams)
gsd <- density(gs)
bsd <- density(bs)
plot(musd,ylim=c(0,max(musd$y)),xlim=c(min(musd$x),max(musd$x)),xlab="",ylab="",main="")
par(new=TRUE)
hist(mus,freq=FALSE,ylim=c(0,max(musd$y)),xlim=c(min(musd$x),max(musd$x)))
plot(lamsd,ylim=c(0,max(lamsd$y)),xlim=c(min(lamsd$x),max(lamsd$x)),xlab="",ylab="",main="")
par(new=TRUE)
hist(lams,freq=FALSE,ylim=c(0,max(lamsd$y)),xlim=c(min(lamsd$x),max(lamsd$x)))
plot(gsd,ylim=c(0,max(gsd$y)),xlim=c(min(gsd$x),max(gsd$x)),xlab="",ylab="",main="")
par(new=TRUE)
hist(gs,freq=FALSE,ylim=c(0,max(gsd$y)),xlim=c(min(gsd$x),max(gsd$x)))
plot(bsd,ylim=c(0,max(bsd$y)),xlim=c(min(bsd$x),max(bsd$x)),xlab="",ylab="",main="")
par(new=TRUE)
hist(bs,freq=FALSE,ylim=c(0,max(bsd$y)),xlim=c(min(bsd$x),max(bsd$x)))

# Sample MSMD plots
SMSMD <- Dur(20000,c(2,0.5,3,0.05),7)
ts.plot(SMSMD,ylab="Duration")
ts.plot(SMSMD[10000:10100],ylab="Duration")

# Transition Probability plots
for(i in 3:8){
  temp <- numeric(7)
  for(j in 1:7){
    temp[j] <- TransP(i/2,5/(2*i),7,j)
  }
  if(i == 3){ts.plot(temp,ylim = c(0,1),xlab="Component",ylab="Transition Probability")}
  else{lines(temp)}
}