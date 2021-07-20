rm(list=ls())

library(Rfast)
library(emplik)
library(limSolve)
library(cvTools)

#setwd("")
source("example1codeCI.R")

###########################example1
dat1 <- function(n=1e4, p=7, case=1){
  beta0  <- rep(0.2, p)
  ind <- 1:p
  
  sigmax  <- 0.2 ^abs(outer(ind, ind, "-"))
  X  <-Rfast::rmvnorm(n, rep(0, p), sigmax)
  switch(case,
         {    ## 1:  normal mean 0
           Y<- X %*% beta0 + rnorm(n)
         }, { ## 2: t10
           Y<- X %*% beta0 + rt(n, df=10)
         }, { ## 3: mix nnorm
           components <- sample(1:2,prob=c(0.5, 0.5),size=n,replace=TRUE)
           mus <- c(1, -1)
           sds <- sqrt(c(1, 1))
           e <- rnorm(n=n,mean=mus[components],sd=sds[components])
           Y<- X %*% beta0 + e
         })

  return(list(y=Y, x=X, beta=beta0))
}

#######################
case <- 1 # 1,2,3
######################
n <- 1e5
p <- 7
m <- 500 # repeated time


#######BLB and SDB setting begin
S.BLB <- 20
S.SDB <- 500
r <- 100
b06 <- ceiling(n^0.6)
b08 <- ceiling(n^0.8)
#######BLB and SDB setting end

aa <- 1000

p.DACEL50 <- matrix(aa, m, p)
p.DACEL100 <- matrix(aa, m, p)
p.DACEL150 <- matrix(aa, m, p)

p.BLB06 <- matrix(aa, m, p)
p.BLB08 <- matrix(aa, m, p)

p.SDB06 <- matrix(aa, m, p)
p.SDB08 <- matrix(aa, m, p)

p.boot <- matrix(aa, m, p)

CI.BLB06 <- matrix(aa, m, p)
CI.BLB08 <- matrix(aa, m, p)

CI.SDB06 <- matrix(aa, m, p)
CI.SDB08 <- matrix(aa, m, p)

CI.boot <- matrix(aa, m, p)

CI.DACEL50  <- matrix(aa, m, p)
CI.DACEL100 <- matrix(aa, m, p)
CI.DACEL150 <- matrix(aa, m, p)


for(i in 1:m){
  dat <- dat1(n=n, p=p, case=case)
  x <- dat$x
  y <- dat$y
  beta0 <- dat$beta
  beta  <- beta0 ###### Empirical size
  #beta <- NULL  ###### Empirical power
  
  
  ##################TB begin
  fit.lm <- lm(y~x-1)
  
  fit.boot <- lm.boot(x=x, y=y, B=100, beta0=beta, alpha=0.05)
  CI.boot[i, ] <- fit.boot$coef.UL[2, ] - fit.boot$coef.UL[1,]
  p.boot[i, ]  <- fit.boot$p.coef
  ##################TB end
  
  #################DACEL begin K=50 
  fit.mog50 <- lm.moel(x=x, y=y,  beta0=beta, K=50, type="random", Intercept=FALSE)
  p.DACEL50[i, ] <- fit.mog50$p.coefel
  CI.DACEL50[i,] <- fit.mog50$lengthel.UL
  #################DACEL end K=50
  
  #################DACEL begin K=100
  fit.mog100 <- lm.moel(x=x, y=y,  beta0=beta, K=100, type="random",  Intercept=FALSE)
  p.DACEL100[i, ] <- fit.mog100$p.coefel
  CI.DACEL100[i,] <- fit.mog100$lengthel.UL
  #################DACEL end K=100
  
  #################DACEL begin K=150
  fit.mog150 <- lm.moel(x=x, y=y,  beta0=beta, K=150, type="random",  Intercept=FALSE)
  p.DACEL150[i, ] <- fit.mog150$p.coefel
  CI.DACEL150[i,] <- fit.mog150$lengthel.UL
  #################DACEL end K=150
  
  
  
  ####################################################BLB and SDB begin
  fit.BLB06 <- BLB.lm(x=x, y=y, S=S.BLB, r=r, b=b06, beta0=beta, alpha=0.05)
  CI.BLB06[i, ] <- fit.BLB06$coef.ULBLB[2, ] - fit.BLB06$coef.ULBLB[1, ]
  p.BLB06[i, ] <- fit.BLB06$p.coefBLB
  
  fit.SDB06 <- SDB.lm(x=x, y=y, S=S.SDB, b=b06, beta0=beta, alpha=0.05)
  CI.SDB06[i, ] <- fit.SDB06$coef.ULSDB[2, ] - fit.SDB06$coef.ULSDB[1, ]
  p.SDB06[i, ] <- fit.SDB06$p.coefSDB
  ####################################################BLB and SDB end
  
  ####################################################BLB and SDB begin
  fit.BLB08 <- BLB.lm(x=x, y=y, S=S.BLB, r=r, b=b08, beta0=beta, alpha=0.05)
  CI.BLB08[i, ] <- fit.BLB08$coef.ULBLB[2, ] - fit.BLB08$coef.ULBLB[1, ]
  p.BLB08[i, ] <- fit.BLB08$p.coefBLB
  
  fit.SDB08 <- SDB.lm(x=x, y=y, S=S.SDB, b=b08, beta0=beta, alpha=0.05)
  CI.SDB08[i, ] <- fit.SDB08$coef.ULSDB[2, ] - fit.SDB08$coef.ULSDB[1, ]
  p.SDB08[i, ] <- fit.SDB08$p.coefSDB
  ####################################################BLB and SDB end

  print(i)
}

name.model <- c( "DACEL50", "DACEL100", "DACEL150",
                 "BLB06", "BLB08", "SDB06","SDB08", "boot", 
                 "CIDACEL50", "CIDACEL100", "CIDACEL150",
                 "CIBLB06", "CIBLB08", "CISDB06", "CISDB08", "CIboot")

result.p <- matrix(aa, length(name.model), p)

rownames(result.p) <- name.model 
colnames(result.p) <- c(paste0('x', 1:p))


##################  empirical size
result.p["DACEL50",1:p] <- apply(p.DACEL50, 2, function(x)(mean(I(x <= 0.05))))
result.p["DACEL100",1:p] <- apply(p.DACEL100, 2, function(x)(mean(I(x <= 0.05))))
result.p["DACEL150",1:p] <- apply(p.DACEL150, 2, function(x)(mean(I(x <= 0.05))))

result.p["BLB06",1:p] <- apply(p.BLB06, 2, function(x)(mean(I(x <= 0.05))))
result.p["BLB08",1:p] <- apply(p.BLB08, 2, function(x)(mean(I(x <= 0.05))))

result.p["SDB06",1:p] <- apply(p.SDB06, 2, function(x)(mean(I(x <= 0.05))))
result.p["SDB08",1:p] <- apply(p.SDB08, 2, function(x)(mean(I(x <= 0.05))))

result.p["boot",1:p] <- apply(p.boot, 2, function(x)(mean(I(x <= 0.05))))

###################CI length
result.p["CIDACEL50", 1:p] <- apply(CI.DACEL50, 2, mean) 
result.p["CIDACEL100", 1:p] <- apply(CI.DACEL100, 2, mean) 
result.p["CIDACEL150", 1:p] <- apply(CI.DACEL150, 2, mean) 

result.p["CIBLB06", 1:p] <- apply(CI.BLB06, 2, mean)
result.p["CIBLB08", 1:p] <- apply(CI.BLB08, 2, mean)

result.p["CISDB06", 1:p] <- apply(CI.SDB06, 2, mean)
result.p["CISDB08", 1:p] <- apply(CI.SDB08, 2, mean)

result.p["CIboot", 1:p] <- apply(CI.boot, 2, mean)

print(case)
result.p






