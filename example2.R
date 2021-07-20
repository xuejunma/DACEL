rm(list=ls())

library(Rfast)
library(emplik)
library(limSolve)
library(cvTools)

#setwd("")
source("example2codeCI.R")

###################### example 2
dat1 <- function(n=1e4, p=7, case=1){
  beta0  <- rep(0.2, p)
  corr  <- 0.5
  sigmax  <- matrix(corr, p, p) + diag(1-corr, p)
  switch(case,
         {    ## case 1:  normal mean 0
           X  <- rmvnorm(n, rep(0, p), sigmax)
         }, { ## case 2: normal mean 1 variance 1
           X  <- rmvnorm(n, rep(1.5, p), sigmax)
         },  { ## case 3: mixnormal means 1 and -1
           X  <- rmvnorm(n/2, rep(1, p), sigmax/4)
           X  <- rbind(rmvnorm(n/2, rep(-1, p), sigmax/4), X)
         }, { ## case 4: t3
           X  <- Rfast::rmvt(n, rep(0, p), sigmax, 3) / 10
         }, { ## case 5: exponential
           X <- matrix(rexp(n*p, 2), n, p)
        }, { ## case 6: mixnormal means -2.14 and -2.9
             X  <- rmvnorm(n/2, rep(-2.14, p), sigmax/4)
             X  <- rbind(rmvnorm(n/2, rep(-2.9, p), sigmax/4), X)
         })
  P  <- 1 - 1 / (1 + exp(X %*% beta0))
  Y  <- rbinom(n, 1, P)
  return(list(y=Y, x=X, beta=beta0))
}

#######################
case <- 1
######################
n <- 1e5
p <- 7
m <- 500 # repeated time

#######subsample setting begin
r1 <- 2000
r2 <- 0.1*n
#######subsample setting end


#######BLB and SDB setting begin
S.BLB <- 20
S.SDB <- 500
r <- 100
b06 <- ceiling(n^0.6)
b08 <- ceiling(n^0.8)
#######BLB and SDB setting end


aa <- 1000

p.mog50 <- matrix(aa, m, p)
p.mog100 <- matrix(aa, m, p)
p.mog150 <- matrix(aa, m, p)

p.mvc <- matrix(aa, m, p)
p.mmse <- matrix(aa, m, p)

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

CI.mog50 <- matrix(aa, m, p)
CI.mog80 <- matrix(aa, m, p)
CI.mog100 <- matrix(aa, m, p)
CI.mog150 <- matrix(aa, m, p)

CI.mvc <- matrix(aa, m, p)
CI.mmse <- matrix(aa, m, p)

###################
for(i in 1:m){
  dat <- dat1(n=n, p=p, case=case)
  x <- dat$x
  y <- dat$y
  beta0 <- dat$beta
  beta  <- beta0 ###### Empirical size
  #beta <- NULL  ###### Empirical power
  
  ##################TB begin
  fit.glm <- glm(y~x-1, family=binomial(link = "logit"))

  fit.boot <- glm.boot(x=x, y=y, B=100, beta0=beta, alpha=0.05)
  CI.boot[i, ] <- fit.boot$coef.UL[2, ] - fit.boot$coef.UL[1,]
  p.boot[i, ]  <- fit.boot$p.coef
  ##################TB  end
  
  
  #################MOG BEGIN K=50
  fit.mog50 <- glm.moel(x=x, y=y,  beta0=beta, K=50, type="random", family=binomial(link = "logit"), Intercept=FALSE)
  p.mog50[i, ] <- fit.mog50$p.coefel
  CI.mog50[i,] <- fit.mog50$lengthel.UL
  #################MOG end K=50
  
  #################MOG BEGIN K=100
  fit.mog100 <- glm.moel(x=x, y=y,  beta0=beta, K=100, type="random", family=binomial(link = "logit"), Intercept=FALSE)
  p.mog100[i, ] <- fit.mog100$p.coefel
  CI.mog100[i,] <- fit.mog100$lengthel.UL
  
  #################MOG BEGIN K=150
  fit.mog150 <- glm.moel(x=x, y=y,  beta0=beta, K=150, type="random", family=binomial(link = "logit"), Intercept=FALSE)
  p.mog150[i, ] <- fit.mog150$p.coefel
  CI.mog150[i,] <- fit.mog150$lengthel.UL
  #################MOG k=150 end
  
  #####################mmse begin
  fit.mmse <- twostep(X=x, Y=y,  beta0=beta, r1=r1, r2=r2, method="mmse")
  
  p.mmse[i, ] <- fit.mmse$p.each
  CI.mmse[i, ] <-  fit.mmse$length.UL
  #####################mmse end
  
  #####################mvc begin
  fit.mvc <- twostep(X=x, Y=y,  beta0=beta, r1=r1, r2=r2, method="mvc")
  
  p.mvc[i, ] <- fit.mvc$p.each
  CI.mvc[i, ] <-  fit.mvc$length.UL
  #####################mvc end
  
  ####################################################BLB and SDB begin
  fit.BLB06 <- BLB.glm(x=x, y=y, S=S.BLB, r=r, b=b06, beta0=beta, alpha=0.05)
  CI.BLB06[i, ] <- fit.BLB06$coef.ULBLB[2, ] - fit.BLB06$coef.ULBLB[1, ]
  p.BLB06[i, ]  <- fit.BLB06$p.coefBLB
  
  fit.SDB06 <- SDB.glm(x=x, y=y, S=S.SDB, b=b06, beta0=beta, alpha=0.05)
  CI.SDB06[i, ] <- fit.SDB06$coef.ULSDB[2, ] - fit.SDB06$coef.ULSDB[1, ]
  p.SDB06[i, ]  <- fit.SDB06$p.coefSDB
  ####################################################BLB and SDB begin
  
  ####################################################BLB and SDB begin
  fit.BLB08 <- BLB.glm(x=x, y=y, S=S.BLB, r=r, b=b08, beta0=beta, alpha=0.05)
  CI.BLB08[i, ] <- fit.BLB08$coef.ULBLB[2, ] - fit.BLB08$coef.ULBLB[1, ]
  p.BLB08[i, ]  <- fit.BLB08$p.coefBLB
  
  fit.SDB08 <- SDB.glm(x=x, y=y, S=S.SDB, b=b08, beta0=beta, alpha=0.05)
  CI.SDB08[i, ] <- fit.SDB08$coef.ULSDB[2, ] - fit.SDB08$coef.ULSDB[1, ]
  p.SDB08[i, ]  <- fit.SDB08$p.coefSDB
  ####################################################BLB and SDB begin

  print(i)
}

name.model <- c( "DACEL50", "DACEL100", "DACEL150","mvc", "mmse", 
                 "BLB06", "BLB08", "SDB06","SDB08", "boot", 
                 "CIDACEL50", "CIDACEL100", "CIDACEL150", "CImvc", "CImmse",
                 "CIBLB06", "CIBLB08", "CISDB06", "CISDB08", "CIboot")

result.p <- matrix(aa, length(name.model), p)

rownames(result.p) <- name.model 


####emprical size
result.p["DACEL50",1:p] <- apply(p.mog50, 2, function(x)(mean(I(x <= 0.05))))
result.p["DACEL100",1:p] <- apply(p.mog100, 2, function(x)(mean(I(x <= 0.05))))
result.p["DACEL150",1:p] <- apply(p.mog150, 2, function(x)(mean(I(x <= 0.05))))

result.p["mmse",1:p]  <- apply(p.mmse, 2, function(x)(mean(I(x <= 0.025))))
result.p["mvc",1:p] <- apply(p.mvc, 2, function(x)(mean(I(x <= 0.025))))

result.p["BLB06",1:p] <- apply(p.BLB06, 2, function(x)(mean(I(x <= 0.05))))
result.p["BLB08",1:p] <- apply(p.BLB08, 2, function(x)(mean(I(x <= 0.05))))

result.p["SDB06",1:p] <- apply(p.SDB06, 2, function(x)(mean(I(x <= 0.05))))
result.p["SDB08",1:p] <- apply(p.SDB08, 2, function(x)(mean(I(x <= 0.05))))

result.p["boot",1:p] <- apply(p.boot, 2, function(x)(mean(I(x <= 0.05))))

###################CI length


result.p["CIDACEL50", 1:p] <- apply(CI.mog50, 2, mean) 
result.p["CIDACEL100", 1:p] <- apply(CI.mog100, 2, mean) 
result.p["CIDACEL150", 1:p] <- apply(CI.mog150, 2, mean) 

result.p["CImvc", 1:p] <- apply(CI.mvc, 2, mean)
result.p["CImmse", 1:p] <- apply(CI.mmse, 2, mean)
result.p["CIBLB06", 1:p] <- apply(CI.BLB06, 2, mean)
result.p["CIBLB08", 1:p] <- apply(CI.BLB08, 2, mean)

result.p["CISDB06", 1:p] <- apply(CI.SDB06, 2, mean)
result.p["CISDB08", 1:p] <- apply(CI.SDB08, 2, mean)

result.p["CIboot", 1:p] <- apply(CI.boot, 2, mean)



print(case)

result.p





