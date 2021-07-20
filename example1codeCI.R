

################ our proposed method 
lm.moel <- function(x, y, beta0=NULL,  K=NULL, type="consecutive", Intercept=FALSE, alph=0.05){#type = C("random", "consecutive")
  n <- length(y)
  p <- dim(x)[2]
  fit.cv <- cvTools::cvFolds(n, K = K, type = type) ##split data into K block 
  
  if(Intercept==FALSE){
    ############################################################## no intercept begin
    if(is.null(beta0)==TRUE){
      beta0 <- rep(0, p)
    }
    
    ### DAC (split and conquer) begin
    coef.full <- matrix(0, K, p)
    coef.test <- matrix(0, K, p)
    for(k in 1: K){
      index.k <- fit.cv$subsets[fit.cv$which==k]
      nk <- length(index.k)
      coef.full[k, ] <- lm(y[index.k]~x[index.k, ] - 1)$coefficients
      coef.test[k, ] <- sqrt(nk) * (coef.full[k, ] - beta0)
    }
    coef.mox <- Rfast::colmeans(coef.full)
    ### DAC edn
    
    ######find UL begin 
    a <- qchisq(p= 1- alph, df=1)
    myfun <- function(theta, x) {
      el.test(x, mu=theta)
    }
    
    coefel.UL <- matrix(0, p, 2)
    colnames(coefel.UL) <- c("Low","Up")
    for(j in 1:p){
      elUL <- findUL(step=0.001, fun=myfun, MLE= mean(coef.full[, j]), level=a, x=coef.full[, j])
      coefel.UL[j, 1] <- elUL$Low
      coefel.UL[j, 2] <- elUL$Up
    }
    lengthel.UL  <- coefel.UL[, 2]- coefel.UL[, 1]
    ######find UL end
    
    #####################DAC-EL begin 
    p.coefel <- NULL
    #p.coefel1 <- NULL # the same as  p.coefel
    for(j in 1:p){
      #p.coefel1[j] <- emplik::el.test(x=coef.full[, j], mu=beta0[j])$Pval 
      p.coefel[j] <- emplik::el.test(x=coef.test[, j], mu=0)$Pval
    }
    
    #p.allel  <- emplik::el.test(x=coef.full, mu=beta0)$Pval ## the same as sqrt(n) *coef.test
    p.allel  <- emplik::el.test(x=coef.test, mu=rep(0, p))$Pval
    ################## DAC-EL end
  
    
    names(coef.mox) <- colnames(x)
    names(p.coefel) <- colnames(x)
    ############################################################## no intercept begin
  }else{
    ######TRUE begin
    ### 
    if(is.null(beta0)==TRUE){
      beta0 <- rep(0, (p +1))
    }
    ### 
    coef.full <- matrix(0, K, (p+1))
    coef.test <- matrix(0, K, (p+1))
    for(k in 1: K){
      index.k <- fit.cv$subsets[fit.cv$which==k]
      nk <- length(index.k)
      coef.full[k, ] <- lm(y[index.k]~x[index.k, ])$coefficients
      coef.test[k, ] <- sqrt(nk) * (coef.full[k, ] - beta0)
    }
    coef.mox <- Rfast::colmeans( coef.full)
    
    
    ######find UL begin  
    a <- qchisq(p= 1- alph, df=1)
    myfun <- function(theta, x) {
      el.test(x, mu=theta)
    }
    
    coefel.UL <- matrix(0, (p+1), 2)
    colnames(coefel.UL) <- c("Low","Up")
    for(j in 1:(p+1)){
      elUL <- findUL(step=0.001, fun=myfun, MLE= mean(coef.full[, j]), level=a, x=coef.full[, j])
      coefel.UL[j, 1] <- elUL$Low
      coefel.UL[j, 2] <- elUL$Up
    }
    
    lengthel.UL  <- coefel.UL[, 2]- coefel.UL[, 1]
    ######find UL end
    
    #####################DAC-EL begin 
    p.coefel <- NULL
    for(j in 1:(p+1)){
      p.coefel[j] <- emplik::el.test(x=coef.test[, j], mu=0)$Pval
    }
    
    p.allel  <- emplik::el.test(x=coef.test[,-1], mu=rep(0, p))$Pval
    ##################DAC-EL end
    
    names(coef.mox) <- c("intercept", colnames(x))
    names(p.coefel) <- c("intercept", colnames(x))
  }

  return(list(coef=coef.mox,K=K, p.coefel=p.coefel, p.allel=p.allel, coefel.UL=coefel.UL,lengthel.UL=lengthel.UL))
}


################################BLB SDB
SDB.lm <- function(x, y, S=NULL , b=NULL, beta0=NULL, alpha=0.05){
  n <- dim(x)[1]
  p <- dim(x)[2]
  if( is.null(S)==TRUE ) S=200
  if( is.null(b)==TRUE ) b <- ceiling(n^0.6)
  if( is.null(beta0)==TRUE ) beta0 <- rep(0,p)
  
  resamples <- stats::rmultinom(S, n, rep(1/b, b))
  res <- lapply(1:S, function(ii) {
    index.S <- sample(1:n, b,  replace = TRUE )
    x.sub <- x[index.S,]
    y.sub <- y[index.S]
    weights <- resamples[, ii] / n
    coef(lm(y.sub~x.sub-1, weights=weights))
  })
  
  coef.SDB <- data.frame(do.call(rbind, res))
  
  #########begin P.VALUE
  p.coef  <- NULL
  
  alpha.L <- alpha / 2
  alpha.U <- 1 - alpha/2
  coef.UL <- apply(coef.SDB, 2, quantile,  probs = c(alpha.L, alpha.U))
  
  for(j in 1:p){
    p.coef[j] <- I( beta0[j] >=  coef.UL[1, j]  ) * I(beta0[j] <= coef.UL[2, j])
  }
  return(list(p.coefSDB=p.coef,  coef.ULSDB= coef.UL))
}


################################BLB begin
BLB.lm <- function(x, y, S=NULL, r=100, b=NULL, beta0=NULL, alpha=0.05){
  n <- dim(x)[1]
  p <- dim(x)[2]
  if( is.null(S)==TRUE ) S=50
  if( is.null(b)==TRUE ) b <- ceiling(n^0.6)
  if( is.null(beta0)==TRUE ) beta0 <- rep(0,p)
  
  
  coef.BLB <- matrix(0, S, p)
  
  ##########################BLB begin
  ################################begin S loop
  for(i in 1:S){
    index.S <- sample(1:n, b,  replace = TRUE )
    x.sub <- x[index.S,]
    y.sub <- y[index.S]
    
    
    ####begin r loop
    resamples <- stats::rmultinom(r, n, rep(1/b, b))
    res <- lapply(1:r, function(ii) {
      weights <- resamples[, ii] / n
      coef(lm(y.sub~x.sub-1, weights= weights))
    })
    
    b.wr <- data.frame(do.call(rbind, res))
    coef.BLB[i, ] <- apply(b.wr, 2, mean)   
  }
  ################################begin S end
 
  #########begin P.VALUE
  p.coefBLB  <- NULL
  
  alpha.L <- alpha / 2
  alpha.U <- 1 - alpha/2
  coef.ULBLB <- apply(coef.BLB, 2, quantile,  probs = c(alpha.L, alpha.U))
  
  for(j in 1:p){
    p.coefBLB[j] <- I( beta0[j] >=  coef.ULBLB[1, j]  ) * I(beta0[j] <= coef.ULBLB[2, j])
  }
  #########end BLB
  
  return(list( p.coefBLB=p.coefBLB, coef.ULBLB=coef.ULBLB))
}
#####################################################################BLB end

###############################################Bootstrap
lm.boot <- function(x, y, B=100, beta0=NULL, alpha=0.05){
  p <- dim(x)[2]
  if( is.null(beta0)==TRUE ) beta0 <- rep(0,p)
  
  res.full <- lapply(1:B, function(ii) {
    index.S <- sample(1:n, n,  replace = TRUE )
    x.sub <- x[index.S,]
    y.sub <- y[index.S]
    coef(lm(y.sub~x.sub-1))
  })
  
  coef.B <- data.frame(do.call(rbind, res.full))
  
  #########begin P.VALUE
  p.coef <- NULL
  
  alpha.L <- alpha / 2
  alpha.U <- 1 - alpha/2
  coef.UL <- apply(coef.B, 2, quantile,  probs = c(alpha.L, alpha.U))
  
  for(j in 1:p){
    p.coef[j] <- I( beta0[j] >=  coef.UL[1, j]  ) * I(beta0[j] <= coef.UL[2, j])
    
  }

  return(list( p.coef=p.coef, coef.UL=coef.UL))
}


