

################our proposed method 
glm.moel <- function(x, y, beta0=NULL,  K=NULL, type="consecutive", family=binomial(link = "logit"), Intercept=FALSE, alph=0.05){#type = C("random", "consecutive")
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
      coef.full[k, ] <- glm(y[index.k]~x[index.k, ]-1, family=family)$coefficients
      coef.test[k, ] <- sqrt(nk) * (coef.full[k, ] - beta0)
    }
    coef.mox <- Rfast::colmeans(coef.full)
    ### DAC end
    

    
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
    
    #####################DACELbegin 
    p.coefel <- NULL
    for(j in 1:p){
      p.coefel[j] <- emplik::el.test(x=coef.test[, j], mu=0)$Pval
    }
    p.allel  <- emplik::el.test(x=coef.test, mu=rep(0, p))$Pval
    ################## DACEL end
    
    names(coef.mox) <- colnames(x)
    names(p.coefel) <- colnames(x)
    ############################################################## no intercept end
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
      coef.full[k, ] <- glm(y[index.k]~x[index.k, ], family=family)$coefficients
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
    
    #####################DAC begin 
    p.coefel <- NULL
    for(j in 1:(p+1)){
      p.coefel[j] <- emplik::el.test(x=coef.test[, j], mu=0)$Pval
    }
    p.allel  <- emplik::el.test(x=coef.test[,-1], mu=rep(0, p))$Pval
    ##################dac end
    
    names(coef.mox) <- c("intercept", colnames(x))
    names(p.coefel) <- c("intercept", colnames(x))
    
  }
  
  
  return(list(coef=coef.mox,p.coefel=p.coefel, p.allel=p.allel, coefel.UL=coefel.UL,lengthel.UL=lengthel.UL))
}

######################################wang method
#################come from the package
getMLE <- function(x, y, w) {
  d <- ncol(x)
  beta <- rep(0, d)
  loop  <- 1
  Loop  <- 100
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1 - 1 / (1 + exp(x %*% beta)))
    H <- t(x) %*% (pr * (1 - pr) * w * x)
    S <- colSums((y - pr) * w * x)
    tryCatch(
{shs <- NA
 shs <- solve(H, S) },
error=function(e){
  cat("\n ERROR :", loop, conditionMessage(e), "\n")})
if (is.na(shs[1])) {
  msg <- "Not converge"
  beta <- loop <- NA
  break
}
beta.new <- beta + shs
tlr  <- sum((beta.new - beta)^2)
beta  <- beta.new
if(tlr < 0.000001) {
  msg <- "Successful convergence"
  break
}
if (loop == Loop)
  warning("Maximum iteration reached")
loop  <- loop + 1
  }
list(par=beta, message=msg, iter=loop)
}

###############add V.simp output
twostep <- function(X, Y, beta0=NULL, r1, r2, method=c("mvc", "mmse", "uni"),alph=0.05) {
  call <- match.call()
  method <- match.arg(method)
  n <- length(Y)
  p <- dim(X)[2]
  ### 
  if(is.null(beta0)==TRUE){
    beta0 <- rep(0, p)
  }
  ### 
  
  ######################## uni begin 
  if (method == "uni") {
    idx.simp <- sample(1:n, r1, T)
    x.simp <- X[idx.simp,]
    y.simp <- Y[idx.simp]
    pinv.simp <- n
    fit.simp <- getMLE(x=x.simp, y=y.simp, w=pinv.simp)
    beta.simp <- fit.simp$par
    msg <- fit.simp$message
    if (fit.simp$message == "Successful convergence") {
      p.simp  <- 1 - 1 / (1 + exp(c(x.simp %*% beta.simp))) ### alg1 pi
      w.simp <- p.simp * (1 - p.simp) ### alg1 wi
      W.simp <- solve (t(x.simp) %*% (x.simp * (w.simp * n)) ) * r1 * n ### alg1 MX^{-1}  n can be canceled 
      Vc.simp <- t(x.simp) %*% (x.simp * (y.simp-p.simp)^2 * n^2) / r1^2 / n^2
      V.simp <- W.simp %*% Vc.simp %*% W.simp
      amse.simp <- sqrt(diag(V.simp))
      
      #########################
      z.each <- (beta.simp - beta0)/ amse.simp
      pp <- pnorm(q=z.each)
      p.each <- pmin(1-pp, pp)
      ##############CORRECT RATION BEGIN
      xbeta <- X %*%    beta.simp 
      
      pi <- 1 - 1 / (1 + exp(xbeta))
      haty <- I(pi >=0.5)
      y.min <- Y - haty
      core <- mean(y.min== 0)######correct
      TF <- sum(y.min== 1) /sum(Y== 1)  ######ture 1 but method 0
      FT <- sum(y.min== -1)/sum(Y== 0)######ture 0 but method 1
      ##############CORRECT RATION END
    }
    else{
      amse.simp = NA
      p.each <- NA
      core <- NA
      TF <- NA
      FT <- NA
    }
    
    
    
    return(list(par=beta.simp, Vc =V.simp, amse=amse.simp, p.each=p.each, msg=msg,core=core, TF=TF,FT=FT, method="uni"))
  }
  ######################## uni end
  
  #######################alogorithm 2 begin
  n1 <- sum(Y) ###the number S1
  n0 <- n - n1 ###the number S1
  PI.prop <- rep(1/(2*n0), n)
  PI.prop[Y==1] <- 1/(2*n1)
  idx.prop <- sample(1:n, r1, T, PI.prop)
  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]
  pinv.prop <- n ### the code can be canceled
  pinv.prop <- 1/PI.prop[idx.prop]
  fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
  beta.prop <- fit.prop$par
  if (is.na(beta.prop[1])) ###eoor
    return(list(opt=NA, msg="first stage not converge"))
  
  if (method == "mmse") {
    P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
    p.prop <- P.prop[idx.prop]
    w.prop <- p.prop * (1 - p.prop)
    W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop)) #########MX-1
    PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE / sum(PI.mMSE)
    ######step2 r1 + r2
    idx.mMSE <- sample(1:n, r2, T, PI.mMSE)
    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    pinv.mMSE <- c(1 / PI.mMSE[idx.mMSE], pinv.prop)
    fit.mMSE <- getMLE(x=x.mMSE, y=y.mMSE, w=pinv.mMSE)
    
    ru <- length(y.mMSE) ##r1+r2
    beta.mMSE <- fit.mMSE$par
    p.mMSE  <- 1 - 1 / (1 + exp(c(x.mMSE %*% beta.mMSE)))
    w.mMSE <- p.mMSE * (1 - p.mMSE)
    W.mMSE <- solve(t(x.mMSE) %*% (x.mMSE * (w.mMSE * pinv.mMSE))) * ru * n
    Vc.mMSE <- t(x.mMSE) %*% (x.mMSE * (y.mMSE-p.mMSE)^2 * pinv.mMSE^2) / ru^2 / n^2
    V.mMSE <- W.mMSE %*% Vc.mMSE %*% W.mMSE #####p.835 (15)
    
    amse <- sqrt(diag(V.mMSE))
    msg <- c(fit.prop$message, fit.mMSE$message)
    
    ###p.each begin
    #z.each <- beta.mMSE / amse
    z.each <- (beta.mMSE - beta0)/ amse
    pp <- pnorm(q=z.each)
    p.each <- pmin(1-pp, pp)
    ###p.each end
    
    ############FindUL begin 20190304
    length.UL <- abs(qnorm(p=alph/2))* amse
    ############FindUL end 20190304
    
    #####Tn for model
    V.mMSEslove <- limSolve::Solve(V.mMSE)
    Tn.mMSE <- t(beta.mMSE - beta0) %*%   V.mMSEslove %*% (beta.mMSE - beta0)
    p.all <- pchisq(Tn.mMSE, df=p)
    
    ##############CORRECT RATION BEGIN
    xbeta <- X %*%  beta.mMSE 
    
    pi <- 1 - 1 / (1 + exp(xbeta))
    haty <- I(pi >=0.5)
    y.min <- Y - haty
    core <- mean(y.min== 0)######correct
    TF <- sum(y.min== 1) /sum(Y== 1)  ######ture 1 but method 0
    FT <- sum(y.min== -1)/sum(Y== 0)######ture 0 but method 1
    ##############CORRECT RATION END
    
    return(list(par=beta.mMSE, Vc=V.mMSE, amse=amse, p.each=p.each, p.all=p.all,msg=msg,core=core, TF=TF,FT=FT,method="mmse", length.UL=length.UL))
  }
  
  if (method == "mvc") {
    P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
    PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
    PI.mVc <- PI.mVc / sum(PI.mVc)
    idx.mVc <- sample(1:n, r2, T, PI.mVc)
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(1 / PI.mVc[idx.mVc], pinv.prop)
    fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    
    ru <- length(y.mVc)
    beta.mVc <- fit.mVc$par
    p.mVc  <- 1 - 1 / (1 + exp(c(x.mVc %*% beta.mVc)))
    w.mVc <- p.mVc * (1 - p.mVc)
    W.mVc <- solve(t(x.mVc) %*% (x.mVc * (w.mVc * pinv.mVc))) * ru * n
    Vc.mVc <- t(x.mVc) %*% (x.mVc * (y.mVc-p.mVc)^2 * pinv.mVc^2) / ru^2 / n^2
    V.mVc <- W.mVc %*% Vc.mVc %*% W.mVc
    ## ## LCC
    ## PI.optU <- abs(Y - P.prop)
    ## PI.optU <- PI.optU / sum(PI.optU)
    ## idx.optU <- sample(1:n, r2, T, PI.optU)
    ## x.optU <- X[c(idx.optU),]
    ## y.optU <- Y[c(idx.optU)]
    ## pinv.optU <- 1
    ## fit.optU <- getMLE(x=x.optU, y=y.optU, w=pinv.optU)
    
    amse <- sqrt(diag(V.mVc))
    msg <- c(fit.prop$message, fit.mVc$message)
    ###p.each begin
    #z.each <- beta.mVc / amse
    z.each <- (beta.mVc - beta0)/ amse
    pp <- pnorm(q=z.each)
    p.each <- pmin(1-pp, pp)
    ###p.each end
    
    ############FindUL begin 20190304
    length.UL <- abs(qnorm(p=alph/2))* amse
    ############FindUL end 20190304
    
    
    #####Tn for model
    V.mVcslove <- limSolve::Solve(V.mVc)
    Tn.mVc <- t(beta.mVc - beta0) %*%   V.mVcslove %*% (beta.mVc - beta0)
    p.all <- pchisq(Tn.mVc, df=p)
    
    ##############CORRECT RATION BEGIN
    xbeta <- X %*%  beta.mVc 
    
    pi <- 1 - 1 / (1 + exp(xbeta))
    haty <- I(pi >=0.5)
    y.min <- Y - haty
    core <- mean(y.min== 0)######correct
    TF <- sum(y.min== 1) /sum(Y== 1)  ######ture 1 but method 0
    FT <- sum(y.min== -1)/sum(Y== 0)######ture 0 but method 1
    ##############CORRECT RATION END
    
    return(list(par=beta.mVc, Vc=V.mVc, amse=amse, p.each=p.each, p.all=p.all, msg=msg, core=core, TF=TF,FT=FT, method="mvc", length.UL=length.UL))
  }
}

################################BLB SDB
SDB.glm <- function(x, y, S=NULL , b=NULL, beta0=NULL, alpha=0.05){
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
    coef(glm(y.sub~x.sub-1, family = "binomial", weights=weights))
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
BLB.glm <- function(x, y, S=NULL, r=100, b=NULL, beta0=NULL, alpha=0.05){
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
      #coef(lm(y.sub~x.sub-1, weights= weights))
      coef(glm(y.sub~x.sub-1, family = "binomial", weights=weights))
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
glm.boot <- function(x, y, B=100, beta0=NULL, alpha=0.05){
  p <- dim(x)[2]
  if( is.null(beta0)==TRUE ) beta0 <- rep(0,p)
  
  res.full <- lapply(1:B, function(ii) {
    index.S <- sample(1:n, n,  replace = TRUE )
    x.sub <- x[index.S,]
    y.sub <- y[index.S]
    coef(glm(y.sub~x.sub-1, family = "binomial"))
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

