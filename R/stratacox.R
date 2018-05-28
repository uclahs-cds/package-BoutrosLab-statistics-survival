#####################################################################################
####### ------ PROXIMAL GRADIENT DESCENT FOR STRATIFIED COX MODELS ------- ##########

#####################################
### ----- HELPER FUNCTIONS ------ ###

# Sigmoid
expit <- function(x) { return( 1/(1+exp(-x)) ) }

##############################
### ----- LIBRARIES ------ ###

library(survival); library(data.table)

########################################################
####### ------ COX/SURVIVAL FUNCTIONS ------- ##########

# Create the Ymatrix (Note! Uses Breslow approximation)
Ymat.fun <- function(time) {
  N <- length(time)
  Ymat <- matrix(0,nrow=N,ncol=N)
  Ymat[upper.tri(Ymat)] <- 1
  diag(Ymat) <- 1
  # Get any duplicate time counts
  tdup <- unique(time[duplicated(time)])
  if (length(tdup)>=1) {
    for (tt in tdup) {
      idx.dup <- which(time==tt)
      Ymat[idx.dup,idx.dup] <- 1
    }
  }
  return(Ymat)
}

# Partial likelihood P matrix
Pfun <- function(eta,Ymat,tYmat) {
  haz <- as.vector(exp(eta))
  rsk <- as.vector(Ymat %*% haz)
  P <- outer(haz,rsk,FUN='/') * tYmat
  return(P)
}

# Partial likelihood
loss.partial <- function(eta,resp,Ymat) { 
  omega <- as.vector(exp(eta))
  Omega <- Ymat %*% omega
  nll <- -sum(resp * log(omega/Omega))
  return(nll)
}

# Grad for the different components
grad.nabla <- function(X,w,resp,P) {
  eta <- (X %*% w)
  Presp <- P %*% resp
  dll <- -(t(X) %*% (resp - Presp))
  return(as.vector(dll))
}

# Min without zeros
zmin <- function(zz) return(min(zz[zz!=0]))


# Function that inserts a column of NAs if a function is missing any
cn.standardize <- function(x,vars) {
  cn <- colnames(x)
  nc <- length(cn)
  # Vars to add
  v2a <- vars[which(!(vars %in% cn))]
  if (length(v2a) >= 1) {
    x <- cbind(matrix(rep(NA,nrow(x)*length(v2a)),nrow=nrow(x)),x)
    colnames(x)[1:length(v2a)] <- v2a
    cn <- c(v2a,cn)
    return(x[,match(vars,cn)])
  } else {
    return(x)
  }
}

# Main function
stratacox <- function(X,So,strat,lam,alpha,tau=0,target=1,w=NULL,
                          standardize=T,ord=T,max.iter=250,tol=1e-8,findlam=F,verbose=F) {
  # ---- (i) Preliminaries ----  #
  stopifnot(nrow(X) == nrow(So))

  # ---- (ii) Set up lists ----  #
  strat <- as.numeric(as.factor(strat))
  So <- data.table(strata=strat,So)
  # order by the death events
  if (ord) {
    if (verbose) { print('Ordering data') }
    X <- data.table(time=So$time,strata=strat,X)
    X <- X[order(strata,time)]
    So <- So[order(strata,time)]
    strat <- So$strata
    X <- as.matrix(X[,-(1:2)]) # Remove the time and strata
  }
  time <- So$time
  # Set up the Ymatrices and event
  ldelta <- split(So$status,strat)
  lYmat <- tapply(So$time,strat,Ymat.fun)
  ltYmat <- lapply(lYmat,t)

  # ---- (iii) Standardize X ----  #
  # Make X a matrix and standardize
  if (standardize) {
    if (verbose) { print('Standardizing data') }
    X <- scale(as.matrix(X))
    vecMu <- attr(X,'scaled:center')
    vecSd <- attr(X,'scaled:scale')
  } else {
    vecMu <- rep(0,ncol(X))
    vecSd <- rep(1,ncol(X))
  }

  # ---- (iv) S and Pi matrices ----  #
  N <- nrow(X)
  p <- ncol(X)
  ustrata <- unique(strat)
  D <- length(ustrata)
  pi <- tapply(strat,strat,length)/N
  pi[-target] <- (1-tau)*pi[-target]
  pi[target] <- 1 - sum(pi[-target])
  Pi <- matrix(rep(pi,each=p),ncol=D)
  # Get the missing lst
  lst.missing <- list()
  for (s in 1:D) {
    lst.missing[[s]] <- !is.na(X[which(strat == ustrata[s])[1],])
  }
  S <- t(do.call('rbind',lst.missing))
  rownames(S) <- NULL
  iotaD <- matrix(1,nrow=D)
  Omega <- (S * Pi) %*% iotaD
  OmegaInv <- as.vector(1/Omega)
  stopifnot(!any(OmegaInv==Inf))
  # Replace all NAs with zeros
  X[is.na(X)] <- 0

  # ---- (v) Gradient functions ----  #
  # Soft-thresholding
  Softfun <- function(w,thresh,int) {
      vec <- ifelse(w-thresh>0,pmax(w-thresh,0),pmin(w+thresh,0))
      if (int) { vec[1] <- w[1] }
      return(vec)
  }
  # Function for Barzilai Borwein Gradient
  bb.grad.fun <- function(step,alpha,lam,X,ow,strat,lYmat,ltYmat,ldelta) {
    N <- nrow(X)
    # Calculate with old weights
    ow.eta <- split(X %*% ow,strat)
    ow.res <- mapply(function(eta,yy,tyy,dd) dd -  (Pfun(eta,yy,tyy) %*% dd),ow.eta,lYmat,ltYmat,ldelta,SIMPLIFY = F)
    ow.nabla <- matrix(NA,nrow=p,ncol=D)
    for (d in 1:D) {
      ow.nabla[,d]  <- -(t(X[strat==d,]) %*% ow.res[[d]])
    }
    ow.grad <- (1/N) * OmegaInv * as.vector((ow.nabla * Pi) %*% iotaD) + (1-alpha)*lam*ow#(1/N) * (OmegaInv * ow.nabla) + (1-alpha)*lam*ow
    # Calculate with new weights
    z <- ow - step*ow.grad
    w <- Softfun(z,step*alpha*lam,int=F)
    w.eta <- split(X %*% w,strat)
    w.res <- mapply(function(eta,yy,tyy,dd) dd -  (Pfun(eta,yy,tyy) %*% dd),w.eta,lYmat,ltYmat,ldelta,SIMPLIFY = F)
    w.nabla <- matrix(NA,nrow=p,ncol=D)
    for (d in 1:D) {
      w.nabla[,d]  <- -(t(X[strat==d,]) %*% w.res[[d]])
    }
    w.grad <- (1/N) * OmegaInv * as.vector((w.nabla * Pi) %*% iotaD) + (1-alpha)*lam*ow
    # Difference
    diff <- sum(abs(w - ow))
    # Return
    ret.list <- list(w=w,ow=ow,w.eta=w.eta,w.grad=w.grad,ow.grad=ow.grad,diff=diff)
    return(ret.list)
  }

  # ---- (vi) Initialize Barzilai-Borwein ----  #
  # print('Barzilai borwein')

  grad.wrapper <- function(w,alpha,lam,X,strat,lYmat,ltYmat,ldelta,findlam) {
    p <- ncol(X)
    N <- nrow(X)
    if (is.null(w)) { w <- rep(0,p) }
    ow <- w
    if (verbose) { print('Initializing BB') }
    # Take a gradient step
    step <- 0.01 #zmin(abs(ow.grad))
    # Initialize gradiet
    bb.desc <- bb.grad.fun(step,alpha,lam,X,ow,strat,lYmat,ltYmat,ldelta)
    w <- bb.desc$w
    w.grad <- bb.desc$w.grad
    w.eta <- bb.desc$w.eta
    ow <- bb.desc$ow
    ow.grad <- bb.desc$ow.grad
    diff <- sum(abs(w - ow))

    if (findlam) { return( max(abs(bb.desc$ow.grad)) ) }

    # sum of log-likelihoods..
    w.loss <- sum(mapply(function(eta,dd,yy) loss.partial(eta,dd,yy),w.eta,ldelta,lYmat) * Pi[1,])/N +
      (1/2)*(1-alpha)*lam*sum(Omega * w^2) + alpha*lam*sum(Omega * abs(w))

    # ---- (vii) Begin gradient descent ----  #
    if (verbose) { print('Running gradient descent') }
    k=1; ll=c(w.loss,rep(NA,max.iter-1)); ss=c(step,rep(NA,max.iter-1))
    while (diff > tol & (k <= max.iter)) {
      k <- k+1
      if (verbose) { if( (k %%10)==0) {print(k)} }
      # Update Uk
      s <- w - ow
      g <- w.grad - ow.grad
      aSD <- sum(s^2)/sum(s*g)
      aMG <- sum(s*g)/sum(g^2)
      step <- max(min(c(aSD,aMG)),1e-4)
      # Descend
      ow <- w
      bb.desc <- bb.grad.fun(step,alpha,lam,X,ow,strat,lYmat,ltYmat,ldelta)
      w <- bb.desc$w
      w.grad <- bb.desc$w.grad
      w.eta <- bb.desc$w.eta
      ow <- bb.desc$ow
      ow.grad <- bb.desc$ow.grad

      w.loss <- sum(mapply(function(eta,dd,yy) loss.partial(eta,dd,yy),w.eta,ldelta,lYmat) * Pi[1,])/N +
        (1/2)*(1-alpha)*lam*sum(Omega * w^2) + alpha*lam*sum(Omega * abs(w))
      # Status
      diff <- sum(abs(w - ow))
      ll[k] <- w.loss
      ss[k] <- step

    }
    ws <- as.vector(w/vecSd)
    names(ws) <- colnames(X)
    ret.list <- list(w=ifelse(abs(ws)<10e-8,0,ws),ll=na.omit(ll),step=na.omit(ss),niter=k)
    return(ret.list)
  }

  if (findlam) {
    lamstar <- grad.wrapper(w,alpha,lam,X,strat,lYmat,ltYmat,ldelta,findlam)
    return(lamstar)
  }

  if (length(lam) == 1) {
    mdl <- grad.wrapper(w,alpha,lam,X,strat,lYmat,ltYmat,ldelta,findlam)
    return(mdl)
  } else {
    # Sort the lam
    lam <- sort(lam)
    w.old <- rep(0,p)
    w.list <- list()
    for (ll in seq_along(lam)) {
      if (verbose) { print(sprintf('lam: %0.3f',lam[ll])) }
      mdl <- grad.wrapper(w.old,alpha,lam[ll],X,strat,lYmat,ltYmat,ldelta,findlam)
      w.old <- mdl$w
      w.list[[ll]] <- w.old
    }
    w.mat <- do.call('cbind',w.list)
    w.mat <- data.frame(vars=rownames(w.mat),w.mat,row.names = NULL)
    colnames(w.mat)[-1] <- paste('lam',seq_along(lam),sep = '.')
    return(w.mat)
  }

}