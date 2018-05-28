# The BoutrosLab.statistics.general package is copyright (c) 2017 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

logrank.perm <- function(So, x, rho=0, strat=NULL, no_cores=1, verbose=TRUE) {
  NSo <- nrow(So)
  pSo <- ncol(So)
  N <- length(x)
  k <- sum(x)
  # Checks
  check1 <- NSo==N
  check2 <- pSo==2
  if (!is.null(strat)) {
    check3 <- length(strat)==length(x)
  } else {
    check3 <- TRUE
  }
  
  # ----- Error checks ------ #
  if (!check1) {
    print('Error! The number of observations in So does not match x')
    stopifnot(check1)
  }
  if (!check2) {
    print('Error! So should have two columns!')
    stopifnot(check2)
  }
  if (!check3) {
    print('Error! Strata group should be the same length as x!')
    stopifnot(check3)
  }
  
  # ----- (1) Functions ------ # 
  
  # (ia) non-stratified log-rank stat
  lr.chi1 <- function(x,So,rho,strat) {
    fit <- survival:::survdiff.fit(y=So,x=x,rho=rho)
    otmp <- fit$observed
    etmp <- fit$expected
    stat <- sum((otmp-etmp)^2)/sum(diag(fit$var))
    return(stat)
  }
  # (ib) stratified log-rank stat
  lr.chi2 <- function(x,So,rho,strat) {
    fit <- survival:::survdiff.fit(y=So,x=x,rho=rho,strat=strat)
    otmp <- apply(fit$observed,1,sum)
    etmp <- apply(fit$expected,1,sum)
    stat <- sum((otmp-etmp)^2)/sum(diag(fit$var))
    return(stat)
  }
  # Assign the function to the call
  if (is.null(strat)) {
    lr.chi <- lr.chi1
  } else {
    lr.chi <- lr.chi2
  }
  
  # (ii) Assign a vector ones to the locations of an index
  idx2ones <- function(idx,N) {
    zz <- rep(0,N)
    zz[idx] <- 1
    return(zz)
  }
  
  # (iii) Combine functions
  idx2stat <- function(idx,N,statfun,So,rho,strat) {
    zz <- rep(0,N)
    zz[idx] <- 1
    return(statfun(x=zz,So=So,rho=rho,strat=strat))
  }
  
  # ----- (2) Run all permutations  ------ # 
  
  # Parallelization
  if (no_cores > 1) {
    if (verbose) print('Registering cores')
    library(doParallel,quietly = T,warn.conflicts = F,verbose = F)
    library(parallel,quietly = T,warn.conflicts = F,verbose = F)
    nc <- detectCores()
    check4 <- (nc-1) >= no_cores
    if (!check4) {
      print('Error! Use one less core that is available!')
      stopifnot(check4)
    }
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
  }
  
  # Record run time
  time.start <- Sys.time()
  
  # Make progress
  count <- choose(N, k)
  if (verbose) print(sprintf('Beginning execution: %i total combinations',count))
  if (count < 1e5 & verbose) print('Warning! <100K combinations will not benefit from multiple cores')
  pcs <- floor(count*seq(1,10)/10)[-10]
  pclabs <- paste(seq(1,9)*10,'% done',sep='')
  out <- rep(NA,length=count)
  
  # Holders
  e <- 0
  h <- k
  a <- 1:k
  y <- seq(N)
  nmmp1 <- N - k + 1
  mp1 <- k + 1
  
  # ---- COMBINATIONS FOR PARALLEL LOOP ---- #
  if (no_cores > 1) {
    if (verbose) print('Creating matrix of indices')
    matidx <- matrix(0,nrow=count,ncol=k)
    test <- for (i in 2:count) {
      # Prepare x
      if (e < N - h) {
        h <- 1
        e <- a[k]
        j <- 1
      } else {
        h <- h + 1
        e <- a[mp1 - h]
        j <- 1:h
      }
      a[k - h + j] <- e + j
      # Store idx
      matidx[i,] <- y[a]
    }
    matidx[1,] <- 1:k
    time.end1 <- Sys.time()
    # Now run the log-rank in paralelle
    if (verbose) print('Matrix completed, running statistic in parallel')
    out <- parApply(cl=cl,X=matidx,MARGIN=1,FUN=idx2stat,N=N,statfun=lr.chi,So=So,rho=rho,strat=strat)
    # End cluster
    stopCluster(cl)
    time.end2 <- Sys.time()
    # Total seconds
    sec.total <- as.numeric(difftime(time.end2,time.start,units='secs'))
    sec.mat <- as.numeric(difftime(time.end1,time.start,units='secs'))
    sec.para <- as.numeric(difftime(time.end2,time.end1,units='secs'))
    if (verbose) print(sprintf('Function took a total of %0.1f seconds, %0.1f (for mat) and %0.1f (for log-rank)',
                               sec.total,sec.mat,sec.para))
  }
  # ---- COMBINATIONS FOR NORMAL LOOP LOOP ---- #
  if (no_cores == 1) {
    out <- rep(NA,length=count)
    for (i in 2:count) {
      # Prepare x
      if (e < N - h) {
        h <- 1
        e <- a[k]
        j <- 1
      } else {
        h <- h + 1
        e <- a[mp1 - h]
        j <- 1:h
      }
      a[k - h + j] <- e + j
      # Get log-rank statistic
      out[i] <- lr.chi(x=idx2ones(y[a],N),So,rho,strat)
      if(any(pcs==i)) {
        if (verbose) print(pclabs[which(pcs==i)])
      }
    }
    out[1] <- lr.chi(x=idx2ones(1:k,N),So,rho,strat)
    # Total seconds
    time.end <- Sys.time()
    sec.combn <- as.numeric(difftime(time.end,time.start,units='secs'))
    if (verbose) print(sprintf('Function took %0.1f seconds',sec.combn))
  }
  # Remove any NAs (this happens when the censored observations fall on the x vector positions)
  out <- na.omit(out)
  
  # ----- (3) p-value and final info ------ # 
  
  # Calculate our baseline statistic
  baseline.stat <- lr.chi(x=x,So=So,rho=rho,strat=strat)
  baseline.pval <- pchisq(baseline.stat,1,lower.tail = F)
  
  # Number of statistics that are larger (w/ equality)
  ngreater <- sum(out > baseline.stat)
  pval <- ngreater/count
  if (pval==0) {
    pval <- 1/(count+1)
  }
  
  # Return
  dat <- data.frame(ntests=count,ngeq=ngreater,exact=pval,logrank=baseline.pval)
  return(dat)
}



