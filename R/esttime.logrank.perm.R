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

esttime.logrank.perm <- function(N, k, no_cores=1, strat=FALSE, tform='secs', verbose=TRUE) {
  # Error checks
  check1 <- tform %in% c('secs','mins','hours')
  if (!(check1)) {
    print('Choose a valid tform: secs, mins, or hours')
    stopifnot(check1)
  }
  
  # ----- (2) Functions ------ # 
  
  lr.chi <- function(x,So,rho) {
    fit <- survival:::survdiff.fit(y=So,x=x,rho=rho)
    otmp <- fit$observed
    etmp <- fit$expected
    stat <- sum((otmp-etmp)^2)/sum(diag(fit$var))
    return(stat)
  }
  
  # (ii) Assign a vector ones to the locations of an index
  idx2ones <- function(idx,N) {
    zz <- rep(0,N)
    zz[idx] <- 1
    return(zz)
  }
  
  # ------- (2) Set-up ------- #
  
  # Load in the microbenchmark package
  # library(microbenchmark,warn.conflicts = F,quietly = T,verbose = F)
  library(rbenchmark,warn.conflicts = F,quietly = T,verbose = F)
  
  # Calculate the number of combinations
  ncombn <- choose(N,k)
  
  # Random survival data
  Tobs <- rexp(N)
  is.event <- rep(1,N)
  So <- matrix(c(Tobs,is.event),ncol=2)
  x <- rep(0,N)
  x[sample(1:N,k)] <- 1
  rho <- 0
  
  # ----- (3) Use rbenchmark ------ # 
    
  nops.max <- 1e4
  # Use benchmark function for sequence of up to 10K calculations (or lower)
  if (ncombn < nops.max) { 
    nops.max <- ncombn
  }
  # Number of sims
  nsim <- 3
  
  # Holders
  time.store <- rep(NA,nsim)
  for (q in 1:nsim) {
    start <- Sys.time()
    out <- rep(NA,length(nops.max))
    e <- 0
    h <- k
    a <- 1:k
    y <- seq(N)
    nmmp1 <- N - k + 1
    mp1 <- k + 1
    for (i in 2:nops.max) {
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
      out[i] <- lr.chi(x=idx2ones(y[a],N),So,rho)
      # if (mod(i,1000)==0) print(i)
    }
    out[1] <- lr.chi(x=idx2ones(1:k,N),So,rho)
    end <- Sys.time()
    qq <- as.numeric(end - start)
    time.store[q] <- qq
  }
  
  # Time for 10K (nops.max) operations
  est.secs <- (ncombn/nops.max) * mean(time.store)
  
  # Use our rough parallelization efficiency:
  nc <- ifelse(no_cores > 6,6,no_cores)
  para.skedge <- list('1'=1.0000000,'2'=0.6041985,'3'=0.3852139,
       '4'=0.3716457,'5'=0.2728588,'6'=0.2769552)
  para.gain <- para.skedge[[as.character(nc)]]
  # Update
  pred.secs.parallel <- est.secs * para.gain
  # Total
  total.secs <- pred.secs.parallel
  
  # Get the cost from the strat for the different cores
  if (strat) {
    strat.skedge <- list('1'=2.081311,'2'=1.858964,'3'=1.819308,
                         '4'=1.737949,'5'=1.706770,'6'=1.683083)
    strat.cost <- strat.skedge[[as.character(nc)]]
    total.secs <- total.secs * strat.cost  
  }
  # Combine
  # Format
  if (tform=='secs') {
    tt.ret <- total.secs
    print(sprintf('Estimated runtime: %0.1f seconds',tt.ret))
  }
  if (tform=='mins') {
    tt.ret <- total.secs/60
    print(sprintf('Estimated runtime: %0.1f minutes',tt.ret))
  }
  if (tform=='hours') {
    tt.ret <- (total.secs/60)/60
    print(sprintf('Estimated runtime: %0.1f hours',tt.ret))
  }
  # Return
  return(tt.ret)
}

