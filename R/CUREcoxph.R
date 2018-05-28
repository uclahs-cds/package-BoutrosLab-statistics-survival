CUREcoxph <- function(survform,cureform,So,data,bootstrap=NULL,verbose=TRUE,subset=NULL,tol=1e-3,max.iter=50) {
  # --- Error checking --- #
  # (i) survform
  if (any(class(survform) %in% c('character','formula')) ) {
    old.survform <- survform
    survform <- as.formula(paste('So',as.character(as.formula(survform))[2],sep='~'))
  } else { print('Error! survform needs to be a character or formula class')}
  # (ii) cureform
  if (any(class(cureform) %in% c('character','formula')) ) {
    cureform <- as.formula(cureform)
  } else { print('Error! cureform needs to be a character or formula class')}
  # Initialize the model
  em.old <- CUREcoxph.EM(NULL,survform,cureform,So,data,subset)
  # Initialize loop
  old.coefs <- em.old$coefs
  diff <- 1
  j=0
  while ((diff > tol) & j <= max.iter) {
    j=j+1
    em.new <- CUREcoxph.EM(em.old$params,survform,cureform,So,em.old$data,subset)
    new.coefs <- em.new$coefs
    diff <- sum(abs(new.coefs - old.coefs)^2)
    old.coefs <- new.coefs
    em.old <- em.new
    if (verbose & j==max.iter) print('Maximum iteration hit')
  }
  # Assign the coefficients
  coef.final <- em.new$coefs
  idx.int <- which(names(coef.final) == '(Intercept)')
  bb.final <- coef.final[1:(idx.int-1)]
  gg.final <- coef.final[idx.int:length(coef.final)]
  # Check for bootsrapping
  if (is.null(bootstrap)) {
    lst.final <- list(beta=bb.final,gamma=gg.final)
  } else {
    em.bs <- CUREcoxphBS(old.survform,cureform,So,data,bootstrap,verbose,subset,tol,max.iter)
    ci.beta <- rbind(bb.final,em.bs$beta)
    rownames(ci.beta)[1] <- 'point.est'
    lst.final <- list(beta=ci.beta,gamma=gg.final)
  }
  return(lst.final)
}

# Inverse logistic link
expit <- function(x) {
  return(1/(1+exp(-x)))
}

################################
# --- FUNCTION FOR EM STEP --- #

CUREcoxph.EM <- function(params,survform,cureform,So,dd,subset,tol=1e-3,mdl=F) {
  cn <- colnames(dd)
  NN <- nrow(dd)
  # Initialize columns is there are none
  if ( sum(cn %in% c('y','pw'))!=2 ) {
    dd <- data.frame(y=1,pw=1,dd)
    idx.cens <- which(So[,2]==0)
    idx.event <- which(So[,2]==1)
    # Assign pw on a uniform scale
    dd$pw[idx.cens] <- rev(seq_along(idx.cens)/length(idx.cens))
    # If greater than last event time, put to zero
    dd$y[which(So[,1] > max(So[,1][So[,2]==1]))] <- 0
  }
  if (is.null(params)) {
    lst.params <- list(exponential=NULL,weibull=NULL) #gompertz=NULL
  } else {
    lst.params <- params
  }
  # Formula for glm model
  pw.cureform <- as.formula(paste('pw',as.character(cureform)[2],sep='~'))
  # - (i) CoxPH model - #
  mdl.cox <- coxph(formula=survform,weights=pw,data=dd)
  # Re-scale to proper level (coxph rescale around etabar)
  bhat <- mdl.cox$coefficients
  eta.cox <- mdl.cox$linear.predictors + sum(mdl.cox$means * bhat)
  # - (ii) Parametric model parameters - #
  fit.exp <- survPH(params=lst.params$exponential,dist='exponential',tt=So[,1],delta=So[,2],eta=eta.cox,w=dd$pw,tol=tol,subset=subset)
  fit.wei <- survPH(params=lst.params$weibull,dist='weibull',tt=So[,1],delta=So[,2],eta=eta.cox,w=dd$pw,tol=tol,subset=subset)
  lst.params <- list(exponential=fit.exp$params,weibull=fit.wei$params) 
  # - (iii) Predicted survival probabilities - #
  dd$psurv <- pmin(fit.exp$sdf$shat,fit.wei$sdf$shat)
  # - (iv) Cure model - #
  mdl.cure <- glm(formula=pw.cureform,family=quasibinomial,data=dd)
  ghat <- coef(mdl.cure)
  # Predicted probabilities
  dd$pprob <- predict(mdl.cure,type='response')
  temp.pred <- predict(mdl.cure,type='link',se.fit=TRUE)
  dd$pprob <- expit(temp.pred$fit)
  # Estimate of probability weights
  dd$pw <- with(dd,So[,2] + (1-So[,2])*(pprob * psurv)/(1-pprob + pprob*psurv))
  dd$pw <- ifelse(dd$pw < 1e-6,1e-6,dd$pw)
  # Update y
  dd$y <- ifelse(dd$pw > 0.5, 1, 0)
  # - (v) Return the data and coefficients - #
  ret.lst <- list(data=dd,coefs=c(bhat,ghat),params=lst.params)
  # Of the models
  if (mdl) {
    ret.lst <- list(cox=mdl.cox,cure=mdl.cure,
                    data=dd[,c('y','pw','psurv','pprob')],
                    params=lst.params)
  }
  return(ret.lst)
}

###################################
# --- BOOTSTRAP FUNCTIONALITY --- #

CUREcoxphBS <- function(survform,cureform,So,data,bootstrap=499,verbose,subset,tol,max.iter) {
  lst.bb <- list()
  lst.gg <- list()
  idx.cens <- which(So[,2]==0)
  idx.event <- which(So[,2]==1)
  bs.idx <- floor(seq(1,bootstrap,length.out = 11))[-c(1,11)]
  for (b in 1:bootstrap) {
    if (any(bs.idx == b) & verbose) {
      print(sprintf('Bootstrap progress: %0.0f%% done',which(bs.idx == b)*10))
    }
    set.seed(b)
    # Block bootstrap
    ridx <- c(sample(idx.cens,replace = T),sample(idx.event,replace = T))
    # ridx <- sample(1:nrow(data),replace = T)
    temp.mdl <- CUREcoxph(survform,cureform,So[ridx],data[ridx,],bootstrap=NULL,verbose,subset,tol,max.iter)
    lst.bb[[b]] <- temp.mdl$beta
    lst.gg[[b]] <- temp.mdl$gamma
  }
  # Get the 95% CI using the quantile approach
  level <- 0.05
  ci.bb <- apply(do.call('rbind',lst.bb),2,quantile,probs=c(level/2,1-level/2))
  # ci.gg <- apply(do.call('rbind',lst.gg),2,quantile,probs=c(level/2,1-level/2))
  # Return in list
  lst.ret <- list(beta=ci.bb) #,gamma=ci.gg
  return(lst.ret)
}


###############################################
# --- (Negative) Log-likelihood functions --- #

nLL.PH <- function(dist,params,delta,tt,eta=NULL,w=NULL) {
  if (!any(dist %in% c('exponential','weibull','gompertz'))) {
    print('Error! Pick a valid distribution'); break
  }
  # Initialize value if null
  if (is.null(eta)) { eta <- rep(0,length(tt)) }
  if (is.null(w)) { w <- rep(1,length(tt)) }
  # Choose one of the distributions
  if (dist=='exponential') {
    gam <- params[1]
    bet <- 1#params[2]
    nLL <- -sum( w*( delta*(bet*eta) + delta*gam - exp(gam + bet*eta)*tt ) )
  }
  if (dist=='weibull') {
    if (params[2] <= 0) { print('Error! Weibull distribution shape parameter > 0') }
    gam <- params[1]
    alph <- params[2]
    bet <- 1#params[3]
    nLL <- -sum( w*( delta*(log(alph)+gam+(bet*eta)+(alph-1)*log(tt)) - exp(gam + bet*eta)*tt^alph ) )
  }
  if (dist=='gompertz') {
    if (length(params)!=2) { print('Error! Gompertz distribution should have two parameters') }
    if (params[1] <= 0) { print('Error! Gompertz scale parameter should be > 0')}
    scale <- params[1]
    shape <- params[2]
    nLL <- -sum( w*( delta*(log(scale)+eta+shape*tt) + (scale*exp(eta))/shape*(1-exp(shape*tt)) ) )
  }
  return(nLL)
}

###########################################################
# --- GRADIENT OF (Negative) Log-likelihood functions --- #

dnLL.PH <- function(dist,params,delta,tt,eta=NULL,w=NULL) {
  if (!any(dist %in% c('exponential','weibull','gompertz'))) {
    print('Error! Pick a valid distribution'); break
  }
  # Initialize value if null
  if (is.null(eta)) { eta <- rep(0,length(tt)) }
  if (is.null(w)) { w <- rep(1,length(tt)) }
  if (dist=='exponential') {
    gam <- params[1]
    bet <- params[2]
    dl.dgam <- -sum( w*( delta - exp(gam + bet*eta)*tt ) )
    dl.dbet <- -sum( w*( delta*eta - exp(gam + bet*eta)*eta*tt ) )
    dl.params <- matrix(c(dl.dgam),ncol=1) #dl.dbet
  }
  if (dist=='weibull') {
    gam <- params[1]
    alph <- params[2]
    bet <- 1 #params[3]
    dl.dgam <- -sum( w*( delta - exp(gam + bet*eta)*tt^alph ) )
    dl.dalph <- -sum( w*( delta*(1/alph + log(tt)) - exp(gam + bet*eta)*tt^alph*log(tt) ) )
    # dl.dbet <- -sum( w*( delta*eta - exp(gam + bet*eta)*eta*tt^alph ) )
    dl.params <- matrix(c(dl.dgam,dl.dalph),ncol=1) #dl.dbet
  }
  if (dist=='gompertz') {
    scale <- params[1]
    shape <- params[2]
    dl.dscale <- -sum( w*( delta/scale - exp(eta)/shape*(exp(shape*tt)-1) ) )
    dl.dshape <- -sum( w*( delta*tt - scale*exp(eta)*( (tt*exp(shape*tt))/shape - (exp(shape*tt)-1)/shape^2 ) ) )
    dl.params <- matrix(c(dl.dscale,dl.dshape),ncol=1)
  }
  # Return
  return(dl.params)
}

############################################################
# --- GRADIENT DESCENT WITH BARZILAI/-BORWEIN AND BTLS --- #

gd.survdist.PH <- function(params=NULL,dist,tt,delta,eta=NULL,w=NULL,tol=1e-5) {
  if (is.null(eta)) { eta <- rep(0,length(tt)) }
  if (is.null(w)) { w <- rep(1,length(tt)) }
  if (dist == 'exponential') {
    params <- log(sum(w * delta) / sum( w * exp(eta) * tt ))
    return(params)
    # if (is.null(params)) {
    #   gam.init <- log(sum(w * delta) / sum( w * exp(eta) * tt ))
    #   bet.init <- 1
    #   params <- c(gam.init,bet.init)
    # }
    # params.check <- function(pp) { return(pp) } # No constraints on exponential
  }
  if (dist=='weibull') {
    if (is.null(params)) {
      gam.init <- log(sum(w * delta) / sum( w * exp(eta) * tt ))
      alph.init <- 1
      # bet.init <- 1
      params.old <- c(gam.init,alph.init) #,bet.init
    }
    params.check <- function(pp) {
      pp[2] <- ifelse(pp[2]<0,1e-20,pp[2])
      return(pp)
    }
  }
  if (dist=='gompertz') {
    if (is.null(params)) {
    }
    params.check <- function(pp) { 
      pp[1] <- ifelse(pp[1]<0,1e-20,pp[1])
      return(pp)
    }
  }
  if (!is.null(params)) {
    params.old <- params
  }
  # - Barzilai-Borwein initialization - #
  grad.old <- dnLL.PH(dist,params.old,delta,tt,eta,w)
  LL.old <- nLL.PH(dist,params.old,delta,tt,eta,w)
  safe.step <- min(abs(params.old/grad.old))/5
  params.new <- params.check(params.old - safe.step*grad.old)
  while (nLL.PH(dist,params.new,delta,tt,eta,w) > LL.old - (safe.step/2)*sum(grad.old^2)) {
    safe.step <- 0.8 * safe.step
    params.new <- params.check(params.old - safe.step*grad.old)
  }
  LL.new <- nLL.PH(dist,params.new,delta,tt,eta,w)
  grad.new <- dnLL.PH(dist,params.new,delta,tt,eta,w)
  sk <- params.new - params.old
  yk <- grad.new - grad.old
  step1 <- sum(sk^2)/sum(sk*yk)
  step2 <- sum(sk*yk)/sum(yk*yk)
  grad.old <- grad.new
  params.old <- params.new
  # - Gradient descent - #
  diff <- 1
  j=0
  while ((diff > tol) & j <= 20) {
    j=j+1
    LL.old <- nLL.PH(dist,params.old,delta,tt,eta,w)
    safe.step <- max(step1,step2)
    params.new <- params.check(params.old - safe.step*grad.old)
    while (nLL.PH(dist,params.new,delta,tt,eta,w) > LL.old - (safe.step/2)*sum(grad.old^2)) {
      safe.step <- 0.8 * safe.step
      params.new <- params.check(params.old - safe.step*grad.old)
    }
    grad.new <- dnLL.PH(dist,params.new,delta,tt,eta,w)
    sk <- params.new - params.old
    yk <- grad.new - grad.old
    step1 <- sum(sk^2)/sum(sk*yk)
    step2 <- sum(sk*yk)/sum(yk*yk)
    diff <- sum(abs(params.new - params.old))
    params.old <- params.new
    grad.old <- grad.new
  }
  # Return
  ret.list <- as.vector(params.new)
  return(ret.list)
}

############################################################################
# --- FUNCTION TO CALCULATE SURVIVAL PROBABILITES FOR EACH OBSERVATION --- #

# params=gd.fit;
SfunPH <- function(params,dist,tt,eta=NULL,w=NULL) {
  # Iniliatize if NULL
  if (is.null(eta)) { eta <- rep(0,length(tt)) }
  if (is.null(w)) { w <- rep(1,length(tt)) }
  # Calculate
  if (dist == 'exponential') {
    gam <- params[1]
    bet <- 1 #params[2]
    shat <- exp(-exp(gam + bet*eta)*tt)
    hhat <- exp(gam + bet*eta)
    tmu <- exp(-(gam + bet*eta))
    tmed <- tmu * log(2)
  }
  if (dist == 'weibull') {
    gam <- params[1]
    alph <- params[2]
    bet <- 1 #params[2]
    hhat <- alph*exp(gam+bet*eta)*tt^(alph-1)
    shat <- exp(-exp(gam+bet*eta)*tt^alph)
    tmu <- exp(-(gam+bet*eta)/alph) * gamma(1+1/alph)
    tmed <- exp(-(gam+bet*eta)/alph) * (log(2))^(1/alph)
  }
  if (dist == 'gompertz') {
  }
  # Return
  dat <- data.frame(shat,hhat,tmu,tmed)
  return(dat)
}

################################################################
# --- MASTER WRAPPER FUNCTION: PARAMETERS, SHAT, TMED, TMU --- #

survPH <- function(params=NULL,dist,tt,delta,eta=NULL,w=NULL,subset=NULL,tol=1e-4) {
  if (is.null(subset)) { subset <- rep(1,length(tt)) }
  gd.fit <- gd.survdist.PH(params,dist,tt,delta,eta,w*subset,tol)
  gd.LL <- nLL.PH(dist,gd.fit,delta,tt,eta,w*subset)
  sdf <- SfunPH(gd.fit,dist,tt,eta,w*subset)
  ret.list <- list(params=gd.fit,LL=gd.LL,sdf=sdf)
  return(ret.list)
}
