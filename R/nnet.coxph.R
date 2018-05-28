#####################################
# ---------- LIBRARIES ------------ #

library(stringr);
library(survival);

#############################################
# ---------- SURVIVAL MATRICES ------------ #

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
};

# Partial likelihood P matrix
Pfun <- function(A,Ymat,tYmat) {
  haz <- as.vector(exp(A))
  rsk <- as.vector(Ymat %*% haz)
  P <- outer(haz,rsk,FUN='/') * tYmat
  return(P)
};

###################################################
# ---------- LOSS/ACCURACY FUNCTIONS ------------ #

conc.acc <- function(A, So) { 
  A <- as.vector(A)
  return(survConcordance(So ~ A)$concordance)
};

# Partial likelihood
partial.loss <- function(A,delta,Ymat) { 
  A <- as.vector(A)
  omega <- as.vector(exp(A))
  Omega <- Ymat %*% omega
  nll <- -sum(delta * (A  - log(Omega)))
  return(nll)
};

partial.backprop <- function(delta,Pmat) {
  deriv <- (delta - Pmat %*% delta)
  return(-matrix(deriv,nrow=1))
};


########################################
# -------- GRADIENT CHECKER ---------- #

# Function that returns the N-dimensional gradient checker
gradcheckfun.surv <- function(X,So,delta,Ymat,tYmat,strat,aterminal,dterminal,loss,weights,lam,leaky,ee=1e-7) {
  # Forward prop to get the cache
  cache <- forward.propagation.surv(X,weights,aterminal,leaky)
  namAL <- str_c('A',length(weights)/2) # weights have Wi bi
  # Update gradient with full batch (we don't want noisy sgd)
  K <- length(delta)
  tevents <- sum(unlist(delta))
  Pimat <- vector('list',K)
  for (k in seq(K)) {
    Pimat[[k]] <- Pfun(cache[[namAL]][strat == k],Ymat[[k]],tYmat[[k]])
  }
  full.grads <- backward.propagation.surv(X, delta, cache, Pimat, weights, 
                                          dterminal, leaky, keep.prob=NULL, lam)
  # Make a copy of the weights and gradients
  temp.weights1 <- temp.weights2 <- weights
  temp.nams <- names(weights)
  temp.dnams <- str_c('d',temp.nams)
  full.grads <- full.grads[temp.dnams]
  temp.grads <- lapply(full.grads,function(ll) ifelse(ll==0,ll,0))
  # Loop through each layer and nudge
  for (ll in temp.nams) {
    print(ll)
    dll <- str_c('d',ll)
    if (class(temp.weights1[[ll]])=='matrix') {
      temp.dims <- dim(temp.weights1[[ll]])
      for (r1 in 1:temp.dims[1]) {
        for (c1 in 1:temp.dims[2]) {
          print(sprintf('r1: %i, c1: %i',r1,c1))
          # Nudge each
          temp.weights1[[ll]][r1,c1] <- weights[[ll]][r1,c1] + ee
          temp.weights2[[ll]][r1,c1] <- weights[[ll]][r1,c1] - ee
          # Calculate the change in the loss function
          loss1 <- loss2 <- 0
          for (k in seq(K)) {
            loss1 <- loss1 + loss(forward.propagation.surv(X,temp.weights1,aterminal,leaky)[[namAL]][strat == k],delta[[k]],Ymat[[k]])
            loss2 <- loss2 + loss(forward.propagation.surv(X,temp.weights2,aterminal,leaky)[[namAL]][strat == k],delta[[k]],Ymat[[k]])
          }
          temp.gradapprox <- ((loss1 - loss2)/(2*ee)) / tevents
          temp.grads[[dll]][r1,c1] <- temp.gradapprox
          # Revert
          temp.weights1[[ll]][r1,c1] <- weights[[ll]][r1,c1]
          temp.weights2[[ll]][r1,c1] <- weights[[ll]][r1,c1]
        }
      }
    } else {
      temp.len <- length(temp.weights1[[ll]])
      for (v1 in seq(temp.len)) {
        # Nudge
        temp.weights1[[ll]][v1] <- weights[[ll]][v1] + ee
        temp.weights2[[ll]][v1] <- weights[[ll]][v1] - ee
        # Calculate the loss
        loss1 <- loss2 <- 0
        for (k in seq(K)) {
          loss1 <- loss1 + loss(forward.propagation.surv(X,temp.weights1,aterminal,leaky)[[namAL]][strat == k],delta[[k]],Ymat[[k]])
          loss2 <- loss2 + loss(forward.propagation.surv(X,temp.weights2,aterminal,leaky)[[namAL]][strat == k],delta[[k]],Ymat[[k]])
        }
        temp.gradapprox <- ((loss1 - loss2)/(2*ee)) / tevents
        temp.grads[[dll]][v1] <- temp.gradapprox
        # Reset
        temp.weights1[[ll]][v1] <- weights[[ll]][v1]
        temp.weights2[[ll]][v1] <- weights[[ll]][v1]
      }
    }
  }
  vec.full.grads <- unlist(full.grads)
  vec.temp.grads <- unlist(temp.grads)
  diff <- sqrt(sum((vec.full.grads-vec.temp.grads)^2))/(sqrt(sum(vec.full.grads^2))+sqrt(sum(vec.temp.grads^2)))
  return(diff)
};

###############################################
# -------- INDEX GENERATORS FOR CV ---------- #

# Function to do a block sample split of the survival data
idx.gen <- function(So,frac,ss=1) {
  set.seed(ss)
  som <- as.matrix(So)
  delta <- som[,2]
  n <- length(delta)
  idx.cens <- which(delta==0)
  n.cens <- length(idx.cens)
  idx.event <- which(delta==1)
  n.event <- length(idx.event)
  test.idx <- sample(idx.event,ceiling(n.event * frac)) # Sample from event
  test.idx <- c(test.idx,sample(idx.cens, ceiling(n * frac) - length(test.idx) )) # sample from cens
  train.idx <- setdiff(1:n,test.idx)
  return(list(train=train.idx,test=test.idx))
};

####################################################
# -------- CONCORDANCE FUNCTION WRAPPER ---------- #

concFun <- function(So,X,w) {
  cnX <- colnames(X)
  idxX <- which(cnX %in% names(w))
  if (length(idxX) == 0) {
    eta <- rep(0,nrow(X))
  } else {
    cnXhit <- cnX[idxX]
    wX <- w[match(cnXhit,names(w))]
    if (length(wX)==1) {
      eta <- as.vector(X[,idxX] * wX)
    } else {
      eta <- as.vector(X[,idxX] %*% wX)
    }
  }
  temp.conc <- survConcordance(So ~ eta)$concordance
  return(temp.conc)
};


################################################
# ---------- ACTIVATION FUNCTIONS ------------ #

relu <- function(x) {
  s <- pmax(x,0) 
  return(s)
};

lrelu <- function(x,l) {
  s <- pmax(x,0) + pmin(x,0)*l
  return(s)
};

linear.activation <- function(x) { return(x) };

################################################
# -------- PARAMETER INITIALIZATION ---------- #

# He initialization
initialize.weights.he <- function(layer.dims,ss,sig=1) {
  set.seed(ss)
  weights <- list()
  L <- length(layer.dims) - 1
  for (i in 1:L) {
    nl <- layer.dims[i+1]
    nl1 <- layer.dims[i]
    weights[[str_c('W',i)]] <- matrix(rnorm(nl * nl1,sd=sig),nrow=nl,ncol=nl1) / sqrt(2/nl1)
    weights[[str_c('b',i)]] <- rep(0,nl)
  }
  return(weights)
};

# Initialize the adam parameters with zeros
initialize.adam <- function(layer.dims) {
  v <- list()
  s <- list()
  L <- length(layer.dims) - 1
  for (i in 1:L) {
    nl <- layer.dims[i+1]
    nl1 <- layer.dims[i]
    v[[str_c('dW',i)]] <- matrix(0,nrow=nl,ncol=nl1)
    v[[str_c('db',i)]] <- rep(0,nl)
    s[[str_c('dW',i)]] <- matrix(0,nrow=nl,ncol=nl1)
    s[[str_c('db',i)]] <- rep(0,nl)
  }
  return(list(v=v,s=s))
};

###########################################
# -------- FORWARD PROPOGATION ---------- #

forward.propagation.surv <- function(X, parameters, aterminal, leaky=0, keep.prob=NULL) {
  # Get the number of layers
  L <- length(parameters)/2
  if (is.null(keep.prob)) { keep.prob <- rep(1,L) }
  # Build cache
  cache <- list()
  for (i in 1:L) {
    # Weights
    temp.W <- parameters[[str_c('W',i)]]
    temp.b <- parameters[[str_c('b',i)]]
    # Linear neurons
    if (i == 1) {
      temp.Z <- (temp.W %*% X) + temp.b
    } else {
      temp.Z <- (temp.W %*% cache[[str_c('A',i-1)]]) + temp.b
    }
    # Activations
    if (i == L) {
      temp.D <- matrix(1,nrow=nrow(temp.Z))
      temp.A <- aterminal(temp.Z)
    } else {
      temp.D <- matrix(runif(length(temp.Z)) <= keep.prob[i],nrow=nrow(temp.Z))
      temp.A <- lrelu(temp.Z,leaky) * temp.D / keep.prob[i]
    }
    # Storage
    cache[[str_c('Z',i)]] <- temp.Z
    cache[[str_c('A',i)]] <- temp.A
    cache[[str_c('D',i)]] <- temp.D
  }
  # Return
  return(cache)
};

########################################
# -------- BACK PROPOGATION ---------- #

backward.propagation.surv <- function(X, delta, cache, PiMat, weights, dterminal, leaky=0, keep.prob=NULL, lam=0) {
  m <- ncol(X)
  ndelta <- sum(unlist(delta))
  L <- length(weights)/2
  gradients <- list()
  if (is.null(keep.prob)) { keep.prob <- rep(1,L) }
  # Reverse
  for (i in seq(L,1)) {
    temp.D <- cache[[str_c('D',i)]]
    if (i == L) { # derivative w.r.t output layer
      temp.dA <- NULL
      temp.dZ <- matrix(unlist(mapply(function(dd,PP) dterminal(dd,PP), delta, PiMat)),nrow=1) #dterminal(delta,PiMat)
    } else { # standard
      temp.dA <- ( t(weights[[str_c('W',i+1,sep='')]]) %*% gradients[[str_c('dZ',i+1)]] ) * temp.D / keep.prob[i]
      temp.dZ <- temp.dA * ifelse((cache[[str_c('A',i)]] > 0),1,leaky) # leaky-relu g'
    }
    if (i == 1) {
      temp.dW <- (1/ndelta) * temp.dZ %*% t(X) + lam*weights[[str_c('W',i)]]
    } else {
      temp.dW <- (1/ndelta) * temp.dZ %*% t(cache[[paste('A',i-1,sep='')]]) + lam*weights[[str_c('W',i)]]
    }
    temp.db <- (1/ndelta) * apply(temp.dZ,1,sum)
    # store
    gradients[[str_c('dA',i)]] <- temp.dA
    gradients[[str_c('dZ',i)]] <- temp.dZ
    gradients[[str_c('dW',i)]] <- temp.dW
    gradients[[str_c('db',i)]] <- temp.db
  }
  return(gradients)
};

#########################################
# -------- UPDATE PARAMETERS ---------- #

update.parameters.surv <- function(parameters, gradients, learning.rate, b1, b2, eps=1e-8) {
  L <- length(parameters$weights) / 2 # number of layers in the neural networks
  # Update rule for each parameter
  if (is.null(parameters$vs)) {
    for (k in seq(1,L)) {
      strW <- paste('W',k,sep='')
      strdW <- paste('dW',k,sep='')
      strb <- paste('b',k,sep='')
      strdb <- paste('db',k,sep='')
      parameters$weights[[strW]] <- parameters$weights[[strW]] - learning.rate * gradients[[strdW]]
      parameters$weights[[strb]] <- parameters$weights[[strb]] - learning.rate * gradients[[strdb]]
    }
  } else { # Using adam update
    for (k in seq(1,L)) {
      strW <- paste('W',k,sep='')
      strdW <- paste('dW',k,sep='')
      strb <- paste('b',k,sep='')
      strdb <- paste('db',k,sep='')
      # Moving average of the gradients (v)
      parameters$vs$v[[strdW]] <- b1*parameters$vs$v[[strdW]] + (1-b1)*gradients[[strdW]]
      parameters$vs$v[[strdb]] <- b1*parameters$vs$v[[strdb]] + (1-b1)*gradients[[strdb]]
      # Moving average of the squared gradients
      parameters$vs$s[[strdW]] <- b2*parameters$vs$s[[strdW]] + (1-b2)*gradients[[strdW]]^2
      parameters$vs$s[[strdb]] <- b2*parameters$vs$s[[strdb]] + (1-b2)*gradients[[strdb]]^2
      # Update the parameters
      parameters$weights[[strW]] <- parameters$weights[[strW]] - learning.rate * (parameters$vs$v[[strdW]]/(sqrt(parameters$vs$s[[strdW]])+eps))
      if (k == L) { # Do not use intercept in the terminal node
        parameters$weights[[strb]] <- parameters$weights[[strb]]
      } else {
        parameters$weights[[strb]] <- parameters$weights[[strb]] - learning.rate * (parameters$vs$v[[strdb]]/(sqrt(parameters$vs$s[[strdb]])+eps))  
      }
      
    }
  }
  # return 
  return(parameters)
};

#################################################
# -------- START FROM CUSTOM WEIGHTS ---------- #

weight.trans <- function(weights,init) {
  # Old weight names
  namWeights <- names(weights)
  # Double check they match
  stopifnot(all(namWeights %in% names(init)))
  # Assign if dimensionals match
  for (nn in namWeights) {
    if (class(init[[nn]]) == 'numeric') {
      if (length(init[[nn]])==length(weights[[nn]])) {
        init[[nn]] <- weights[[nn]]
      }
    }
    if (class(init[[nn]]) == 'matrix') {
      if ((nrow(init[[nn]])==nrow(weights[[nn]])) & (ncol(init[[nn]])==ncol(weights[[nn]]))) {
        init[[nn]] <- weights[[nn]]
      }
    }
  }
  # Return
  return(init)
};

##########################################
# -------- NEURAL NET WRAPPER ---------- #
 
nnet.coxph <- function(X, So, layers, strat = NULL, X.test=NULL, So.test=NULL, strat.test=NULL,
                       learning.rate=0.01, eta=1, num.epochs=1500, leaky=0,
                       weights = NULL, sgd=T, sig=1,
                       lam=0, keep.prob = NULL, adam=F, b1=0.9, b2=0.999, eps=1e-8, ss=1,
                       standardize=T, verbose=T, checkgrads=F, coxbench=F, 
                       aterminal=linear.activation, dterminal=partial.backprop, loss=partial.loss) {
  
  # --- Initialize holders --- #
  if (standardize) {
    if (verbose) print('Scaling')
    vecMu <- apply(X,1,mean,na.rm=T)
    vecSd <- apply(X,1,sd,na.rm=T)
    vecMu <- ifelse(is.na(vecMu),0,vecMu)
    vecSd <- ifelse(is.na(vecSd) | vecSd==0,1,vecSd)
    X <- sweep(sweep(X,1,vecMu,'-'),1,vecSd,'/')
  } else {
    vecMu <- rep(0,nrow(X))
    vecSd <- rep(1,nrow(X))
  }
  # Scale the X-test
  if (!is.null(X.test)) { X.test <- sweep(sweep(X.test,1,vecMu,'-'),1,vecSd,'/') }
  # Check for NAs after standardization
  if (any(is.na(X))) { X[is.na(X)] <- 0 }
  if (!is.null(X.test)) {
    if(!is.null(X.test) & any(is.na(X.test))) { X.test[is.na(X.test)] <- 0 }
  }
  stopifnot(all(rownames(X) == rownames(X.test)))
  
  train.costs <- rep(NA,floor(num.epochs/25))
  train.concs <- test.costs <- test.concs <- train.costs
  
  # --- Set up parameters --- #
  m <- ncol(X)
  n0 <- nrow(X)
  layer.dims <- c(n0,layers) 
  L <- length(layer.dims) - 1 # Take off input layer
  # Parameter weights
  init <- initialize.weights.he(layer.dims,ss,sig)
  if (is.null(weights)) {
    weights <- init
  } else {
    weights <- weight.trans(weights,init)
  }
  # Adam-gradient descent
  if (adam) {
    vs <- initialize.adam(layer.dims)
  } else {
    vs = NULL
  }
  parameters <- list(weights=weights, vs=vs)
  
  # Names of the Weight parameters
  namWeights <- str_c('W',1:L)
  namAL <- str_c('A',L)
  namGrads <- c(str_c('dW',1:L),str_c('db',1:L))

  # --- Strata --- #
  if (is.null(strat)) {
    strat <- rep(1,m)
  } else {
    strat <- as.numeric(factor(strat,levels=unique(strat)))
  }
  K <- length(unique(strat))
  # Make sure patients are ordered by strata
  stopifnot(sum(diff(strat))==(K-1))
  
  # --- Set up all objects into lists --- #
  lst.event <-  vector(mode='list',length=K)
  lst.stime <- vector(mode='list',length=K)
  lst.delta <- vector(mode='list',length=K)
  lst.Pimat <- vector(mode='list',length=K)
  lst.Ymat <- vector(mode='list',length=K)
  lst.tYmat <- vector(mode='list',length=K)
  for (k in 1:K) {
    idx.strat <- which(strat == k)
    lst.stime[[k]] <- So[idx.strat,1]
    lst.delta[[k]] <- So[idx.strat,2]
    lst.event[[k]] <- sum(lst.delta[[k]])
    lst.Ymat[[k]] <- Ymat.fun(lst.stime[[k]])
    lst.tYmat[[k]] <- t(lst.Ymat[[k]])
  }
  stopifnot(all(unlist(lapply(lst.stime,function(ll) all(diff(ll)>=0)))))
  # Total delta
  vec.events <- unlist(lst.event)
  tevents <- sum(vec.events)
  prob.events <- cumsum(vec.events/tevents) ## Relative sampling prob
  oneofK <- 1:K
  # Repeat for the tests
  if (!is.null(X.test)) {
    K.test <- length(unique(strat.test))
    strat.test <- as.numeric(as.factor(strat.test))
    lst.test.Ymat <- vector(mode='list',length=K.test)
    lst.test.delta <- vector(mode='list',length=K.test)
    for (k in seq(K.test)) {
      idx.test.strat <- which(strat.test == k)
      lst.test.Ymat[[k]] <- Ymat.fun(So.test[idx.test.strat,1])
      lst.test.delta[[k]] <- So.test[idx.test.strat,2]
    }
    tevents.test <- sum(unlist(lst.test.delta))  
  }
  
  # --- GET THE LOSS FROM COXPH --- #
  if (coxbench) {
    mdl.coxph <- coxph(So ~ . + strata(strat),data=data.frame(t(X)))
    cost.cox <- 0
    for (k in 1:K) { cost.cox <- cost.cox + loss(predict(mdl.coxph)[which(strat == k)],lst.delta[[k]],lst.Ymat[[k]]) }
    cost.cox <- cost.cox/tevents + (lam/2)*sum(coef(mdl.coxph)^2) # Divide by number of events and add on L2 reg
    conc.cox <- survConcordance(So ~ predict(mdl.coxph) + strata(strat))$conc
    # Repeat for test
    if(!is.null(X.test)) {
      cox.pred.test <- as.vector(t(X.test) %*% coef(mdl.coxph))
      conc.cox.test <- survConcordance(So.test ~  cox.pred.test + strata(strat.test))$conc
      cost.cox.test <- 0
      for (k in 1:K.test) { cost.cox.test <- cost.cox.test + 
        loss(cox.pred.test[which(strat.test == k)],lst.test.delta[[k]],lst.test.Ymat[[k]])
      }
      cost.cox.test <- cost.cox.test/tevents.test
    }
  } else { cost.cox = conc.cox= NA }
  
  print('Beginning the descent')
  # --- Gradient descent --- #
  for (i in 1:num.epochs) {
    # Forward prop
    cache <- forward.propagation.surv(X, parameters$weights, aterminal, leaky, keep.prob)
    # Calculate the cohort-level Pfun
    for (k in 1:K) {
      # Update the Pfun
      lst.Pimat[[k]] <- Pfun(cache[[namAL]][strat == k],lst.Ymat[[k]],lst.tYmat[[k]])
    }
    # Backward propagation
    grads <- backward.propagation.surv(X, lst.delta, cache, lst.Pimat, parameters$weights, 
                                            dterminal, leaky, keep.prob, lam)
    # Update parameters
    parameters <- update.parameters.surv(parameters, grads, learning.rate, b1, b2, eps) # grads
    # Update learning rate
    learning.rate <- eta * learning.rate
    # Get an update every 25th epoch
    if (((i %% 100) == 0) & verbose ) {
      # Get the current L2
      l2 <- sum(unlist(lapply(parameters$weights[namWeights],function(mm) sum(mm^2))))
      # Training Loss
      train.cost <- 0
      for (k in 1:K) { train.cost <- train.cost + loss(cache[[namAL]][strat == k],lst.delta[[k]],lst.Ymat[[k]]) }
      train.cost <- train.cost/tevents + (lam/2)*l2
      train.conc <- survConcordance(So ~ as.vector(cache[[namAL]]) + strata(strat))$conc
      train.costs[i/100] <- train.cost
      train.concs[i/100] <- train.conc
      # Print training
      print(sprintf('--epoch %i-- TRAINING: nnet-loss: %0.3f, cox-loss: %0.3f, nnet-conc: %0.3f, cox-conc %0.3f',
                    i,train.cost,cost.cox,train.conc,conc.cox))
      # If there is a test set get accuracy
      if(!is.null(X.test)) {
        test.cache <- forward.propagation.surv(X.test,parameters$weights,aterminal,leaky,keep.prob)
        test.cache.ahat <- as.vector(test.cache[[namAL]])
        test.cost <- 0
        for (k in 1:K.test) { test.cost <- test.cost + 
          loss(test.cache.ahat[strat.test == k],lst.test.delta[[k]],lst.test.Ymat[[k]])
        }
        test.cost <- test.cost/tevents.test
        test.conc <- survConcordance(So.test ~ test.cache.ahat  + strata(strat.test))$conc
        print(sprintf('--epoch %i-- TESTING: nnet-loss: %0.3f, cox-loss: %0.3f, nnet-conc: %0.3f, cox-conc %0.3f',
                      i,test.cost,cost.cox.test,test.conc,conc.cox.test))
        test.costs[i/25] <- test.cost
        test.concs[i/25] <- test.conc
      }
    }
  }
  # Check gradient
  if (checkgrads) {
    gdiff <- gradcheckfun.surv(X,So,lst.delta,lst.Ymat,lst.tYmat,strat,aterminal,dterminal,loss,
                               parameters$weights,lam,leaky)
    print(sprintf('Gradient check pass: %s, gradient error: %0.20f',gdiff < 1e-6,gdiff)); stopifnot(gdiff < 1e-6)
  }
  # Final training predictions
  ahat <- as.vector(forward.propagation.surv(X,parameters$weights,aterminal)[[namAL]])
  # data.frame of costs
  cost.df <- data.frame(train.costs,train.concs,test.costs,test.concs)
  # return
  ret.list <- list(weights=parameters$weights,aterminal=aterminal,
                   ahat=ahat,#coxweights=coef(mdl.coxph),
                   leaky=leaky,
                   vecMu=vecMu,vecSd=vecSd,
                   cost.df=cost.df)
  return(ret.list)
};


###############################################
# -------- POST-ESIMATION WRAPPERS ---------- #

# Prediction wrapper
predict.snet <- function(Xraw,nnetmdl) {
  m <- ncol(Xraw)
  vecMu <- nnetmdl$vecMu
  vecSd <- nnetmdl$vecSd
  # Scale X
  X <- sweep(sweep(Xraw,1,nnetmdl$vecMu,'-'),1,nnetmdl$vecSd,'/')
  # Forward prop
  cache <- forward.propagation.surv(X,nnetmdl$weights,nnetmdl$aterminal, nnetmdl$leaky)
  zhat <- cache[[str_c('Z',length(cache)/3)]]
  ahat <- nnetmdl$aterminal(zhat)
  return(ahat)
};


##############################################
# -------- HYPER PARAMETER TUNING ---------- #

tune.learnrate <- function(X,So,layers,strat=NULL,X.test,So.test,lam=0) {
  nrate <- 6
  rate.seq <- 10^(-seq(1,nrate))
  conc.store <- rep(NA,nrate)
  for (k in seq_along(rate.seq)) {
    temp.mdl <- snet.model(X,So,layers,strat,verbose = F,weights = NULL,lam=lam,ss=1,
                           learning.rate = rate.seq[k],
                           num.epochs=250,adam=T,standardize = T,sig = 1/4,leaky = 0.01)
    temp.conc <- survConcordance(So.test ~ as.vector(predict.snet(X.test,temp.mdl)))$conc
    conc.store[k] <- temp.conc
  }
  temp.slice <- data.frame(rate=rate.seq[which.max(conc.store)],conc=max(conc.store))
  return(temp.slice)
};