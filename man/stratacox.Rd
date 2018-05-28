\name{stratacox}
\alias{stratacox}
\title{Elastic-Net Stratified Cox Model}
\description{Implements a Stratified Cox model with elastic-net regularization and a transfer-weighting scheme}
\usage{
stratacox(X,So,strat,lam,alpha,tau=0,target=1,w=NULL,
                          standardize=T,ord=T,max.iter=250,tol=1e-8,findlam=F,verbose=F)
}
\arguments{
        \item{X}{A model.matrix type object}
        \item{So}{A matrix from a as.matrix(Surv()) object}
        \item{strat}{A vector of integers/characters specifying each row to a respective strata}
        \item{lam/alpha}{Elastic net regularization parameters (see glmnet)}
        \item{tau}{Relative weighting to be put on target dataset, where 0 is fair weighting (dataset size) and 1 uses only the target's likelihood}
        \item{target}{Integer specifying which of the strata levels is to be used as the target}
        \item{w}{Initial weighting}
        \item{standardize}{Boolean whether to standardize the features}
        \item{ord}{Boolean whether the dataset has been ordered by strata/time already}
        \item{max.iter}{Maximum number of gradient steps to use}
        \item{tol}{L2 norm difference used to end gradient dsecent}
        \item{findlam}{Boolean whether to find the largest lambda to set all features equal to zero}
        \item{verbose}{Boolean whether to print output during training}
}
\value{Returns a list with: w (coefficient weights), ll (log-likelihood iterations), step (Barzilai-Borwein step size iterations), and niter (number of steps until convergence).}
\author{Erik Drysdale}

\examples{
# # libraries
# library(BoutrosLab.statistics.survival)
# library(survival)
# library(glmnet)
# library(forcats)
# library(data.table)
# # Can stratified analysis help?
# dir.base <- '/.mounts/labs/cpcgene/private/projects/GermlineSomaticAssociations/genome-wide/input/all/all'
# setwd(dir.base)
# # Load IGRT/RadP datasets
# dat <- rbind(read.table('2017-11-16_IGRT-patient-matrix.txt',header=T),
#              read.table('2017-11-16_RadP-patient-matrix.txt',header=T))
# dat$gleason <- fct_recode(dat$gleason,'3+4'='3+5','4+3'='4+4')
# dat$gleason <- as.factor(ifelse(is.na(dat$gleason),'3+4',as.character(dat$gleason)))
# # Survival and design
# Xdat <- data.table(model.matrix(~primary_treatment+time_to_bcr+bcr+gleason+t_category+psa+age+ets_consensus+
#                                   kataegis_score+kataegis_binary+pga,data=dat)[,-1])
# colnames(Xdat)[1:3] <- c('strat','time','event')
# colnames(Xdat) <- gsub('\\+','\\_',colnames(Xdat))
# Xdat[,strat := ifelse(strat==0,'IGRT','RadP')]
# Xdat <- Xdat[order(strat,time)]
# n.IGRT <- sum(Xdat$strat == 'IGRT')
# n.RadP <- sum(Xdat$strat == 'RadP')
# idx.IGRT <- which(Xdat$strat == 'IGRT')
# idx.RadP <- which(Xdat$strat == 'RadP')
# # Get baseline distribution
# nsim <- 1000
# cn <- c('coxph','strata','wstrata')
# store <- data.frame(matrix(NA,nrow=nsim,ncol=length(cn)))
# colnames(store) <- cn
# for (k in 1:nsim) { 
#   if (k %% 25 == 0) print(k)
#   set.seed(k)
#   # Split data
#   idx.IGRT.train <- sort(sample(idx.IGRT,floor(n.IGRT*0.8)))
#   idx.IGRT.test <- setdiff(idx.IGRT,idx.IGRT.train)
#   dat.train <- Xdat[c(idx.IGRT.train,idx.RadP)]
#   dat.test <- Xdat[c(idx.IGRT.test)]
#   # Formulas
#   form.coxph <- as.formula(paste0('Surv(time,event)~',paste(colnames(dat.test)[-(1:3)],collapse='+')))
#   form.strata <- as.formula(paste0('Surv(time,event)~strata(strat)+',paste(colnames(dat.test)[-(1:3)],collapse='+')))
#   # Fit models
#   mdl.coxph <- coxph(form.coxph,data=dat.train,ties = 'breslow',subset=(strat=='IGRT'))
#   mdl.naivestrata <- coxph(form.strata,data=dat.train,ties = 'breslow')
#   mdl.wstrata <- stratacox(X=as.matrix(dat.train[,-(1:3)]),
#                                So=as.matrix(with(dat.train[,2:3],Surv(time,event))),
#                                strat=dat.train$strat,
#                                lam=0.01,alpha=1,tau=0.1,target=1,ord = F,
#                                standardize = T,tol = 1e-5,max.iter = 200)
#   # Make predictions
#   eta.coxph <- predict(mdl.coxph,newdata=dat.test)
#   eta.naivestrata <- predict(mdl.naivestrata,newdata=dat.test)
#   eta.wstrata <- as.vector(as.matrix(dat.test[,-(1:3)]) %*% mdl.wstrata$w)
#   # Predictions
#   conc.coxph <- survConcordance(Surv(time,event)~eta.coxph,data=dat.test)$conc
#   conc.naivestrata <- survConcordance(Surv(time,event)~eta.naivestrata,data=dat.test)$conc
#   conc.wstrata <- survConcordance(Surv(time,event)~eta.wstrata,data=dat.test)$conc
#   # Store
#   store[k,] <- c(conc.coxph,conc.naivestrata,conc.wstrata)
# }
# # Gain by using strata's
# t.test(store[,2],store[,1])
# # Gain from weight+regularization
# t.test(store[,3],store[,2])
}
