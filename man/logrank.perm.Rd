\name{logrank.perm}
\alias{logrank.perm}
\title{Log-Rank Permutation Test}
\description{Calculates log-rank statistic for all (N k) possible outcomes and associated p-value under null of no survival differences}
\usage{
logrank.perm(So, x, rho=0, strat=NULL, no_cores=1, verbose=TRUE);
}
\arguments{
        \item{So}{Matrix-formatted Surv() object}
        \item{x}{Binary vector for two groups}
        \item{rho}{A scalar parameter that controls the type of test (see survdiff())}
        \item{strat}{Argument to be passed onto strata()}
        \item{no_cores}{Number of cores to use}
        \item{verbose}{Prints updates during function execution}
}
\value{Returns a data.frame with the number of tests, the number which exceed the baseline, the logrank p-value (logrank) and the exact permutation-based p-value (exact)}
\author{Erik Drysdale}
\examples{
# # WARNING! FOR LARGE NUMBERS OF CALCULATIONS, REQUEST SUFFICIENT RAM OR PROGRAM MAY CRASH! #
# # Example 
# library(survival)
# library(BoutrosLab.statistics.survival)
# idx <- c(1:8,16:23)
# Sob <- Surv(time=aml$time,event=aml$status)[idx,]
# So <- as.matrix(Sob)
# # Encode one of the groups
# x <- ifelse(aml$x=="Maintained",1,0)[idx]
# # Run the classical log-rank test
# lrtest <- survdiff(Sob ~ x, rho=0)
# # Permutation version
# permtest <- logrank.perm(So,x,rho=0,no_cores=2,verbose=F)
# # Compare p-values
# print(sprintf("p-value from log-rank %0.3f and p-value from exact permutation test %0.3f",permtest["logrank"],permtest["exact"]))
# # Compare finite-sample bias to traditional log-rank test
# nsim <- 500; N <- 25; k <- 2
# storage <- data.frame(matrix(NA,nrow=nsim,ncol=2))
# colnames(storage) <- c('chisq','exact')
# # Begin loop
# for (i in 1:nsim) {
#   set.seed(i)
#   idx <- sample(1:N,k)
#   x <- rep(0,N)
#   x[idx] <- 1
#   strat <- sample(rep(c("a","b"),ceiling(N/2)))[1:N]
#   Tobs <- rexp(N,1)
#   is.event <- rep(1,N)
#   Sob <- Surv(time=Tobs,event=is.event)
#   So <- as.matrix(Sob)
#   lrtest <- survdiff(Surv(time=Tobs,event=is.event)~x+strata(strat))
#   pval.lrtest <- pchisq(lrtest$chisq,1,lower.tail=F)
#   permtest <- logrank.perm(So=So,x=x,strat = strat,verbose=F)
#   pval.permtest <- permtest$exact
#   storage[i,] <- c(pval.lrtest,pval.permtest)
# }
# coverage <- apply(storage,2,function(cc) mean(cc <= 0.05,na.rm = F))
# print(coverage)
}
