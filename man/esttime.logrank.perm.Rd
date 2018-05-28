\name{esttime.logrank.perm}
\alias{esttime.logrank.perm}
\title{Estimate run time for Log-Rank Permutation Test}
\description{Calculates the estimated run time for performing a log-rank permutation test with N patients with k variants}
\usage{
esttime.logrank.perm(N,k,no_cores=1,strat=FALSE,tform='secs',verbose=TRUE);
}
\arguments{
        \item{N}{The number of patients}
        \item{k}{The number of patients in one of the groups}
        \item{no_cores}{Number of cores that will be used to parallelize (efficiency peaks at 6)}
        \item{strat}{Boolean whether there will be a stratification}
        \item{tform}{The returned time value (default minutes). Options include: secs, mins, hours}
        \item{verbose}{Boolean whether to print internal values}
}
\value{Returns the expected run-time. Note, ignores time required to initialize clusters and create holder matrix. Warning, estimate is based on a linear extrapolation, so small errors will magnify for very large number of calculations.}
\author{Erik Drysdale}
\examples{
# library(survival)
# N <- 50
# k <- 4
# tt <- rexp(N,1)
# delta <- rep(1,N)
# So <- as.matrix(Surv(tt,delta))
# x <- sample(c(rep(1,k),rep(0,N-k)))
# ss <- sample(rep(c('a','b'),ceiling(N/2)))[1:N]
# # 50 choose 4, no strata, one core
# esttime.logrank.perm(N=N,k=k,no_cores=1,strat=FALSE,tform='secs',verbose=T)
# print(system.time({logrank.perm(So,x,strat=NULL,verbose=F)}))
# # 50 choose 4, no strata, two core core
# esttime.logrank.perm(N=N,k=k,no_cores=2,strat=FALSE,tform='secs',verbose=T)
# print(system.time({logrank.perm(So,x,strat=NULL,no_cores = 2,verbose=F)}))
# # 50 choose 4, with strata, one core
# esttime.logrank.perm(N=N,k=k,no_cores=1,strat=TRUE,tform='secs',verbose=T)
# system.time({logrank.perm(So,x,strat=ss,verbose=F)})
# # 50 choose 4, with strata, two core
# esttime.logrank.perm(N=N,k=k,no_cores=2,strat=TRUE,tform='secs',verbose=T)
# system.time({logrank.perm(So,x,strat=ss,no_cores=2,verbose=F)})
}
