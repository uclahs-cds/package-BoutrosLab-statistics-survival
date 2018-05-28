\name{CUREcoxph}
\alias{CUREcoxph}
\title{Proportional Hazards Mixture Cure Model (MCM)}
\description{Estimates a PH-MCM with EM algorithm. Returns coefficients from survival (beta) and cure (gamma) models. For reference to notation see http://post.queensu.ca/~pengp/papers/Peng2014_4.pdf.}
\usage{
CUREcoxph(survform,cureform,So,X,bootstrap=NULL,verbose=TRUE,subset=NULL,tol=1e-3,max.iter=50);
}
\arguments{
        \item{survform}{A character or formula that is passed to survival model}
        \item{cureform}{A character or formula that is passed to cure model}
        \item{So}{A Surv() object}
        \item{data}{A data.frame with columns that contain the survform + cureform variables}
        \item{bootstrap}{Number of bootstrap iterations to conduct; if NULL, does not execute}
        \item{verbose}{Prints updates during function execution}
        \item{subset}{Subsets dat to this index}
        \item{tol}{L2 parameter tolerance}
        \item{max.iter}{Maximum number of iterations during EM algorithm}
}
\value{Returns a list with parameters from the survival (beta) and cure (gamma) models, along with 95 percent CI if the bootstrap argument is used.}
\author{Erik Drysdale}

\examples{
# # Simple example 
# library(survival)
# library(BoutrosLab.statistics.survival)
# df <- survival::cancer
# So <- with(df,Surv(time,status==2))
# # Normal Cox model
# mdl.coxph <- coxph(So ~ age + factor(sex),data=df)
# # Cure model
# survform <- '~age+factor(sex)'
# cureform <- '~age+factor(sex)'
# mdl.cure <- CUREcoxph(survform,cureform,So,df,max.iter = 250)
# coef.compare <- data.frame(coxph=c(NA,coef(mdl.coxph)),CUREcox=c(NA,mdl.cure$beta),CUREprob=mdl.cure$gamma)
# rownames(coef.compare) <- c('intercept','age','sex')
# print(round(coef.compare,3))
# 
# # Run Bootstrap
# mdl.bs <- CUREcoxph(survform,cureform,So,df,bootstrap=99,verbose=F,max.iter=250)
# # Show 95 percent CI for PH coefficients (age is insignificant)
# print(mdl.bs$beta)
}
