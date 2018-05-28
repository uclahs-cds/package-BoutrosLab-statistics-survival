\name{ph.fails}
\alias{ph.fails}
\title{Utility function to test Proportional Hazards (PH) assumption}
\description{A utility function that tests the Proportional Hazards (PH) assumption.}
\usage{
ph.fails(cox.model, cox.zph.threshold = 0.1, pvalues = FALSE);
}
\arguments{
  \item{cox.model}{The coxph object for which the PH assumption will be checked.}
  \item{cox.zph.threshold}{The threshold of pvalue for the cox.zph test below which warnings and/or residual plots will be generated.  Defaults to 0.1}
  \item{pvalues}{If TRUE, returns a 1-row matrix with the pvalues from the cox.zph test.  If FALSE, returns a TRUE/FALSE indicator of whether or not any of the pvalues are smaller than cox.zph.threshold.  Defaults to FALSE.}
}
\value{If pvalues == TRUE, returns a matrix with pvalues from the cox.zph test.  If pvalues == FALSE, returns a boolean indicator of whether or not the PH assumption fails.}
\author{Nathalie Moon}
\seealso{\code{\link{fit.coxmodel}}}
\keyword{hplot}
