\name{logrank.analysis}
\alias{logrank.analysis}
\title{Perform a logrank test}
\description{Perform a logrank test and return summary statistics}
\usage{
logrank.analysis(survival.object, groups);
}
\arguments{
  \item{survival.object}{Survival object}
  \item{groups}{Factor variable indicating which entries in 'values' belong to which group for the logrank comparison}
  }
\value{Data frame with relevant statistics}
\author{Nathalie Moon}
\examples{

set.seed(100);

survtime <- c(
	sample(1:75, 100, replace = TRUE),
	sample(1:99, 100, replace = TRUE)
	);

survstat <- c(
	sample( c(0,0,0,1,1), 100, replace = TRUE),
	sample( c(0,0,0,1), 100, replace = TRUE)
	);

survobj <- Surv(survtime, survstat);
two.groups <- factor(c( rep('A', 100), rep('B', 100)), levels = c('A','B'));
three.groups <- factor(c( rep('A', 100), rep('B', 50), rep('C', 50) ), levels = c('A','B','C'));

logrank.analysis(survobj, two.groups);
logrank.analysis(survobj, three.groups);

}
\keyword{htest}
