\name{calculate.number.at.risk}
\alias{calculate.number.at.risk}
\title{Calculate the number of patients at risk}
\description{From a survival object and a patient-grouping vector, calculate how many patients are at risk at a vector of time-points.}
\usage{
calculate.number.at.risk(survobj, cut.points);
}
\arguments{
  \item{survobj}{An object of class Surv (from the survival package)}
  \item{cut.points}{A vector of time-points at which the number of patients at risk are calculated}
}
\value{A vector of length length(cut.points) giving the number of patients at risk at each time-point}
\author{Paul C. Boutros}
\examples{
# create some synthetic survival data
survtime <- sample(seq(0.1,10,0.1), 100, replace = TRUE);
survstat <- sample(c(0,1), 100, replace = TRUE);
survobj <- Surv(survtime, survstat);

# number of patients at risk at two-year intervals
calculate.number.at.risk(
	survobj = survobj,
	cut.points = seq(0,10,2)
	);
}
\keyword{survival}
