\name{power.analysis}
\alias{power.analysis}
\title{Power analysis for an expected hazard ratio}
\description{This function has been moved to BoutrosLab.statistics.power and renamed power.cox.univariate. It performs power analysis for expected hazard ratio and number of events}
\usage{
power.analysis(event.number, HR, alpha = 0.05);
}
\arguments{
  \item{event.number}{Number indicating the number of events in a cohort}
  \item{HR}{Hazard ratio}
  \item{alpha}{Statistical level, defaults to 0.05}  
}
\value{power calculated for the given event.number and HR}
\author{Maud HW. Starmans}
\examples{
HRseries <- c(100:250)/100;
Powerseries <- power.analysis(100,HRseries);
}
\keyword{survival}
