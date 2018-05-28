\name{calculate.median.followUp.time}
\alias{calculate.median.followUp.time}
\title{Calculate the median follow up time for a patient cohort}
\description{Calculate the median follow up time for a patient cohort, excluding patients that have already died.}
\usage{
calculate.median.followUp.time(fu.time, vital.status, method);
}
\arguments{
  \item{fu.time}{A vector with follow-up time}
  \item{vital.status}{A vector with vital status for each patient. Must be numeric (1 = death).}
  \item{method}{Either 1 or 2, where 1 uses only alive patients to calculate median follow up time, and 2 uses the Kaplan-Meier estimate of potential follow-up. See http://cemsiis.meduniwien.ac.at/fileadmin/msi_akim/CeMSIIS/KB/volltexte/Schemper_Smith_1996_CCT.pdf }
}
\value{Median follow-up time}
\author{Emilie Lalonde}
\examples{
# create some synthetic survival data
survtime <- sample(seq(0.1,10,0.1), 100, replace = TRUE);
vitalstat <- sample(c(0,1), 100, replace = TRUE);

# number median follow-up time
calculate.median.followUp.time(
	fu.time = survtime,
	vital.status = vitalstat,
	method = 1
	);
calculate.median.followUp.time(
	fu.time = survtime,
	vital.status = vitalstat,
	method = 2
	);
}
\keyword{survival}
