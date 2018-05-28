\name{multi.point.HR.table}
\alias{multi.point.HR.table}
\title{
multi.point.HR.table
}
\description{
A function to estimate HR at different time points (time windows) when Cox PH assumption fails.
}
\usage{
multi.point.HR.table(all.groups = NULL, all.survtime = NULL, all.survstat = NULL,
     truncation.thresholds = c(5, 10), covariates = NULL)
}
\arguments{
  \item{all.groups}{Grouping of patients (A vector containing survival groups}
  \item{all.survtime}{A vector containing survival times}
  \item{all.survstat}{A vector containing survival status}
  \item{truncation.thresholds}{A vector of end points of time windows for which HR is estimated}
  \item{covariates}{A data frame with covariates to be included in the cox model}
}
\value{
A matrix of HR, 95L CI, 95U CI, Wald P, n (samples) and cox.zph P for the full time frame, and for each time window established using \code{truncation.thresholds}
}
\examples{
# set up survival data
set.seed(100);

survtime.ph <- c(
	sample(1:100, 100, replace = TRUE),
	sample(1:60, 50, replace = TRUE),
	sample(40:100, 50, replace = TRUE)
	);

survstat.ph <- c(
	sample( c(0,1), 100, replace = TRUE),
	sample( c(0,0,0,0,1), 50, replace = TRUE),
	sample( c(0,1,1,1,1), 50, replace = TRUE)
	);
groups.ph <- factor(c(rep('A', 100), rep('B', 100)), levels = c('A','B'));

# try example without covariates, specifying cut points
multi.point.HR.table(
	all.groups = groups.ph,
	all.survtime = survtime.ph,
	all.survstat = survstat.ph,
	truncation.thresholds = c(60)
	);

# now try example with covariates
gender.ph <- factor(sample(c("Male","Female"), 200, replace = TRUE), levels = c("Male","Female"));
score.ph <- rgamma(n=200, shape = 20, rate=1); # this is just an arbitrary "continuous" variable
multi.point.HR.table(
	all.groups = groups.ph,
	all.survtime = survtime.ph/12, 
	all.survstat = survstat.ph,
	covariates = data.frame(gender = gender.ph, score = score.ph)
	);

}
\keyword{survival}
