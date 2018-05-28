\name{fit.coxmodel}
\alias{fit.coxmodel}
\title{Fit a Cox proportional hazards model}
\description{Fit a Cox model (possibly with some linear adjustments) and return key statistics about the fit.}
\usage{
fit.coxmodel(groups, survobj, other.data = NULL, stratification.factor, 
    stratification.value, rounding = 3, return.cox.model);
}
\arguments{
  \item{groups}{Grouping of patients (passed directly to coxph, so factors & continuous variables are okay)}
  \item{survobj}{An object of class Surv (from the survival package) -- patient ordering needs to be identical as for groups}
  \item{other.data}{A data-frame (or matrix?) of variables to be controlled in the Cox model. If null, no adjustment is done.  No interactions are fit.}
  \item{stratification.factor}{A vector (or factor) containing a categorical variable on which the cox model will be stratified.  If there are too many levels, it is possible that the model will not converge.  Use return.cox.model = TRUE to verify if model looks reasonable.}
  \item{stratification.value}{A pre-determined value used to split the stratification factor. Must be within the range of the variable. (i.e. strata(stratification.factor > stratification.value)). }
  \item{rounding}{Specifies the number of digits of precision to be returned.  Defaults to 3.}
  \item{return.cox.model}{If TRUE, returns coxph object.  If FALSE, returns summary statistics only.  Defaults to FALSE}
}
\value{A vector containing: HR, lower 95\% CI of HR, upper 95\% CI of HR, P-value (for groups), number of samples (total with group assignments, although some may not be included in fit for other reasons so this is an upper-limit).}
\author{Paul C. Boutros}
\examples{
set.seed(151)
survtime <- sample(seq(0.1,10,0.1), 100, replace = TRUE);
survstat <- sample(c(0,1), 100, replace = TRUE);
survobj <- Surv(survtime, survstat);
groups <- sample(c('A','B'), 100, replace = TRUE);
gender <- factor(sample(c('Male','Female'), 100, replace = TRUE), levels=c("Male","Female"));
age <- factor(sample(c('Young','MiddleAged','Old'), 100, replace = TRUE), 
            levels = c('Young', 'MiddleAged', 'Old'));
numbered.age <- sample(seq(40,90), 100, replace = TRUE);

fit.coxmodel(
	groups = as.factor(groups),
	survobj = survobj
	);

fit.coxmodel(
	groups = as.factor(groups),
	survobj = survobj,
	other.data = data.frame(age, gender),
	return.cox.model = TRUE	
	)

fit.coxmodel(
	groups = as.factor(groups),
	survobj = survobj,
	other.data = data.frame(age),
	stratification.factor = gender,
	return.cox.model = TRUE
	)

fit.coxmodel(
	groups = as.factor(groups),
	survobj = survobj,
	other.data = data.frame(gender),
	stratification.factor = numbered.age,
	stratification.value = 60,
	return.cox.model = TRUE
	)
}
\keyword{survival}
