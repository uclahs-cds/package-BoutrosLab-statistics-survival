\name{dichotomize.dataset}
\alias{dichotomize.dataset}
\title{Dichotomize a single dataset}
\description{Split a dataset into two groups by median-dichotomization}
\usage{
dichotomize.dataset(x);
}
\arguments{
  \item{x}{A vector of values to be dichotomized}
}
\value{A vector of the data dichotomized onto a 0/1 (low/high) scale.}
\author{Paul C. Boutros}
\examples{
tmp <- data.frame(y = rnorm(100));
tmp$x <- dichotomize.dataset(tmp$y);
boxplot(y ~ x, tmp);
}
\keyword{survival}
