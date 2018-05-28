\name{get.survival.rate}
\alias{get.survival.rate}
\title{Find the event-free rate}
\description{Find the event-free rates of groups of patients using the Kaplan Meier method.}
\usage{
get.survival.rate(surv.obj, groups, cut.points);
}
\arguments{
  \item{surv.obj}{A Surv object}
  \item{groups}{A factor indicating patient grouping}
  \item{cut.points}{A vector of time points to be checked}
}
\value{A matrix of levels(groups) x cut.points containing the event-free rates for each combination.}
\author{Emilie Lalonde}
