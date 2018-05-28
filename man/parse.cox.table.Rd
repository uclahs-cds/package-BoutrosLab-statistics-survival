\name{parse.cox.table}
\alias{parse.cox.table}
\title{
parse.cox.table
}
\description{
A useful function to format output of coxph if the desired model is too complicated for fit.coxmodel.
}
\usage{
parse.cox.table(cox.data, table.rownames = NULL, zph.p = NULL)
}
\arguments{
  \item{cox.data}{Output of call to coxph}
  \item{table.rownames}{Optional vector of rownames for output table}
  \item{zph.p}{An optional vector of p-values from cox.zph call. If NULL (default), cox.zph will be called in the function.}
}
\value{
A table with HR, 95percent CI, p-value and cox.zph p-value per variable in the cox model. Last row contains wald and log-rank p value.
}
