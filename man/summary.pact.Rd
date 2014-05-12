\name{summary.pact}
\alias{summary.pact}
\title{Summarize a predictive model fit}
\usage{
\method{summary}{pact}(object, ...)
}
\arguments{
  \item{object}{The object returned from \code{'pact.fit'}}

  \item{...}{Additional arguments for 'summary'}
}
\value{
The summary statistics for the regression model fitted by
\code{pact.fit} is returned
}
\description{
summary method for objects of class 'pact'
}
\details{
Returns the summary statistics for the regression model of
the 'pact' object in a tabular format
}
\examples{
data(prostateCancer)
Y <- prostateCancer[,3:4]
X <- prostateCancer[,5:9]
Treatment <- prostateCancer[,2]
p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="cox")
summary(p)
}
\author{
Jyothi Subramanian and Richard Simon \cr Maintainer: Jyothi
Subramanian <\email{subramanianj01@gmail.com}>
}

