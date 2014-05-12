\name{print.pact}
\alias{print.pact}
\title{Print an object of class 'pact'}
\usage{
\method{print}{pact}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{The object returned from \code{'pact.fit'}}

  \item{digits}{significant digits in the print}

  \item{...}{Additional print arguments}
}
\value{
The classification function is printed
}
\description{
print method for objects of class 'pact'
}
\details{
The call that produced the object is printed, followed by
the classification function from \code{pact.fit} for
calculating the predictive scores for new subjects
}
\examples{
data(prostateCancer)
Y <- prostateCancer[,3:4]
X <- prostateCancer[,5:9]
Treatment <- prostateCancer[,2]
p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="cox")
print(p)
}
\author{
Jyothi Subramanian and Richard Simon \cr Maintainer: Jyothi
Subramanian <\email{subramanianj01@gmail.com}>
}

