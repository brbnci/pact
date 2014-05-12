\name{pact.fit}
\alias{pact.fit}
\title{Fits a predictive model to the full dataset}
\usage{
pact.fit(Y, X, Treatment, family = c("binomial", "cox"))
}
\arguments{
  \item{Y}{Response variable. For \code{family='binomial'},
  Y should be a factor with two levels. For
  \code{family='cox'}, Y should be a two-column matrix with
  columns named 'time' and 'status'. The latter is a binary
  variable, with '1' indicating death, and '0' indicating
  right censored.}

  \item{X}{A dataframe of the predictor variables to be
  used for model development. All variables in the
  dataframe are used in the model.}

  \item{Treatment}{The treatment assignement indicator. A
  factor with two levels. '0' indicating control (C) and
  '1' indicating experimental (E) treatment.}

  \item{family}{Type of the response variable. See above.
  Possible values are 'binomial' or 'cox'.}
}
\value{
An object of class 'pact' which is a list with the
following components

\item{reg}{The fitted regression model}

\item{family}{The response variable type used}

\item{Y}{The response variable used}

\item{X}{The dataframe of predictor variables used}

\item{Treatment}{The treatment assignment indicator used}

\item{call}{The call that produced the return object}
}
\description{
\code{pact.fit} Fits a predictive model to the full
dataset. Currently supports Cox and logistic regression
models for 'survival' and 'binary' response types
respectively.
}
\details{
A Cox proportional hazards (PH) regression or a logistic
regression model is developed for data with survival and
binary response respectively. Data from subjects in both
'experimental' (E) and 'control' (C) groups from a RCT is
used for model development. Main effect of treatment, main
effect of covariates and all treatment by covariate
interaction terms are considered in the model development.
}
\examples{
data(prostateCancer)
Y <- prostateCancer[,3:4]
X <- prostateCancer[,5:9]
Treatment <- prostateCancer[,2]
pact.fit(Y=Y, X=X, Treatment=Treatment, family="cox")
}
\author{
Jyothi Subramanian and Richard Simon \cr Maintainer: Jyothi
Subramanian <\email{subramanianj01@gmail.com}>
}
\keyword{pact,}
\keyword{pact.fit}

