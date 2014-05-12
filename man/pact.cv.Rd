\name{pact.cv}
\alias{pact.cv}
\title{Cross-validation for pact}
\usage{
pact.cv(p, nfold)
}
\arguments{
  \item{p}{An object of class 'pact'}

  \item{nfold}{The number of folds (k) for the k-fold
  cross-validation. k equal to the sample size would mean a
  leave-one-out cross-validation}
}
\value{
A list with the following components

\item{PredScore}{The cross-validated scores for each
subject (a vector)}

\item{Y}{The response variable used}

\item{X}{The dataframe of predictor variables used}

\item{Treatment}{The treatment assignment indicator used}

\item{family}{The response variable type}

\item{call}{The call that produced this output}
}
\description{
k-fold cross-validation for evaluation of the model
developed in \code{pact.fit}
}
\details{
Evaluation of the model developed in \code{pact.fit} using
k-fold cross-validation. In each fold of the
cross-validation, a model is developed from the
observations in the training set. The estimated
coefficients of the regression model developed using
training set are used to make predictions for the left out
observations (test set). This is repeated for all the
folds. Scores are thus obtained for all the subjects in the
dataset. The function \code{\link{eval.pact.cv}} provides
the evaluation method options for the cross-validated
scores.
}
\examples{
data(prostateCancer)
Y <- prostateCancer[,3:4]
X <- prostateCancer[,5:9]
Treatment <- prostateCancer[,2]
p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="cox")
cv <- pact.cv(p, nfold=5)
}
\author{
Jyothi Subramanian and Richard Simon \cr Maintainer: Jyothi
Subramanian <\email{subramanianj01@gmail.com}>
}
\keyword{pact,}
\keyword{pact.cv}

