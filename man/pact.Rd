\docType{package}
\name{pact}
\alias{pact}
\alias{pact-package}
\title{Predictive Analysis of Clinical Trials}
\description{
The \code{pact} package implements a prediction-based
approach to the analysis of data from randomized clinical
trials (RCT). Based on clincial response and covariate data
from a RCT comparing a new experimental treatment E versus
a control C, the purpose behind the functions in
\code{pact} is to develop and internally validate a model
that can identify subjects likely to benefit from E rather
than C. Currently, 'survival' and 'binary' response types
are permitted.
}
\details{
Package: pact \cr Type: Package \cr Version: 0.1 \cr Date:
2014-04-30 \cr Author: Dr. Jyothi Subramanian and Dr.
Richard Simon \cr Maintainer: Jyothi Subramanian
<\email{subramanianj01@gmail.com}> \cr License: GPL-2
\cr\cr \code{pact.fit} fits a predictive model to a data
from a RCT. Currently, 'survival' and 'binary' response
types are supported. An object of class 'pact' is returned.
\code{print}, \code{summary} and \code{predict} methods are
available for objects of class 'pact'. The function
\code{pact.cv} takes as an input the object returned by
\code{pact.fit} and performs an evaluation of the model
using k-fold cross-validation. Evaluations of the
cross-validated predictions are performed by the function
\code{eval.pact.cv}. \cr\cr Finally, the function
\code{overall.analysis} also takes an object of class
'pact' as input and computes some summary statistics f the
comparison of treatments E and C.
}
\examples{
### Survival response
set.seed(1)
data(prostateCancer)
Y <- prostateCancer[,3:4]
X <- prostateCancer[,5:9]
Treatment <- prostateCancer[,2]
p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="cox")
print(p)
overall.analysis(p)
cv <- pact.cv(p, nfold=5)
eval.pact.cv(cv, method="continuous", plot.score=TRUE, perm.test=FALSE, nperm=100)

### Binary response
set.seed(2)
data(EORTC10994)
Y <- EORTC10994[,4]
X <- EORTC10994[,c(2,5:7)]
Treatment <- EORTC10994[,3]
p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="binomial")
print(p)
overall.analysis(p)
cv <- pact.cv(p, nfold=5)
eval.pact.cv(cv, method="discrete", g=log(1), perm.test=FALSE, nperm=100)
}

