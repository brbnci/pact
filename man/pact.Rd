% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pact.R
\docType{package}
\name{pact}
\alias{pact}
\alias{pact-package}
\title{Predictive Analysis of Clinical Trials}
\description{
The \code{pact} package implements a prediction-based approach to the 
analysis of data from randomized clinical trials (RCT). Based on clincial response and 
covariate data from a RCT comparing a new experimental treatment E versus a control C, 
the purpose behind the functions in \code{pact} is to develop and internally validate
a model that can identify subjects likely to benefit from E rather than C. Currently, 
'survival' and 'binary' response types are permitted.
}
\details{
Package: pact
\cr Type: Package
\cr Version: 0.5.0
\cr Date: 2016-04-14
\cr Author: Dr. Jyothi Subramanian and Dr. Richard Simon
\cr Maintainer: Jyothi Subramanian <\email{subramanianj01@gmail.com}>
\cr License: GPL-3
\cr\cr \code{pact.fit} fits a predictive model to data from RCT. Currently, 'survival' and 
'binary' response types are supported. Analysis of high dimensional covariate data is supported.
If known and available, a limited number of prognostic covariates can also be specified and 
fixed to remain in the predictive model. An object of class 'pact' is returned. 
\code{print}, \code{summary} and \code{predict} methods are available for objects of 
class 'pact'.
Additionally, the function \code{pact.cv} takes as an input the object returned by \code{pact.fit} 
and computes predictive scores for each subject through k-fold cross-validation. 
Evaluations of the cross-validated predictions are performed by the function \code{eval.pact.cv}. 
\cr\cr  Finally, the function \code{overall.analysis}
also takes an object of class 'pact' as input and computes some summary statistics 
for the comparison of treatments E and C.
}
\examples{
### Survival response
set.seed(10)    
data(prostateCancer)     
Y <- prostateCancer[,3:4]
Xf <- prostateCancer[,7:8]  ## Prognostic covariates fixed to always be in the model
Xv <- prostateCancer[,c(5:6,9)]
Treatment <- as.factor(prostateCancer[,2])
p <- pact.fit(Y=Y,Xf=Xf,Xv=Xv,Treatment=Treatment,family="cox",varSelect="univar")
print(p)
overall.analysis(p)
cv <- pact.cv(p, nfold=5)
eval.pact.cv(cv, method="continuous", plot.score=TRUE, perm.test=FALSE, nperm=100)

### Binary response
set.seed(10)
data(EORTC10994)
Y <- as.factor(EORTC10994[,4])
## No prognostic covariates (Xf) specified
Xv <- EORTC10994[,c(2,5:7)]
Treatment <- as.factor(EORTC10994[,3])
p <- pact.fit(Y=Y,Xv=Xv,Treatment=Treatment,family="binomial",varSelect="none")
print(p)
overall.analysis(p)
cv <- pact.cv(p, nfold=5)
eval.pact.cv(cv, method="discrete", g=log(1), perm.test=FALSE, nperm=100)

### High dimensional data, survival response
\dontrun{
set.seed(10)    
data(GSE10846)     
Y <- GSE10846[,1:2]
Xv <- GSE10846[,-c(1:3)]
Treatment <- as.factor(GSE10846[,3])
p <- pact.fit(Y=Y,Xv=Xv,Treatment=Treatment,family="cox",varSelect="lasso",penalty.scaling=2)
print(p)
overall.analysis(p)
cv <- pact.cv(p, nfold=5)
eval.pact.cv(cv, method="continuous", plot.score=TRUE, perm.test=FALSE)
}

}

