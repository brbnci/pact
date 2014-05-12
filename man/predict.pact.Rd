\name{predict.pact}
\alias{predict.pact}
\title{Predictions from a predictive model fit}
\usage{
\method{predict}{pact}(object, newx, newTreatment, ...)
}
\arguments{
  \item{object}{The object returned from \code{pact.fit}}

  \item{newx}{Dataframe of predictor values for the new
  subjects for whom predictions are to be made}

  \item{newTreatment}{The treatment assignment indicator
  for the new subjects}

  \item{...}{Other arguments to 'predict'}
}
\value{
A numeric vector containing the predicted scores for the
new subjects from the fitted model is returned
}
\description{
Predicts the scores for new subjects from a previously
developed object of class 'pact'
}
\details{
Returns the scores for new subjects from an object of class
'pact', given their coveriate values and treatment
assignment.
}
\examples{
### Survival response
data(prostateCancer)
Y <- prostateCancer[1:400,3:4]
X <- prostateCancer[1:400,5:9]
Treatment <- prostateCancer[1:400,2]
p <- pact.fit(Y=Y, X=X, Treatment=Treatment, family="cox")

newx <- prostateCancer[401:410,5:9]
newTreatment <- prostateCancer[401:410,2]
predict(p, newx, newTreatment)

### Binary response
data(EORTC10994)
Y <- EORTC10994[1:120,4]
X <- EORTC10994[1:120,c(2,5:7)]
Treatment <- EORTC10994[1:120,3]
p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="binomial")

newx <- EORTC10994[121:125,c(2,5:7)]
newTreatment <- EORTC10994[121:125,3]
predict(p, newx, newTreatment)
}
\author{
Jyothi Subramanian and Richard Simon \cr Maintainer: Jyothi
Subramanian <\email{subramanianj01@gmail.com}>
}

