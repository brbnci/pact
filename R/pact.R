#########################################################################
#								                                                        #
# Main R function for predictive analysis of clinical trials (PACT)     #
# This function handles the low dimensional case for both survival or   #
# a binary response variable.                                           #
# 				                                                              #
#                                                                       #
# See the Readme file (PACT_Readme.doc) for details on the methodology  #
# and usage details.			                                              #
# Authors: Richard Simon and Jyothi Subramanian	                        #
#			                                                                  #           
#      	                                               		              #
# Nov-8-2012		                                                        #
#   o  First version.       	                                          #
# Nov-21-2012							                                              #
#   o Second version.					                                          #
#								                                                        #
# Feb-21-2013                                                           #
#   o Third version - for Shiny app development                         #
# Nov-14-2013							                                              #
#   o Version for R package development				                          #
#########################################################################


#' @title Predictive Analysis of Clinical Trials
#'
#' @description The \code{pact} package implements a prediction-based approach to the 
#' analysis of data from randomized clinical trials (RCT). Based on clincial response and 
#' covariate data from a RCT comparing a new experimental treatment E versus a control C, 
#' the purpose behind the functions in \code{pact} is to develop and internally validate
#' a model that can identify subjects likely to benefit from E rather than C. Currently, 
#' 'survival' and 'binary' response types are permitted. 
#' 
#'@details
#' Package: pact
#' \cr Type: Package
#' \cr Version: 0.1
#' \cr Date: 2014-04-30
#' \cr Author: Dr. Jyothi Subramanian and Dr. Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' \cr License: GPL-2
#' \cr\cr \code{pact.fit} fits a predictive model to a data from a RCT. Currently, 'survival' and 
#' 'binary' response types are supported. An object of class 'pact' is returned. 
#' \code{print}, \code{summary} and \code{predict} methods are available for objects of 
#' class 'pact'.
#' The function \code{pact.cv} takes as an input the object returned by \code{pact.fit} 
#' and performs an evaluation of the model using k-fold cross-validation. 
#' Evaluations of the cross-validated predictions are performed by the function \code{eval.pact.cv}. 
#' \cr\cr  Finally, the function \code{overall.analysis}
#' also takes an object of class 'pact' as input and computes some summary statistics 
#' f the comparison of treatments E and C.
#' @docType package
#' @name pact
#' @examples
#' ### Survival response
#' set.seed(1)    
#' data(prostateCancer)     
#' Y <- prostateCancer[,3:4]
#' X <- prostateCancer[,5:9]
#' Treatment <- prostateCancer[,2]
#' p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="cox")
#' print(p)
#' overall.analysis(p)
#' cv <- pact.cv(p, nfold=5)
#' eval.pact.cv(cv, method="continuous", plot.score=TRUE, perm.test=FALSE, nperm=100)
#' 
#' ### Binary response
#' set.seed(2)
#' data(EORTC10994)
#' Y <- EORTC10994[,4]
#' X <- EORTC10994[,c(2,5:7)]
#' Treatment <- EORTC10994[,3]
#' p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="binomial")
#' print(p)
#' overall.analysis(p)
#' cv <- pact.cv(p, nfold=5)
#' eval.pact.cv(cv, method="discrete", g=log(1), perm.test=FALSE, nperm=100)
#' 
NULL

#' @title Fits a predictive model to the full dataset
#'
#' @description
#' \code{pact.fit} Fits a predictive model to the full dataset. Currently supports  
#' Cox and logistic regression models for 'survival' and 'binary' response types 
#' respectively.
#' 
#' @details
#' A Cox proportional hazards (PH) regression or a logistic regression model is 
#' developed for data with survival and binary response respectively. Data from subjects 
#' in both 'experimental' (E) and 'control' (C) groups from a RCT is used for model
#' development. Main effect of treatment, main effect of covariates and all treatment by 
#' covariate interaction terms are considered in the model development.  
#' 
#' @param Y Response variable. For \code{family='binomial'}, Y should be a factor with two 
#' levels. For \code{family='cox'}, Y should be a two-column matrix with columns named 'time' 
#' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' 
#' indicating right censored. 
#' @param X A dataframe of the predictor variables to be used for model development. 
#' All variables in the dataframe are used in the model.
#' @param Treatment The treatment assignement indicator. A factor with two levels.
#' '0' indicating control (C) and '1' indicating experimental (E) treatment.
#' @param family Type of the response variable. See above. Possible values are 'binomial'
#' or 'cox'.
#' 
#' @return An object of class 'pact' which is a list with the following components
#' @return \item{reg}{The fitted regression model}
#' @return \item{family}{The response variable type used}
#' @return \item{Y}{The response variable used}
#' @return \item{X}{The dataframe of predictor variables used}
#' @return \item{Treatment}{The treatment assignment indicator used}
#' @return \item{call}{The call that produced the return object}
#' 
#' @keywords pact, pact.fit
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' X <- prostateCancer[,5:9]
#' Treatment <- prostateCancer[,2]
#' pact.fit(Y=Y, X=X, Treatment=Treatment, family="cox")

# Fits a predictive model to the full data set and returns an object of class 'pact', 
# which is actually a list 

pact.fit <- function(Y, X, Treatment, family = c("binomial", "cox")) {
  family <- match.arg(family)
  this.call <- match.call()
  if (!inherits(X, "data.frame")) 
    stop("X should be a data frame")
 
  dimx <- dim(X)
  nrowx <- ifelse(is.null(dimx), length(X), dimx[1])
  dimy <- dim(Y)
  nrowy <- ifelse(is.null(dimy), length(Y), dimy[1])
  if (nrowx != nrowy) 
    stop(paste("number of observations in Y (", nrowy, ") not equal to the number of rows of X (", 
               nrowx, ")", sep = ""))
  nobs <- nrowy
  
  Treatment <- as.factor(Treatment)  
  if (!all(match(levels(Treatment), c(0,1),0))) 
    stop("'Treatment' should be a factor with two levels '0' and '1'. See Usage.")
  if (length(Treatment) != nobs)
    stop(paste("number of observations in Treatment (", length(Treatment), ") not equal to the number of rows of X (", 
               nrowx, ")", sep = ""))
  
  res <- switch(family, 
    cox = .pact.fit.survival(Y, X, Treatment),
    binomial = .pact.fit.binary(Y, X, Treatment))
  res$Y <- Y
  res$X <- X
  res$Treatment <- Treatment
  res$call <- this.call
  class(res) <- "pact"
  res
}
 
### Internal function. 'pact.fit' for the "cox" family.
.pact.fit.survival <- function(Y, X, Treatment) {
  if (!all(match(c("time", "status"), dimnames(Y)[[2]], 0))) 
    stop("Cox model requires a matrix with columns 'time' (>0) and 'status'  (binary) 
         as a response")
  if (any(Y[, "time"] < 0)) 
    stop("negative event times encountered; not permitted for Cox family")
  
  SurvObj <- Surv(Y[,"time"], Y[,"status"])
  Covar <- X
  dimx <- dim(X)
  nCovar <- ifelse(is.null(dimx), 1, dimx[2])
 
  T <- Treatment
  tag <- 0
  tryCatch(CoxReg <- coxph(SurvObj ~ T + . + T:., Covar, na.action="na.exclude"), 
           error=function(err) tag <<- 1)

  if (tag == 0) {
      #### Predictions (Resubstitution) for the training set, only to find the model for future application
    LinPred <- predict(CoxReg, type="lp")
      #### Now predicting for inverted treatment assignment
    T <- as.factor(ifelse((T == 1), 0, 1))
    InvLinPred <- predict(CoxReg, cbind(T,Covar), type="lp")
      #### (log HR for T=1 - log HR for T=0)
    PredScore <- ifelse((T == 0), LinPred - InvLinPred, InvLinPred - LinPred)
    missing.predscore <- sum(is.na(PredScore)) > 0	#### Any missing predscores?
  } else  { ### if Cox model failed to fit on the full dataset, no point in doing anything further
    PredScore <- NA
    stop("Predictive model cannot be developed using the given dataset")
  }
  res <- list(reg=CoxReg, family="cox")
  res
} 

### Internal function. 'pact.fit' for the "binomial" family.
.pact.fit.binary <- function(Y, X, Treatment) {  
  nc = dim(Y)

  if (is.null(nc)) {
    Y = as.factor(Y)
    ntab = table(Y)
    classnames = names(ntab)
  }
  else {
    stop("For binomial family, Y must be a factor with two levels, higher dimensions for Y currently not supported")
  }

  T <- Treatment
  Response <- Y
  Covar <- X
  dimx <- dim(X)
  nCovar <- ifelse(is.null(dimx), 1, dimx[2])
  
  tag <- 0
  tryCatch(LogReg <- glm(Response ~ T + . + T:., Covar, family = binomial(link = "logit"), na.action="na.exclude")
               , error=function(err) tag <<- 1)
  if (tag == 0) {
      ### Predictions (Resubstitution) for the training set, only to find the model and cut-point for future application
    LinPred <- predict(LogReg, type="link")
      #### Now predicting for inverted treatment assignment
    T <- as.factor(ifelse((T == 1),0,1))
    InvLinPred <- predict(LogReg, cbind(T,Covar), type="link")
      #### (log odds for T=1 - log odds for T=0)
    PredScore <- ifelse((T == 0), LinPred - InvLinPred, InvLinPred - LinPred)
    missing.predscore <- sum(is.na(PredScore)) > 0	#### Any missing pred scores?
  } else  { ### if logistic model failed to fit, no point doing anything further
    PredScore <- NA
    stop("Predictive model cannot be developed using the given dataset")
  }
  res <- list(reg=LogReg, family="binomial")
  res
}

#' @title Print an object of class 'pact' 
#'
#' @description
#' print method for objects of class 'pact'
#' 
#' @details
#' The call that produced the object is printed, followed by the classification 
#' function from \code{pact.fit} for calculating the predictive scores for new 
#' subjects
#' 
#' @method print pact
#' 
#' @param x The object returned from \code{'pact.fit'}
#' @param digits significant digits in the print 
#' @param ... Additional print arguments
#' 
#' @return The classification function is printed
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' X <- prostateCancer[,5:9]
#' Treatment <- prostateCancer[,2]
#' p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="cox")
#' print(p)

print.pact <- function(x, digits = max(3, getOption("digits") - 3), ...) {	### x = object of class 'pact'
  reg <- x$reg
  nInteractions <- (nrow(summary(reg)$coef)-1)/2

  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("\nfamily: ", x$family, "\n\n")
  if (x$family == "cox") { 
    ClassFn <- paste("(",round(summary(reg)$coef[1,1],digits),")")
    for( i in 2:(nInteractions+1)) {
      ClassFn <- paste(ClassFn, "\n", "  +  (",round(summary(reg)$coef[nInteractions+i,1],digits)," )*",row.names(summary(reg)$coef)[i],sep="")
    }
    cat(paste("\nClassification function for classifying future subjects:\n", "f = ", ClassFn, "\n\n"), sep="")
  } else if (x$family == "binomial") { ### For binary response (intercept also exists in the model output for binary response)
    ClassFn <- paste("(",round(summary(reg)$coef[1,1],digits),")")
    ClassFn <- paste(ClassFn, "\n", "  +  (",round(summary(reg)$coef[2,1],digits)," )")
    for( i in 3:(nInteractions+2)) {
      ClassFn <- paste(ClassFn, "\n", "  +  (",round(summary(reg)$coef[nInteractions+i,1],digits)," )*",row.names(summary(reg)$coef)[i],sep="")
    }
    cat(paste("\nClassification function for classifying future subjects:\n", "f = ", ClassFn, "\n\n"), sep="")
  }
}

#' @title Summarize a predictive model fit
#'
#' @description
#' summary method for objects of class 'pact'
#' 
#' @details
#' Returns the summary statistics for the regression model of the 'pact' object in a tabular 
#' format
#' 
#' @method summary pact
#' 
#' @param object The object returned from \code{'pact.fit'}
#' @param ... Additional arguments for 'summary'
#' 
#' @return The summary statistics for the regression model fitted by \code{pact.fit}
#' is returned
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' X <- prostateCancer[,5:9]
#' Treatment <- prostateCancer[,2]
#' p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="cox")
#' summary(p)
#' 
summary.pact <- function(object, ...) {	### object = object of class 'pact'
  s <- summary(object$reg)
  s
}

#' @title Predictions from a predictive model fit
#'
#' @description
#' Predicts the scores for new subjects from a previously developed object of class 'pact'
#' 
#' @details
#' Returns the scores for new subjects from an object of class 'pact', given their coveriate values and 
#' treatment assignment. 
#' 
#' @method predict pact
#' 
#' @param object The object returned from \code{pact.fit} 
#' @param newx Dataframe of predictor values for the new subjects for whom predictions 
#' are to be made
#' @param newTreatment The treatment assignment indicator for the new subjects
#' @param ... Other arguments to 'predict'
#' 
#' @return A numeric vector containing the predicted scores for the new subjects 
#' from the fitted model is returned
#'
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' ### Survival response
#' data(prostateCancer)
#' Y <- prostateCancer[1:400,3:4]
#' X <- prostateCancer[1:400,5:9]
#' Treatment <- prostateCancer[1:400,2]
#' p <- pact.fit(Y=Y, X=X, Treatment=Treatment, family="cox")
#' 
#' newx <- prostateCancer[401:410,5:9]
#' newTreatment <- prostateCancer[401:410,2]
#' predict(p, newx, newTreatment)
#' 
#' ### Binary response
#' data(EORTC10994)
#' Y <- EORTC10994[1:120,4]
#' X <- EORTC10994[1:120,c(2,5:7)]
#' Treatment <- EORTC10994[1:120,3]
#' p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="binomial")
#' 
#' newx <- EORTC10994[121:125,c(2,5:7)]
#' newTreatment <- EORTC10994[121:125,3]
#' predict(p, newx, newTreatment)

predict.pact <- function(object, newx, newTreatment, ...) {  ### object = object of class 'pact'
  if (missing(newx)) {
    stop("You need to supply a value for 'newx'")
  }
  dimx <- dim(newx)
  nobs <- ifelse(is.null(dimx), length(newx), dimx[1])
  
  if (missing(newTreatment)) {
    stop("You need to supply a value for 'newTreatment'")
  }
  T <- as.factor(newTreatment)
  if (!all(match(levels(T), c(0,1),0))) 
    stop("'newTreatment' should be a factor with two levels '0' and '1'.")
  if (length(newTreatment) != nobs)
    stop(paste("number of observations in newTreatment (", length(newTreatment), ") not equal to the number of observations in newx (", 
               nobs, ")", sep = ""))
  reg <- object$reg  
  Covar <- newx
  nCovar <- ifelse(is.null(dimx), 1, dimx[2])
  
  if (object$family == "cox") {
    LinPred <- predict(reg, cbind(T,Covar), type="lp")
        #### Now predicting for inverted treatment assignment
    T <- as.factor(ifelse((T == 1),0,1))
    InvLinPred <- predict(reg, cbind(T,Covar), type="lp")
        #### (log HR for T=1 - log HR for T=0)
    PredScore <- ifelse((T == 0), LinPred - InvLinPred, InvLinPred - LinPred)
  }
  else if (object$family == "binomial") { ### For binary response 
    LinPred <- predict(reg, type="link")
      #### Now predicting for inverted treatment assignment
    T <- as.factor(ifelse((T == 1),0,1))
    InvLinPred <- predict(reg, cbind(T,Covar), type="link")
      #### (log odds for T=1 - log odds for T=0)
    PredScore <- ifelse((T == 0), LinPred - InvLinPred, InvLinPred - LinPred)
  }
  PredScore
}

### End of main 'pact' fucntions
### Other functions

### Evaluation of an object of class 'pact' using cross-validation

#' @title Split a dataset into k parts for k-fold cross-validation
#'
#' @description
#' Split a dataset into k parts for k-fold cross-validation. This function is used in 
#' \code{pact.cv} to create the cross-validation splits
#' 
#' @param n The sample size 
#' @param k The number of folds. k=n would mean a leave-one-out cross-validation
#'  
#' @return A integer vector of same length as n. Each observation is an integer taking a value between 1 to k 
#' denoting the fold it belongs to.
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' KfoldCV(15,3)
#' KfoldCV(15,15)
 
KfoldCV <- function(n,k) { ### Partition a dataset to k parts for k-fold CV
  f <- ceiling(n/k)
  s <- sample(rep(1:k, f), n)
  s
}

#' @title Cross-validation for pact
#'
#' @description
#' k-fold cross-validation for evaluation of the model developed in \code{pact.fit}
#' 
#' @details
#' Evaluation of the model developed in \code{pact.fit} using k-fold cross-validation.
#' In each fold of the cross-validation, a model is developed from the observations in the 
#' training set. The estimated coefficients of the regression model developed using training set
#' are used to make predictions for the left out observations (test set). This is repeated for 
#' all the folds. Scores are thus obtained for all the subjects in the dataset. The function
#' \code{\link{eval.pact.cv}} provides the evaluation method options for the cross-validated 
#' scores.
#' 
#' @param p An object of class 'pact'
#' @param nfold The number of folds (k) for the k-fold cross-validation. k equal to the sample size 
#' would mean a leave-one-out cross-validation
#' 
#' @return A list with the following components
#' @return \item{PredScore}{The cross-validated scores for each subject (a vector)}
#' @return \item{Y}{The response variable used}
#' @return \item{X}{The dataframe of predictor variables used}
#' @return \item{Treatment}{The treatment assignment indicator used}
#' @return \item{family}{The response variable type}
#' @return \item{call}{The call that produced this output}
#' 
#' @keywords pact, pact.cv
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' X <- prostateCancer[,5:9]
#' Treatment <- prostateCancer[,2]
#' p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="cox")
#' cv <- pact.cv(p, nfold=5)

pact.cv <- function(p, nfold) {
  if(!inherits(p,"pact"))
    stop("First argument must be an object of class 'pact'.")
  family <- p$family
  this.call <- match.call()

  Y <- p$Y
  X <- p$X
  Treatment <- as.factor(p$Treatment)  
  
  dimx <- dim(X)
  nrowx <- ifelse(is.null(dimx), length(X), dimx[1])
  dimy <- dim(Y)
  nrowy <- ifelse(is.null(dimy), length(Y), dimy[1])
  nobs <- nrowy
  
  nfold <- as.integer(nfold)
  if ((nfold > nobs) || (nfold < 2))
    stop("nfold should be an integer >= 2 but cannot be greater than the number of observations")
  
  #### K-fold CV ####
  res <- switch(family, 
                cox = .pact.cv.survival(Y, X, Treatment, nfold),
                binomial = .pact.cv.binary(Y, X, Treatment, nfold))
  res$call <- this.call
  class(res) <- "pact.cv"
  res
}

### Internal cross-validation function. For survival response.
.pact.cv.survival <- function(Y, X, Treatment, nfold) { 
  ### Evaluation of predictive model using CV, survival response
  T <- Treatment
  SurvObj <- Surv(Y[, "time"],Y[, "status"])
  Covar <- X
  dimx <- dim(X)
  nobs <- ifelse(is.null(dimx), length(X), dimx[1])
  nCovar <- ifelse(is.null(dimx), 1, dimx[2])
  
  #### Generate IDs for K-fold cross-validation
  
  CV <- KfoldCV(nobs, nfold)
  
  #### K-fold CV
  
  Temp.CVest <- lapply(1:nfold,function(x) { 
    Ind.train <- which(CV != x)
    Ind.test <- which(CV == x)
    
    T <- Treatment[Ind.train]
    tag <- 0
    
    ##### Develop a Cox PH regression model using only the training set
    
    tryCatch(CoxReg <- coxph(SurvObj[Ind.train] ~ T + . + T:., Covar[Ind.train,,drop=FALSE], 
                             na.action="na.exclude"), error=function(err) tag <<- 1)
    if (tag == 0) {
      #### Now predicting for the test set cases using CoxReg and by repeating the treatment inversion steps
      T <- Treatment[Ind.test]
      LinPred.test <- predict(CoxReg,cbind(T,Covar[Ind.test,,drop=FALSE]),"lp")
      T <- as.factor(ifelse((T == 1),0,1))
      InvLinPred.test <- predict(CoxReg,cbind(T,Covar[Ind.test,,drop=FALSE]),"lp")
      PredScore.test <- ifelse((T == 0), LinPred.test - InvLinPred.test, InvLinPred.test - LinPred.test)
    } else { #### In case the Cox regression model failed to converge
      PredScore.test <- NA
    }
    ifelse((tag == 0), return(PredScore.test), return(NA))
  })
  
  PredScore <- NULL 
  for(i in (1:nfold)) {
    PredScore[CV == i] <- Temp.CVest[[i]]
  }
  out.cv <- list(PredScore=PredScore, Y=Y, X=X, Treatment=Treatment, family="cox", nfold=nfold)
  out.cv
}

### Internal cross-validation function. For binary response.
.pact.cv.binary <- function(Y, X, Treatment, nfold) { 
  ### Evaluation of predictive model using CV, binary response
  T <- Treatment
  Response <- Y
  Covar <- X
  dimx <- dim(X)
  nobs <- ifelse(is.null(dimx), length(X), dimx[1])
  nCovar <- ifelse(is.null(dimx), 1, dimx[2])
  
  #### Generate IDs for K-fold cross-validation
  
  CV <- KfoldCV(nobs, nfold)

  #### K-fold CV
  
  Temp.CVest <- lapply(1:nfold,function(x) { 
    Ind.train <- which(CV != x)
    Ind.test <- which(CV == x)

    T <- Treatment[Ind.train]
    tag <- 0
    
    ##### Develop a Cox PH regression model using only the training set
    tryCatch(LogReg <- glm(Response[Ind.train] ~ T + . + T:., Covar[Ind.train,,drop=FALSE], 
                           family = binomial(link = "logit"), na.action="na.exclude"), error=function(err) tag <<- 1)
 
    if (tag == 0) {
      #### Now predicting for the test set cases using CoxReg and by repeating the treatment inversion steps
      T <- Treatment[Ind.test]
      LinPred.test <- predict(LogReg, cbind(T,Covar[Ind.test,,drop=FALSE]),type="link")
      T <- as.factor(ifelse((T == 1),0,1))
      InvLinPred.test <- predict(LogReg, cbind(T,Covar[Ind.test,,drop=FALSE]),"link")
      PredScore.test <- ifelse((T == 0), LinPred.test - InvLinPred.test, InvLinPred.test - LinPred.test)
    } else { #### In case the Cox regression model failed to converge
      PredScore.test <- NA
    }
    ifelse((tag == 0), return(PredScore.test), return(NA))
  })
  
  PredScore <- NULL 
  for(i in (1:nfold)) {
    PredScore[CV == i] <- Temp.CVest[[i]]
  }
  out.cv <- list(PredScore=PredScore, Y=Y, X=X, Treatment=Treatment, family="binomial", nfold=nfold)
  out.cv
}


#' @title Evaluation functions for cross-validated predictions
#'
#' @description
#' Methods for the evaluation of the cross-validated predictive scores obtained from \code{\link{pact.cv}}
#' 
#' @details
#' Currently two methods are defined for the evaluation of the scores obtained from \code{pact.cv}. In 
#' \code{method='discrete'} a user specified cut-off score is used to classify the subjects into groups
#' 'benefit' or 'do not benefit' from new treatment. In each of the 'benefit' and 'do not benefit' groups 
#' the actual responses in the control (C) and the experimental (E) groups are compared. For the 'cox' family,  
#' the 'score' for a subject represents the predicted change in the log hazard when 
#' the subject is treated with E as against C (with lower values denoting benefit with E). In the case of the 
#' 'binomial' family, the 'score' represents the predicted change in the log odds of a response when the 
#' subject is treated with E as against C (with higher values denoting benefit with E). 
#' For the 'cox' family, examples of the cut-point \code{g} could be \code{g=log(1)} with score < g
#' meaning benefit with E. Or one could be more stringent and have \code{g} correspond to 
#' a 30\% reduction in hazard (\code{g=log(0.70)}). 
#' For the 'binomial' family, \code{g=log(1.20)} with score > g meaning sensitive to E would mean that
#' subjects predicted to receive at least 20\% increase in log odds of response with E are 
#' classified as benefitting from E.  
#' \cr\cr In \code{method='continuous'} the cross-validated scores are kept continuous. 
#' At the end of cross-validation, a Cox proportional hazards (PH) regression or a logistic regression 
#' model (respectively for 'survival' and 'binary' response) is developed that includes 
#' the main effect of treatment, main effect of cross-validated score, and treatment*score interaction.
#' For survival response, this model is used to generate the Kaplan Meier survival curves for each treatment 
#' at the at 20\%, 40\%, 60\% and 80\% quantiles of predictive scores (\code{plot.score = TRUE}). 
#' The model is also used to compute the estimated probability of surviving beyond a landmark time 
#' specified in \code{plot.time} as a function of treatment and (cross-validated) score (if 
#' \code{plot.time = NULL}, this plot is not produced). For binary response, the output from evaluation
#' is a plot of the probability of response as a functions of the predictive score and Treatment.
#' \cr\cr If \code{perm.test=TRUE}, permutation based significance tests are performed on appropriate 
#' test statistics and p-values are computed. See 'Value' and the package vignette and for more details 
#' on the permutation tests.
#'
#' @param out.cv The object from \code{pact.cv} 
#' @param method The evaluation method. Currently two options, \code{method='discrete'} or
#' \code{method='continuous'}, are available. See 'Details'.
#' @param g The cut-point for grouping scores into subsets 'benefit' and 'no benefit' from 
#' new treatment. Ignored for \code{method='continuous'}.
#' @param plot.score Used only for plots if \code{method='continuous'} is chosen for survival response. 
#' Logical representing whether survival 
#' curves at specific quantiles of cross-validated scores are to be drawn. See 'Details'.
#' @param plot.time Used only for plots if \code{method='continuous'} is chosen for survival response. 
#' Probability of survival greater than \code{plot.time} is plotted as a function of cross-validated 
#' score and Treatment. See 'Details'.
#' @param perm.test Logical. If \code{perm.test=TRUE}, a permutation based test of significance
#' is conducted for statistics computed from the cross-validated scores. See 'Value' and the package 
#' vignette and for more details on the permutation tests.
#' @param nperm The number of permutations for the permutation test. Ignored if \code{perm.test=FALSE}
#' 
#' @return The return object is a list whose components depend on the family ('cox' or 'binomial') and the 
#' chosen evaluation method ('continuous' or 'discrete')
#' @return \item{LR.Benefit}{For \code{family='cox'} and \code{method='discrete'}. The log-rank statistic for the 
#' survival difference between E and C for the 'benefit' from E group.}
#' @return \item{LR.NoBenefit}{For \code{family='cox'} and \code{method='discrete'}. The log-rank statistic for the 
#' survival difference between E and C for the 'do not benefit' from E group.}
#' 
#' @return \item{RR.T.Benefit}{For \code{family='binomial'} and \code{method='discrete'}. 
#' The response rate for subjects getting E in the 'benefit' from E group.}
#' @return \item{RR.C.Benefit}{For \code{family='binomial'} and \code{method='discrete'}. 
#' The response rate for subjects getting C in the 'benefit' from E group}
#' @return \item{RR.T.NoBenefit}{For \code{family='binomial'} and \code{method='discrete'}. 
#' The response rate for subjects getting E in the 'do not benefit' from E group}
#' @return \item{RR.C.NoBenefit}{For \code{family='binomial'} and \code{method='discrete'}. 
#' The response rate for subjects getting C in the 'do not benefit' from E group}.
#' 
#' @return \item{pval.Benefit}{If \code{perm.test=TRUE}, p-value from permutation test. 
#' For \code{family='cox'} and \code{method='discrete'}, permutation based p-value for LR.Benefit. 
#' For \code{family='binomial'} and \code{method='discrete'}, permutation based p-value for 
#' difference in response rates for E and C for the subset predicted 'benefit' from E.}
#' @return \item{pval.NoBenefit}{If \code{perm.test=TRUE}, p-value from permutation test. 
#' For \code{family='cox'} and \code{method='discrete'}, permutation based p-value for LR.NoBenefit. 
#' For \code{family='binomial'} and \code{method='discrete'}, permutation based p-value for 
#' difference in response rates for E and C for the subset predicted 'no benefit' from E.} 
#' 
#' @return \item{reg}{For \code{method='continuous'}, the regression model with treatment, predictive score and
#' treatment x predictive score interaction}
#' @return \item{pval.twosided}{For \code{method='continuous'}. Two-sided (non-directional) permutation based p-value for the treatment x predictive score
#' interaction coefficient}
#' @return \item{pval.onesided}{For \code{method='continuous'}. One-sided (directional, greater) permutation based p-value for the treatment x predictive score
#' interaction coefficient}
#' 
#' @return \item{call}{The function call}
#' 
#' @return Additional plots for both \code{method='discrete'} as well as \code{method='continuous'}. 
#' See package vignette.
#' 
#' 
#' @keywords pact, pact.cv
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' ### Survival response
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' X <- prostateCancer[,5:9]
#' Treatment <- prostateCancer[,2]
#' p <- pact.fit(Y=Y, X=X, Treatment=Treatment, family="cox")
#' cv <- pact.cv(p, nfold=5)
#' \dontrun{eval.pact.cv(cv, method="discrete", g=log(0.80), perm.test=TRUE, nperm=100)}  ## At least 20% predicted reduction in HR classified as 'sensitive'
#' eval.pact.cv(cv, method="continuous", plot.score=TRUE, perm.test=FALSE, nperm=100)
#' 
#' ### Binary response
#' data(EORTC10994)
#' Y <- EORTC10994[,4]
#' X <- EORTC10994[,c(2,5:7)]
#' Treatment <- EORTC10994[,3]
#' p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="binomial")
#' cv <- pact.cv(p, nfold=5)
#' eval.pact.cv(cv, method="discrete", g=log(1), perm.test=TRUE, nperm=100)
#' 
#' 
### Evaluation functions for cross-validation

eval.pact.cv <- function(out.cv, method=c("discrete","continuous"), g=log(1), 
                         plot.score=TRUE, plot.time=NULL, perm.test=FALSE, nperm=100) {
  if(!inherits(out.cv,"pact.cv"))
    stop("First argument must be an object of class 'pact.cv'. Run pact.cv first.")
  method <- match.arg(method)
  family <- out.cv$family
  this.call <- match.call()
  if (family == "cox") {
    if (method == "discrete") {
      if (g < 0)
        cat(paste("\n\n Finding subset of subjects predicted to benefit from E with atleast \n",(1-exp(g))*100,
                  " % reduction in HR\n\n", sep=""))
      if (g > 0)
        cat(paste("\n\n Finding subset of subjects predicted no benefit from E with atleast \n",(exp(g)-1)*100,
                    " % increase in HR\n\n", sep=""))
      eval.cv <- .eval.pact.cv.survival.discrete(out.cv, g, perm.test, nperm)
    }
    else {
      if (method == "continuous") 
        eval.cv <- .eval.pact.cv.survival.continuous(out.cv, plot.score, plot.time, perm.test, nperm)
      else 
        stop("'method' can only be 'discrete' or 'continuous'")
    }
  }
  else {
    if (family == "binomial") {
      if (method == "discrete") {
        if (g > 0)
          cat(paste("\n\n Finding subset of subjects predicted to benefit from E with atleast \n",(exp(g)-1)*100,
                    " % increase in odds of response\n\n", sep=""))
        if (g <= 0)
          cat(paste("\n\n Finding subset of subjects predicted no benefit from E with atleast \n",(1-exp(g))*100,
                    " % decrease in odds of response\n\n", sep=""))
        eval.cv <- .eval.pact.cv.binary.discrete(out.cv, g, perm.test, nperm)
      }
      else {
        if (method == "continuous") 
          eval.cv <- .eval.pact.cv.binary.continuous(out.cv, perm.test, nperm)
        else 
          stop("'method' can only be 'discrete' or 'continuous'")
      }
    }
    else
      stop("Error: family can only be 'cox' or 'binomial'")
  }
  eval.cv$call <- this.call
  eval.cv
}


### Internal cross-validation evaluation function. For survival response. Discretized scores.
.eval.pact.cv.survival.discrete <- function(out.cv, g, perm.test, nperm) {
  predscore <- out.cv$PredScore   ### predscore represents change in log HR (E to C)
  strata <- ifelse(predscore < g, 1, 0) ### 1=(predicted) benefit from T, 0=(predicted) no benefit from E
  Y <- out.cv$Y
  SurvObj <- Surv(Y[, "time"],Y[, "status"])
  T <- out.cv$Treatment
  LR.Benefit <- survdiff(SurvObj ~ T, subset=(strata == 1))$chisq  ### LR statistic between control and new treatment groups for cases predicted to benefit from E
  LR.NoBenefit <- survdiff(SurvObj ~ T,subset=(strata == 0))$chisq  ### LR statistic between control and new treatment groups for cases predicted to not benefit from E
  
  par(mfrow=c(1,2), font.lab=2, cex.lab=1, cex.main=0.8)
  sf <- survfit(SurvObj ~ T, subset=(strata == 1)) ### KM plot by treatment for cases predicted to benefit from new
  plot(sf, col=c("red","blue"), lwd=2, main="Subset predicted benefit from E",
       xlab="Time", ylab="Proportion alive")
  legend("bottomleft", c("Standard Treatment", "New Treatment"),col=c("red","blue"),lty=1,lwd=2,bty="n",cex=0.8)

  sf <- survfit(SurvObj ~ T, subset=(strata == 0)) ### KM plot by treatment for cases predicted to not benefit from new
  plot(sf, col=c("red","blue"), lwd=2, main="Subset predicted no benefit from E",
       xlab="Time", ylab="Proportion alive")
  legend("bottomleft", c("Standard Treatment", "New Treatment"),col=c("red","blue"),lty=1,lwd=2,bty="n",cex=0.8)
  
  if (perm.test) {  ### if permutation testing of LR statistic is desired
    X <- out.cv$X
    dimx <- dim(X)
    nobs <- ifelse(is.null(dimx), length(X), dimx[1])
    nCovar <- ifelse(is.null(dimx), 1, dimx[2])
    nfold <- out.cv$nfold

    permute.LR.Benefit <- vector("numeric", nperm)
    permute.LR.NoBenefit <- vector("numeric", nperm)
    fail.count1 <- 0
    fail.count2 <- 0
    cat("Start Permutations \n\n")
    for(i in (1:nperm)) {
      #### Generate a permutation of the treatment labels
      permute.T <- T[sample(nobs,nobs)]
      #### Cross-validated predictions for the permuted data
      permute.PredTmnt.CV <- .pact.cv.survival(Y, X, permute.T, nfold)
      permute.strata <- ifelse(permute.PredTmnt.CV$PredScore < g, 1, 0)  
      tryCatch(permute.LR.Benefit[i] <- survdiff(SurvObj ~ permute.T, subset=(permute.strata == 1))$chisq, error=function(err) fail.count1 <<- fail.count1+1)
      tryCatch(permute.LR.NoBenefit[i] <- survdiff(SurvObj ~ permute.T, subset=(permute.strata == 0))$chisq, error=function(err) fail.count2 <<- fail.count2+1)
    }
    cat("End Permutations \n\n")
    ### Very small or very large values of 'g' would result in too few subjects in a strata and in turn to a failed survfit model
    if (fail.count1*100/nperm > 5)  
      warning("Permutation results may not be valid for the subset predicted to benefit: g may be too small")
    pval.LR.Benefit <- (1+sum((permute.LR.Benefit > LR.Benefit), na.rm=TRUE))/(1+nperm-fail.count1)
    if (fail.count2*100/nperm > 5)
      warning("Permutation results may not be valid for the subset predicted to not benefit: g may be too large")
    pval.LR.NoBenefit <- (1+sum((permute.LR.NoBenefit > LR.NoBenefit), na.rm=TRUE))/(1+nperm-fail.count2)
    eval.cv.res <- list(LR.Benefit=LR.Benefit,LR.NoBenefit=LR.NoBenefit,
                        pval.Benefit=pval.LR.Benefit,pval.NoBenefit=pval.LR.NoBenefit)
  } else ### No permulations, no p-values
    eval.cv.res <- list(LR.Benefit=LR.Benefit,LR.NoBenefit=LR.NoBenefit)
  eval.cv.res  
}

### Internal cross-validation evaluation function. For binary response. Discretized scores.
.eval.pact.cv.binary.discrete <- function(out.cv, g, perm.test, nperm) {
  predscore <- out.cv$PredScore   ### predscore represents change in log odds (E to C)
  strata <- ifelse(predscore > g, 1, 0) ### 1=(predicted) benefit from E, 0=(predicted) no benefit from E
  Response <- out.cv$Y
  T <- out.cv$Treatment
  Sens.subset <- strata == 1
  NotSens.subset <- strata == 0
  
  #### Difference in response rate between control and new treatment groups for the 
  #### sensitive subset (subset predicted to benefit from E)
  
  Response.T.Benefit <- sum(Sens.subset & Response == 1 & T == 1, na.rm=TRUE)
  n.T.Benefit <- sum(Sens.subset & T == 1, na.rm=TRUE)
  Response.C.Benefit <- sum(Sens.subset & Response == 1 & T == 0, na.rm=TRUE)
  n.C.Benefit <- sum(Sens.subset & T == 0, na.rm=TRUE)
  test.Benefit <- prop.test(x=c(Response.T.Benefit, Response.C.Benefit), n=c(n.T.Benefit,n.C.Benefit))
  Benefit.teststat <- test.Benefit$statistic
  
  #### Difference in response rate between control and new treatment groups for the 
  #### not sensitive subset (subset predicted to not benefit from E)
  
  Response.T.NoBenefit <- sum(NotSens.subset & Response == 1 & T == 1, na.rm=TRUE)
  n.T.NoBenefit <- sum(NotSens.subset & T == 1, na.rm=TRUE)
  Response.C.NoBenefit <- sum(NotSens.subset & Response == 1 & T == 0, na.rm=TRUE)
  n.C.NoBenefit <- sum(NotSens.subset & T == 0, na.rm=TRUE)
  
  test.NoBenefit <- prop.test(x=c(Response.T.NoBenefit, Response.C.NoBenefit), n=c(n.T.NoBenefit,n.C.NoBenefit))
  NoBenefit.teststat <- test.NoBenefit$statistic 
  
  if (perm.test) {  ### if permutation testing of test statistics are desired
    X <- out.cv$X
    dimx <- dim(X)
    nobs <- ifelse(is.null(dimx), length(X), dimx[1])
    nCovar <- ifelse(is.null(dimx), 1, dimx[2])
    nfold <- out.cv$nfold
    
    permute.teststat.Benefit <- vector("numeric", nperm)
    permute.teststat.NoBenefit <- vector("numeric", nperm)
    fail.count1 <- 0
    fail.count2 <- 0
    cat("Start Permutations \n\n")
    for(i in (1:nperm)) {
      #### Generate a permutation of the treatment labels
      permute.T <- T[sample(nobs,nobs)]
      #### Cross-validated predictions for the permuted data
      permute.PredTmnt.CV <- .pact.cv.binary(Response, X, permute.T, nfold)
      permute.strata <- ifelse(permute.PredTmnt.CV$PredScore > g, 1, 0) 
      
      permute.sens.subset <- permute.strata == 1
      permute.NotSens.subset <- permute.strata == 0
      
      Response.T.Benefit <- sum(permute.sens.subset & Response == 1 & permute.T == 1,na.rm=TRUE)
      n.T.Benefit <- sum(permute.sens.subset & permute.T == 1,na.rm=TRUE)
      Response.C.Benefit <- sum(permute.sens.subset & Response == 1 & permute.T == 0,na.rm=TRUE)
      n.C.Benefit <- sum(permute.sens.subset & permute.T == 0,na.rm=TRUE)
      
      test <- prop.test(x=c(Response.T.Benefit, Response.C.Benefit), n=c(n.T.Benefit,n.C.Benefit))
      permute.teststat.Benefit[i] <- test$statistic ## the value of Pearson's chi-squared test statistic.
      
      Response.T.NoBenefit <- sum(permute.NotSens.subset & Response == 1 & permute.T == 1,na.rm=TRUE)
      n.T.NoBenefit <- sum(permute.NotSens.subset & permute.T == 1,na.rm=TRUE)
      Response.C.NoBenefit <- sum(permute.NotSens.subset & Response == 1 & permute.T == 0,na.rm=TRUE)
      n.C.NoBenefit <- sum(permute.NotSens.subset & permute.T == 0,na.rm=TRUE)
      
      test <- prop.test(x=c(Response.T.NoBenefit, Response.C.NoBenefit), n=c(n.T.NoBenefit, n.C.NoBenefit))
      permute.teststat.NoBenefit[i] <- test$statistic ## the value of Pearson's chi-squared test statistic.
    }
    cat("End Permutations \n\n")
    pval.Benefit <- (1+sum(permute.teststat.Benefit > Benefit.teststat))/(1+nperm)
    pval.NoBenefit <- (1+sum(permute.teststat.NoBenefit > NoBenefit.teststat))/(1+nperm)
    eval.cv.res <- list(RR.E.Benefit=test.Benefit$estimate[1], RR.C.Benefit=test.Benefit$estimate[2],
                        RR.E.NoBenefit=test.NoBenefit$estimate[1], RR.C.NoBenefit=test.NoBenefit$estimate[2],
                        pval.Benefit=pval.Benefit, pval.NoBenefit=pval.NoBenefit)
  } else ### No permutations, no p-values
    eval.cv.res <- list(RR.E.Benefit=test.Benefit$estimate[1], RR.C.Benefit=test.Benefit$estimate[2],
                        RR.E.NoBenefit=test.NoBenefit$estimate[1], RR.C.NoBenefit=test.NoBenefit$estimate[2])
  eval.cv.res  
}

.eval.pact.cv.survival.continuous <- function(out.cv, plot.score, plot.time, perm.test, nperm) {
  predscore <- out.cv$PredScore   ### predscore represents change in log HR (E to C)
  Y <- out.cv$Y
  SurvObj <- Surv(Y[, "time"],Y[, "status"])
  T <- out.cv$Treatment
  nobs <- length(T)
  
  CoxReg <- coxph(SurvObj ~ T + predscore + T*predscore, na.action="na.exclude")

  if ((plot.score)) { ## plot at 20%, 40%, 60% and 80% quantiles of scores
    score <- quantile(predscore, probs=c(0.2,0.4,0.6,0.8), na.rm=TRUE)
    lenq <- 4
### one plot per quantile as function of time
    newdata <- data.frame(T=as.factor(c(rep(c(1,0),lenq))), predscore=rep(score,each=2))
    nr <- 2   ### no. of graph rows
    mat <- matrix(c(1:4,5,5),nrow=3,ncol=2,byrow=TRUE)
    layout(mat=mat, heights = c(rep.int(0.5,2),0.1))
    par(mar=c(3, 4, 4, 2)+0.1, font.lab=2, cex.lab=0.9, cex.main=0.9, cex.axis=0.9)
    sf <- survfit(CoxReg, newdata=newdata) 
    for (i in seq(1, (2L*lenq), 2)) {
      plot(sf[c(i:(i+1))], lwd=2, main=paste("Score = ",round(newdata[i,2],2),sep=""), 
         xlab="Time", ylab="Proportion alive",col=c("blue","red"))
    }
    par(mai=c(0,0,0,0))   ### dummy plot for legend
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("center",legend=c("New Treatment (E)","Control (C)"),col=c("blue","red"),
         lty=1,lwd=2,xpd=NA,bty="n")
  }
  
  max.time <- max(Y[, "time"], na.rm=TRUE)
  if (!is.null(plot.time)){
    if (plot.time > max.time) 
      stop("Requested 'time' for plot exceeds maximum follow-up time in dataset")
    else {
      #dev.new()
      par(mfrow=c(1,1), mar=c(4, 4, 4, 5)+0.1, font.lab=2, cex.lab=0.8, cex.main=0.9, cex.axis=0.8)
      q <- seq(quantile(predscore,0.005,na.rm=TRUE), quantile(predscore,0.995,na.rm=TRUE), length.out=200)   ### generate a grid of score values
      lenq <- length(q)
      newdata <- data.frame(T=as.factor(c(rep(c(1,0),lenq))), predscore=rep(q,each=2))
      sf <- summary(survfit(CoxReg, newdata=newdata))
      prob.E <- unlist(lapply(1:lenq, function(x) {
        sfn.E <- stepfun(sf$time, c(1, sf$surv[,(2*x-1)]), right=TRUE, f=1)
        val.E <- sfn.E(plot.time)
        val.E
        }))
      prob.C <- unlist(lapply(1:lenq, function(x) {
        sfn.C <- stepfun(sf$time, c(1, sf$surv[,(2*x)]), right=TRUE, f=1)
        val.C <- sfn.C(plot.time)
        val.C
        }))
      plot(prob.E ~ q, type="l", lwd=2, col="blue", xlab="Score", 
           ylab=paste("P[Survival > ",plot.time,"]",sep=""),
           main=paste("Probability of surviving beyond landmark time"),
           ylim=c(0,1))
      lines(prob.C ~ q, lwd=2, col="red")
      legend("bottomright",legend=c("New Treatment (E)","Control (C)"),col=c("blue","red"),
             lty=1,lwd=2,xpd=NA,bty="n",cex=0.8)
    }
  }

  if (perm.test) {  ### if permutation testing is desired
    intercoef <- CoxReg$coef[3] ### The true interaction coeff
    permute.intercoef <- vector("numeric", nperm)
    cat("Start Permutations\n\n")
    for(i in (1:nperm)) {
      #### Generate a permutation of the treatment labels
      permute.T <- T[sample(nobs,nobs)]
      #### Model for permuted data
      permute.CoxReg <- coxph(SurvObj ~ permute.T + predscore + permute.T*predscore, na.action="na.exclude")
      permute.intercoef[i] <- permute.CoxReg$coef[3]
    }
    pval.twosided <- (1+sum(abs(permute.intercoef) > abs(intercoef), na.rm=TRUE))/(1+nperm)    
    pval.onesided <- (1+sum(permute.intercoef > intercoef, na.rm=TRUE))/(1+nperm)
    cat("End Permutations\n\n")
  }
  if (perm.test)  
    eval.cv.res <- list(reg=summary(CoxReg),pval.twosided=pval.twosided,pval.onesided=pval.onesided)
  else ## No permutations, no p-values
    eval.cv.res <- list(reg=summary(CoxReg))
  eval.cv.res  
}

.eval.pact.cv.binary.continuous <- function(out.cv, perm.test, nperm) {
  predscore <- out.cv$PredScore   ### predscore represents change in log odds (E to C)
  Response <- out.cv$Y
  T <- out.cv$Treatment
  nobs <- length(T)
  
  LogReg <- glm(Response ~ T + predscore + T*predscore, 
                family = binomial(link = "logit"), na.action="na.exclude")
  
  ### Probabiliy plot
  par(mfrow=c(1,1), mar=c(4, 4, 4, 5)+0.1, font.lab=2, cex.lab=0.8, cex.main=0.8, cex.axis=0.8)
    ### generate a grid of score values
  q <- seq(quantile(predscore,0.005,na.rm=TRUE), quantile(predscore,0.995,na.rm=TRUE), length.out=200)   
  lenq <- length(q)
  newdata.E <- data.frame(T=as.factor(rep(1,lenq)), predscore=q)
  newdata.C <- data.frame(T=as.factor(rep(0,lenq)), predscore=q)
    ### Need to get prob.E and prob.C (prob of response with E and C respectively) for newdata
  prob.E <- predict(LogReg,newdata.E,type="response")
  prob.C <- predict(LogReg,newdata.C,type="response")
  plot(prob.E ~ q, type="l", lwd=2, col="blue", xlab="Score", ylab=paste("Probability of response"),
       main=paste("Probability of Response"), ylim=c(0,1))
  lines(prob.C ~ q, lwd=2, col="red")
  legend("bottomright",legend=c("New Treatment (E)","Control (C)"),col=c("blue","red"),
             lty=1,lwd=2,xpd=NA,bty="n",cex=0.8)
  
  if (perm.test) {  ### if permutation testing is desired
    intercoef <- LogReg$coef[4] ### The true interaction coeff
    permute.intercoef <- vector("numeric", nperm)
    cat("Start Permutations \n\n")
    for(i in (1:nperm)) {
      #### Generate a permutation of the treatment labels
      permute.T <- T[sample(nobs,nobs)]
      #### Model for permuted data
      permute.LogReg <- glm(Response ~ permute.T + predscore + permute.T*predscore, 
                            family = binomial(link = "logit"), na.action="na.exclude")
      
      permute.intercoef[i] <- permute.LogReg$coef[4]
    }
    pval.twosided <- (1+sum(abs(permute.intercoef) > abs(intercoef), na.rm=TRUE))/(1+nperm)    
    pval.onesided <- (1+sum(permute.intercoef > intercoef, na.rm=TRUE))/(1+nperm)
    cat("End Permutations \n\n")
  }
  if (perm.test)  
    eval.cv.res <- list(reg=summary(LogReg),pval.twosided=pval.twosided,pval.onesided=pval.onesided)
  else
    eval.cv.res <- list(reg=summary(LogReg))
  eval.cv.res  
}

#' @title Overall statistics and inference 
#'
#' @description
#' Produces some statistics for the overall (non-predictive) comparison of the E and C 
#' for the same dataset for which the predictive model \code{pact.fit} was developed. 
#' 
#' @details
#' Statistics for the overall comparison of the E and C is produced for the the 
#' data from a randomized clinical trial. The input is an object of class \code{pact}. 
#' 
#' @param p An object of class \code{pact}
#' 
#' @return An list with the following components. As a side effect, these are also printed
#' on screen
#' 
#' @return \item{family}{The response variable type used}
#' @return \item{nobs}{The sample size}
#' @return \item{n.E}{Number of subjects getting the new treatment}
#' @return \item{n.C}{Number of subjects getting the control treatment}
#' @return \item{nCovar}{The number of predictor variables}
#' @return \item{Covar.names}{The names of predictor variables}
#' @return \item{LR}{The log-rank statistic for the overal difference in survival between
#' E and C groups (for family="cox")}
#' @return \item{LR.pval}{The p-value for LR based on the log-rank test (for family="cox")}
#' @return \item{RR.E}{The response rate for group treated with E (new treatment) (for family="binomial")}
#' @return \item{RR.C}{The response rate for group treated with C (Control) (for family="binomial")}
#' @return \item{RRdiff.pval}{The chi-square test based pvalue for the difference in response rates 
#' (for family="binomial")}
#' 
#' @keywords pact
#' @export
#' @author Jyothi Subramanian and Richard Simon
#' \cr Maintainer: Jyothi Subramanian <\email{subramanianj01@@gmail.com}>
#' @examples
#' data(prostateCancer)
#' Y <- prostateCancer[,3:4]
#' X <- prostateCancer[,5:9]
#' Treatment <- prostateCancer[,2]
#' p <- pact.fit(Y=Y, X=X, Treatment=Treatment, family="cox")
#' overall.analysis(p)
#' 
#' ### Binary response
#' data(EORTC10994)
#' Y <- EORTC10994[,4]
#' X <- EORTC10994[,c(2,5:7)]
#' Treatment <- EORTC10994[,3]
#' p <- pact.fit(Y=Y,X=X,Treatment=Treatment,family="binomial")

### Function for overall analysis

overall.analysis <- function(p) { ### p=an object of class 'pact'  
  if(!inherits(p,"pact"))
  stop("Argument must be an object of class 'pact'.")
  
  Covar <- p$X
  Y <- p$Y
  T <- p$Treatment
  family <- p$family
  
  dimx <- dim(Covar)
  nrowx <- ifelse(is.null(dimx), length(Covar), dimx[1])
  nCovar <- ifelse(is.null(dimx), 1, dimx[2])
  covar.names <- paste(colnames(Covar), collapse=", ")
  
  dimy <- dim(Y)
  nrowy <- ifelse(is.null(dimy), length(Y), dimy[1])
  nobs <- nrowy
  
  n.T <- sum(T == 1,na.rm=TRUE)
  n.C <- sum(T == 0,na.rm=TRUE)
  
  if (family == "binomial") {
    Response.T <- sum(Y == 1 & T == 1,na.rm=TRUE)
    Response.C <- sum(Y == 1 & T == 0,na.rm=TRUE)
    Overall.test <- prop.test(x=c(Response.T, Response.C), n=c(n.T,n.C))
    RR.T <- Overall.test$estimate[1]
    RR.C <- Overall.test$estimate[2]
    pval <- Overall.test$p.value
    output <- paste("Description of the problem:\n",
                    "--------------------------\n\n",
                    "family: ", family,
                    "\n\nNumber of subjects: ", nobs,
                    "\n\nNumber of subjects assigned new treatment: ", n.T,
                    "\n\nNumber of subjects assigned standard treatment: ", n.C,
                    "\n\nNumber of covariates: ", nCovar,
                    "\n\nNames of covariates used in the analyis: ",covar.names, 
                    "\n\nOverall Analysis \n",
                    "--------------------\n\n",
                    "Response rate in the group administered new treatment: ",
                    round(RR.T,2),
                    "\n\nResponse rate in the group administered standard treatment: ",
                    round(RR.C,2),
                    "\n\n(Chi-square test based) P-value for the test of significance of the difference in response rates: ",
                    round(pval,4),
                    "\n\n",
                    sep="")
    cat(output)
    outlist <- list(family=family, nobs=nobs, n.E=n.T, n.C=n.C, nCovar=nCovar, 
                    covar.names=covar.names, RR.E=RR.T, RR.C=RR.C, RRdiff.pval=pval)
  } else {
    SurvObj <- Surv(Y[, "time"],Y[, "status"])
    sf <- survfit(SurvObj ~ Treatment)
    LR.full <- survdiff(SurvObj ~ Treatment)$chisq
    LR.pval <- 1 - pchisq(LR.full, 1)
    
    output <- paste("\nDescription of the problem:\n",
                    "----------------------------\n\n",
                    "family: ", family,
                    "\n\nNumber of subjects: ", nobs,
                    "\n\nNumber of subjects assigned new treatment: ", n.T,
                    "\n\nNumber of subjects assigned standard treatment: ", n.C,
                    "\n\nNumber of covariates: ", nCovar,
                    "\n\nNames of covariates used in the analyis: ",covar.names, 
                    "\n\nOverall Analysis:\n",
                    "--------------------\n\n",
                    "Log-rank statistic: ", round(LR.full,3),
                    "\n\np-value (Log-rank test): ", round(LR.pval,2),
                    "\n\nOverall Kaplan-Meier plot in the plot window \n\n",sep="")
    cat(output)
    
    par(font.lab=2, cex.lab=1.2)
    plot(sf, col=c("red","blue"), lwd=2, main="Overall Analysis",
         xlab="Time", ylab="Proportion alive")
    legend("bottomleft", c("Standard Treatment", "New Treatment"), col=c("red","blue"), lty=1, lwd=2, bty="n", cex=1)
    outlist <- list(family=family, nobs=nobs, n.E=n.T, n.C=n.C, nCovar=nCovar, 
                    covar.names=covar.names, LR=LR.full, LR.pval=LR.pval)
  }
  ##outlist
}

#' Prostate cancer dataset
#' 
#' A dataset containing survival information for 485 subjects with prostate cancer. 
#' See 'Details' for the variables in this dataset.
#' 
#' \itemize{
#'   \item ID. Subject identifier 
#'   \item Treatment. Treatment received. '0' for control and '1' for new
#'   \item time. Survival time in months
#'   \item status. Censoring status. '0' for censored, '1' for died
#'   \item age. Age (in years) at diagnosis
#'   \item pf. Performance status (Normal Activity or Limited Activity)
#'   \item sz. Size of the primary tumor (cm2)
#'   \item sg. Index of a combination of tumor stage and histologic grade 
#'   \item ap. Serum phosphatic acid phosphatase levels
#'   }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 485 rows and 9 variables
#' @name prostateCancer
NULL

#' EORTC10994 dataset
#' 
#' A dataset containing treatment, response and covariate information for 125 subjects 
#' with breast cancer. See 'Details' for the variables in this dataset.
#' 
#' \itemize{
#'   \item ID. Subject identifier 
#'   \item Treatment. Treatment received. '0' for control and '1' for experimental
#'   \item Response. Binary response to treatment - '0' for 'non-responder' and '1' for 'responder'
#'   \item Age. Age (in years) at diagnosis
#'   \item TumorSize. size of tumor
#'   \item Node. Node positive (Yes) or node negative (No)
#'   \item ERBB2Log2. log2 transformed ERBB2 levels
#'   }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 125 rows and 7 variables
#' @name EORTC10994
NULL