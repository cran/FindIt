\name{FindIt}
\alias{FindIt}
\alias{summary.FindIt}


\title{FindIt for Estimating Heterogeneous Treatment Effects}
\description{
  \code{FindIt} takes a binary or continuous outcome and a set of pre-treatment covariates, and returns a model with the most predictive treatment-covariate interactions returned.
}
\usage{
FindIt(model.treat, model.main, model.int, data = NULL, search.lambdas = TRUE,
lambdas = NULL, wts = 1, scale.c = 1, scale.int = 1, fit.glmnet = TRUE,
make.twoway = TRUE, make.allway = TRUE, reference.main = NULL, type = "binary",
treat.type = "multiple", threshold = 0.999999, make.reference = TRUE, nway = 2)
}
\arguments{
  \item{model.treat}{The formula that specifies outcome and treatment variables.}
  \item{model.main}{The formula that specifies pre-treatment covariates.}
  \item{model.int}{The formula specifying pre-treatment covariates that are interacted with treatment assignments when \code{treat.type="single"}.}	
  \item{data}{an optional data frame, list or environment (or object
              coercible by 'as.data.frame' to a data frame) containing the
              variables in the model.  If not found in 'data', the
              variables are taken from 'environment(formula)', typically
              the environment from which 'lm' is called.}

  \item{search.lambdas}{Whether to search for the tuning parameters for the LASSO constraints.  If FALSE, \code{lambdas} must be supplied.}
  \item{lambdas}{Tuning parameters to be given to \code{FindIt}; only used if \code{search.lambdas=FALSE}}	
  \item{wts}{An optional set of scaling weights. The default is 1.}
  \item{scale.c}{A set of weights for recaling the pre-treatment covariates. \code{maketwoway} is useful for generating these.}
  \item{scale.int}{A set of weights for recaling the covariates to be interacted with treatment variables. \code{maketwoway} is useful for generating these.}
  \item{fit.glmnet}{Whether to fit using the coordinate descent method in glmnet (TRUE) or the regularization path method of LARS (FALSE).}
  \item{make.twoway}{
  If \code{make.twoway=TRUE}, maketwoway automatically works within \code{FindIt}.
}
  \item{make.allway}{
  If \code{make.allway=TRUE}, makeallway automatically works within \code{FindIt}.
}
  \item{reference.main}{
  If \code{make.allway=FALSE} and researchers generate a matrix of all possible interactions between factorial treatments, \code{reference} from \code{makeallway} function is better to be passed to \code{FindIt} through this argument.
}
  \item{type}{"binary" for binary outcome variables; 
              "continuous" for continuous outcome variables.}
  \item{treat.type}{"single" for interactions between a single variable and multiple covariates; "multiple" if \code{treat} is a matrix of multiple treatments.}
  \item{threshold}{
  Threshold to drop correlated columns when \code{makeallway} is used.
}
  \item{make.reference}{
  Whether to make a reference matrix to check which columns are dropped when \code{makeallway=TRUE}.
}
  \item{nway}{
  The argument passed to \code{makeallway}. In general, it is recommended to use the number of factorial treatments. The current version covers up to four way interactions.
}





}
\details{Implements the alternating line search algorithm for estimating the tuning parameters, as described in Imai and Ratkovic (2013).}
\value{
  \item{coefs}{A named vector of coefficients}
  \item{coefs.orig}{Coefficients on the original scale, if scale.c and scale.t was used}
  \item{fit}{The fitted values on an SVM scale}
  \item{names.out}{Names of the coefficients}
  \item{y}{The vector of observed outcomes}
  \item{X.c}{Matrix of pre-treatment covariates}
  \item{X.t}{Matrix of treatments or treatment-covariate interactions}
  \item{GCV}{GCV statistic at the minimum}
  \item{ATE}{For a single treatment, the estimated ATE.  For multiple treatments, the estimated ATE for each unique treatment combination}
  \item{lambdas}{Tuning parameters used for the fit}
  \item{reference}{When \code{treat.type="multiple"}, after making all interaction terms, columns with no variation or columns perfectly correlated with one of other columns are automatically dropped. \code{reference} shows which columns are kept and dropped.}
}
\references{Imai, Kosuke and Marc Ratkovic. 2013.  ``Estimating Treatment Effect Heterogeneity in Randomized Program Evaluation.'' Annals of Applied Statistics. \url{http://imai.princeton.edu/research/files/svm.pdf}}

\author{Naoki Egami, Marc Ratkovic and Kosuke Imai.}

%\seealso{\link{summary.FindIt, predict.FindIt}}

\examples{

### Example: heterogeneous effect estimation in the National Supported Work Study Data


## Load the LaLonde data
data(LaLonde)
attach(LaLonde)

y <- re78>re75 

## Generate independent varables.
white <- (black==0&hisp==0)
log.re75 <- log(1+re75)
u75 <- re75==0

## Generate pre-treatment covariates.
## The same covariates are interacted with treatment variables.
X.orig <- cbind(age,educ,black,hisp,white,marr,nodegr,log.re75,u75)
## Compute weights.
wts.extrap <- psid.wts/mean(psid.wts)
## Make the data frame.
data.single <- as.data.frame(cbind(y, treat, X.orig))


## Run to find the LASSO parameters
## F1  <-FindIt(model.treat= y ~ treat,
##              model.main= ~ age+educ+black+hisp+white+
##              marr+nodegr+log.re75+u75,
##              model.int= ~ age+educ+black+hisp+white+
##              marr+nodegr+log.re75+u75,
##              data = data.single,             
##              search.lambdas=TRUE,
##              fit.glmnet=TRUE,
##              type="binary",
##              treat.type="single") 
## With weights.
## F2  <-FindIt(model.treat= y ~ treat,
##              model.main= ~ age+educ+black+hisp+white+
##              marr+nodegr+log.re75+u75,
##              model.int= ~ age+educ+black+hisp+white+
##              marr+nodegr+log.re75+u75,
##              data = data.single, 
## 	     wts = wts.extrap,            
##              search.lambdas=TRUE,
##              fit.glmnet=TRUE,
##              type="binary",
##              treat.type="single") 

## Fit with uncovered lambda parameters.
F1  <-FindIt(model.treat= y ~ treat,
             model.main= ~ age+educ+black+hisp+white+
             marr+nodegr+log.re75+u75,
             model.int= ~ age+educ+black+hisp+white+
             marr+nodegr+log.re75+u75,
             data = data.single, 
	     lambdas =c(-3.8760,-4.0025),            
             search.lambdas=FALSE,
             fit.glmnet=TRUE,
             type="binary",
             treat.type="single") 

## With weights.
## F2  <-FindIt(model.treat= y ~ treat,
##              model.main= ~ age+educ+black+hisp+white+
##              marr+nodegr+log.re75+u75,
##              model.int= ~ age+educ+black+hisp+white+
##              marr+nodegr+log.re75+u75,
##              data = data.single, 
##              wts = wts.extrap,
##              lambdas=c(-7.376,-6.000),	            
##              search.lambdas=FALSE,
##              fit.glmnet=TRUE,
##              type="binary",
##              treat.type="single") 

## Returns a summary of coefficients and model fit.
summary(F1)
## summary(F2)

## Returns all the estimated treatment effects. 
pred1 <- predict(F1)
## Top10
pred1$data[1:10,]
## Bottom 10
pred1$data[713:722,]

## Visualize all the estimated treatment effects.
## plot(pred1)

detach(LaLonde)
	
### Example: heterogeneous effect estimation in the New Haven Get-Out-the-Vote Experiment
###
\dontrun{
data(GerberGreen)
Data1 <- GerberGreen
## Subset to single family households
Data1<-subset(Data1,persons==1)
attach(Data1)

## Make matrix of treatments. 
X.lin <- cbind(persngrp,phnscrpt,mailings,appeal)

## Generate outcome variable and matrix of pre-treatment covariates.
y <- voted98
X.orig <- cbind(age,majorpty,vote96.1,vote96.0)

## Run to search for lambdas.
## F.3<- FindIt(model.treat= y2 ~ X.lin,
##                 nway=4,
##                 model.main= ~ X.orig,
##                 make.allway=TRUE,
##                 treat.type="multiple",
##                 type="binary",
##                 search.lambdas=TRUE)

## Fit, given selected lambdas.
F3<- FindIt(model.treat= y ~ X.lin,
             nway=4,
             model.main= ~ X.orig,
             make.allway=TRUE,
             treat.type="multiple",
             type="binary",
             search.lambdas=FALSE,
 	     lambdas=c(-15.000,-6.237))

## You also can use makeallway outside FindIt
## X.t.0 <- makeallway(X.lin,nway=4)$Final
## colnames.t <- colnames(X.t.0)
## names.keep <- c("persngrp_1",paste("mailings_",1:3,sep=""), paste("phnscrpt_",1:6,sep=""))
## keeps <- unique(unlist(sapply(names.keep,F=function(x) grep(x,colnames.t))))
## X.t <- X.t.0[,keeps]
## rm(X.t.0)
## F.3<- FindIt(model.treat= y ~ X.t,
##              nway=4,
##              model.main= ~ X.orig,
##              make.allway=FALSE,
##              treat.type="multiple",
##              type="binary",
##              search.lambdas=TRUE)

## Returns coefficient estimates.
summary(F3)

## Returns predicted values for each treatment combination.
pred3 <- predict(F3,unique=TRUE)
## Top 10
pred3$data[1:10,]
## Bottom 10
pred3$data[63:72,]

## Visualize predicted values for each treatment combination.
## plot(pred3)
}
}