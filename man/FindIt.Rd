\name{FindIt}
\alias{FindIt}
\alias{summary.FindIt}
\alias{highlow}

\title{FindIt for Estimating Heterogeneous Treatment Effects}
\description{
  \code{FindIt} takes a binary outcome and a set of pre-treatment covariates, and returns a model with the most predictive treatment-covariate interactions returned.
}
\usage{
	  FindIt(y,X.c, treat ,type="single", rescale.c=FALSE, search.lambdas=TRUE,
	  lambdas=NULL, wts=1, scale.c=1, scale.t=1, n.highlow=10)
}
\arguments{
  \item{y}{A vector of outcomes that takes values of +1 or -1.}
  \item{treat}{A vector of treatment assignments that takes values of 0 or 1.  If \code{type="multiple"}, a matrix of multiple treatments, as may arise in a factorial design. This matrix should take values 0 and 1.}
  \item{type}{"single" for interactions between a single variable and multiple covariates; "multiple" if \code{treat} is a matrix of multiple treatments.}
  \item{X.c}{A covariate matrix.}
  \item{rescale.c}{Should \code{X.c} be rescaled internally?  Defaults to FALSE; see \code{maketwoway} and the example below for implementing our function that automatically generates two-way interactions.}
  \item{search.lambdas}{Whether to search for the tuning parameters for the LASSO constraints.  If FALSE, \code{lambdas} must be supplied.}
  \item{lambdas}{Tuning parameters to be given to \code{FindIt}; only used if \code{search.lambdas=FALSE}}
  \item{wts}{An optional set of scaling weights.}
  \item{scale.c}{A set of weights for recaling the pre-treatment covariates. \code{maketwoway} is useful for generating these.}
  \item{scale.t}{A set of scaling weights for the treatment covariates.}
  \item{n.highlow}{When type is "single," returns the indices for observations with the highest and lowest treatment effect estimates.}

}
\details{Implements the alternating line search algorithm for estimating the tuning parameters, as described in Imai and Ratkovic (2012).}
\value{
  \item{coefs}{A named vector of coefficients}
  \item{coefs.orig}{Coefficients on the original scale, if scale.c and scale.t was used}
  \item{fit}{The fitted values on an SVM scale}
  \item{highlow}{Treatment effect estimates and indices for observations with the 10 highest and lowest treatment effects}
  \item{scale}{Vector of weights for scaling the coefficients}
  \item{names.out}{Names of hte coefficients}
  \item{y}{The vector of observed outcomes}
  \item{X.c}{Matrix of pre-treatment covariates}
  \item{X.t}{Matrix of treatments or treatment-covariate interactions}
  \item{GCV}{GCV statistic at the minimum}
  \item{ATE}{For a single treatment, the estimated ATE.  For multiple treatments, the estimated ATE for each unique treatment combination}
  \item{lambdas}{Tuning parameters used for the fit}
  \item{n.highlow}{Number of indices for assessing observations with highest and lowest effect estimates.}
}
\references{Imai, Kosuke and Marc Ratkovic.  ``Estimating Treatment Effect Heterogeneity in Randomized Program Evaluation.'' Forthcoming, Annals of Applied Statistics. \url{http://www.princeton.edu/~ratkovic/svm.pdf}}

\author{Marc Ratkovic and Kosuke Imai.}

%\seealso{\link{summary.FindIt}}

\examples{
\dontrun{	
###These examples reproduce the results in Imai and Ratkovic (forthcoming).
	
###
### Example: heterogeneous effect estimation in the National Supported Work Study Data
###

##Load the LaLonde data
data(LaLonde)
attach(LaLonde)

	
	y<- 2*(re78>re75)-1

	#Generate independent varables.
	white<-(black==0&hisp==0)
	log.re75<-log(1+re75)
	u75<-re75==0

	#Generate matrices of main effects and scaling factors.  Center on weighted means, then interact.
	X.orig<-cbind(age,educ,black,hisp,white,marr,nodegr,log.re75,u75)

	#Construct named matrix of all sensible two-way interactions.
	X<-maketwoway(X.orig)
	wts.extrap<-psid.wts/mean(psid.wts)
	X2<-maketwoway(X.orig,wts=wts.extrap)


	#Run to find the LASSO parameters
	#F1<-FindIt(y,X.c=X$X,treat, scale.c=X$scale.X, search.lambdas=TRUE) #Run to find the LASSO parameters
	#F2<-FindIt(y=y, X.c=X2$X, treat=treat, type="single", rescale.c=F, scale.c=X2$scale.X, search.lambdas=TRUE, wts=wts.extrap)
	
	#Fit with uncovered lambda parameters.
	F1<-FindIt(y,X.c=X$X,treat, scale.c=X$scale.X, search.lambdas=FALSE,lambdas=c(-4.000 ,-4.127) )
	F2<-FindIt(y=y, X.c=X2$X, treat=treat, type="single", rescale.c=FALSE, scale.c=X2$scale.X, search.lambdas=FALSE, lambdas=c( -5.9295, -4.6500), wts=wts.extrap)

	#Returns a summary of coefficients and model fit.
	summary(F1)
	summary(F2)
	#Return the estimated ATE.
	F1$ATE
	F2$ATE

	#Returns the ten highest and lowest treatment effects, as well as the covariate 
	#profiles on the original covariate scale (X.orig)
	highlow(F1, X0=X.orig)
	
detach(LaLonde)


###
### Example: heterogeneous effect estimation in the New Haven Get-Out-the-Vote Experiment
###

data(GerberGreen)

	#Subset to single family households
	Data1<-subset(Data1,persons==1)
	attach(Data1)

	#Make matrix of treatments.
	X.lin<-cbind(persngrp,phnscrpt,mailings,appeal)
	X.t.0<-makeallway(X.lin)
	colnames.t<-colnames(X.t.0)
	names.keep<-c("persngrp_1",paste("mailings_",1:3,sep=""), paste("phnscrpt_",1:6,sep=""))
	keeps<-unique(unlist(sapply(names.keep,F=function(x) grep(x,colnames.t))))
	X.t<-X.t.0[,keeps]
	rm(X.t.0)
	
	#Generate outcome variable and matrix of pre-treatment covariates.
	y<-(2*voted98-1)
	X.orig<-cbind(age,majorpty,vote96.1,vote96.0)
	X<-maketwoway(X.orig,center=T)

	#Run to search for lambdas.
	#F3<-FindIt(y=y, X.c=X$X, treat=X.t, type="multiple", rescale.c=FALSE , scale.c=X$s, search.lambdas=TRUE)
	
	#Fit, given selected lambdas.
	F3<-FindIt(y=y,X.c=X$X, treat=X.t, type="multiple", rescale.c=FALSE, scale.c=X$s, search.lambdas=FALSE, lambda=c(-9.9265,-6.6245))

	#Returns coefficient estimates.
	summary(F3)

	#Returns treatment effect estimates.
	F3$ATE
}

}