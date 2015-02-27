\name{Immigration}
\docType{data}
\alias{Immigration}
\alias{F.conjoint}
\alias{Unifdata}
\title{Data from conjoint analysis in Hainmueller and Hopkins (2014) and Hainmueller, Hopkins and Yamamoto (2014).}
\description{This data set gives the outcomes a well as treatment assignments the conjoint analysis in Hainmueller and Hopkins (2014) and Hainmueller, Hopkins and Yamamoto (2014).}
\usage{data}

\format{A data frame consisting of 6 columns and 6980 observations (5 profiles for each 1396 respondents). 
	\tabular{llll}{
	  outcome \tab integer \tab whether a profile is chosen \tab 0,1 \cr
	  Education \tab factor \tab education  \tab 7 levels \cr
	  Gender \tab factor \tab gender  \tab \code{male} or \code{female} \cr
	  Origin \tab factor \tab origin  \tab 10 levels \cr
	  Experience \tab factor \tab job experience  \tab 4 levels \cr
	  Plans \tab factor \tab job plans  \tab 4 levels 
	}
}

\source{
Data from the conjoint analysis in Hainmueller and Hopkins (2014) and Hainmueller, Hopkins and Yamamoto (2014). Because of pairings, we randomly select one profile within each pair and define a binary outcome variable for that pair, which is equal to 1 if this profile is chosen and to 0 if the other profile is selected. Columns contain 5 factors from the original conjoint analysis. The details of levels of each factor are described in Hainmueller and Hopkins (2014) and Egami and Imai (2015).  
}

\references{
Hainmueller, J. and Hopkins, D. J. 2014. ``The hidden american immigration consensus: A conjoint analysis of attitudes toward immigrants.'' American Journal of Political Science Forthcoming.

Hainmueller, J., Hopkins, D. J., and Yamamoto, T. 2014. ``Causal inference in conjoint analysis: Understanding multidimensional choices via stated preference experiments.'' Political Analysis, Vol.22, No.1, pp. 1-30.

Imai, Kosuke and Marc Ratkovic. 2013. ``Estimating Treatment Effect Heterogeneity in Randomized Program Evaluation.'' Annals of Applied Statistics, Vol.7, No.1(March), pp. 443-470. \url{http://imai.princeton.edu/research/files/svm.pdf}

Egami, Naoki and Kosuke Imai. 2015. ``Causal Interaction in High-Dimension.'' Working paper. \url{http://imai.princeton.edu/research/files/int.pdf}
}

\examples{
################################################### 
## Conjoint Analysis: Causal Interaction.
################################################### 
data(Immigration)

\dontrun{
## The SVM classifier with a lasso constraint is estimated.
## The model includes all three-way and two-way interactions as 
## well as main effects of five factorial treatments.

## Run to search for lambdas.
F.conjoint<- FindIt(model.treat= outcome ~
                    Education+Gender+Origin+Experience+Plans,
                    nway=3,
                    data = Immigration,
                    type="binary",
                    treat.type="multiple")

## Make the full factorial design matrix as the target population. 
full <- full.FindIt(F.conjoint)

## Compute the predicted potential outcomes for the target population.
Unifdata <- predict(F.conjoint,newdata=full,sort=FALSE)$data
}
## load pre-computed F.conjoint and Unifdata
data(F.conjoint)
data(Unifdata)

## Compute AMTEs, ATCEs, and AMTIEs

## Range of each factor interaction.
## The range of the AMTEs for all factors
compare1 <- INT(F.conjoint,target.data=Unifdata,compare=TRUE,order=1)
compare1

## The range of the two-way AMTIEs for all two-way factor interactions
compare2 <- INT(F.conjoint,target.data=Unifdata,compare=TRUE,order=2)
compare2

\dontrun{
## The range of the three-way AMTIEs for all three-way factor interactions
compare3 <- INT(F.conjoint,target.data=Unifdata,compare=TRUE,order=3)
compare3
}

## Compute AMTIEs within factor interactions.
## Origin x Experience
out.OP <- INT(F.conjoint,target.data=Unifdata,
              column=c("Origin","Experience"),
              base=c("India","No.job"), order=2)
out.OP

\dontrun{
## Education x Gender x Origin
out.EGO <- INT(F.conjoint,target.data=Unifdata,
               column=c("Education","Gender","Origin"),
               base=c("No.formal","female","India"), order=3)
out.EGO
}
}


\keyword{datasets}