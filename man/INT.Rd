\name{INT}
\alias{INT}
\title{
   Estimating the AMTE, the ATCE and the AMTIE
}
\description{
   Estimating the average marginal treatment effect (AMTE), the average treatment combination effect (ATCE) and the average marginal treatment interaction effect (AMTIE) under the assumption of independently randomized factorial treatments. Researchers need to change the baseline condition by \code{relevel} the \code{target.data}. 
}
\usage{
INT(object, target.data, column, dist = "target", base, 
    sort = TRUE, compare = FALSE, order = 2)
}
\arguments{
  \item{object}{
   An output object from \code{FindIt}.
}
  \item{target.data}{
  If \code{dist= "target"}, this data is used to compute the AMTEs, ATCEs and AMTIEs. This data should represent the target population in terms of the treatment distribution and samples. The first column contains outcome variables and the other columns should be factorial treatments of interest. 
}
  \item{column}{
  Researchers can explore the AMTIEs within any factor interactions by specifying the names of factors through this argument.  
}
  \item{dist}{
  If \code{dist="target"}, \code{INT} computes the AMTEs, ATCEs and AMTIEs for the target population specified by \code{target.data}. If \code{dist="sample"}, it computes them for the data used to fit the model. If \code{dist="unique"}, it computes them for the unique samples in the data used to fit the model. 
}
  \item{base}{
  An optional argument that specifies the baseline condition. This argument can be used only when \code{column} is used, and the order of \code{base} should follow the order of \code{column}.}
  \item{sort}{
  If \code{sort=TRUE}, the AMTE, the ATCE and the AMTIE are sorted in descending order.
}
  \item{compare}{
  If \code{compare=TRUE}, \code{INT} computes the range of all factor interactions of the order specified with \code{order}. Researchers can identify the most important factor interaction for the specified order. 
}
  \item{order}{
  \code{INT} computes the ATCEs and AMTIEs of this order. Can be 2 or 3.
}
}
\details{
  \code{INT} estimates the AMTEs, ATCEs and AMTIEs under the assumption of independently randomized factorial treatments.
}
\value{
  \item{Range of AMTIE}{The range of the estimated AMTIEs}
  \item{AMTIE}{Estimated AMTIEs for treatment combinations}
  \item{ATCE}{Estimated ATCEs and AMTIEs for treatment combinations}
  \item{Sum of AMTEs}{Sum of the estimated AMTEs of each factor and AMTIEs for treatment combinations}

}
\references{Imai, Kosuke and Marc Ratkovic. 2013.  ``Estimating Treatment Effect Heterogeneity in Randomized Program Evaluation.'' Annals of Applied Statistics, Vol.7, No.1(March), pp. 443-470. \url{http://imai.princeton.edu/research/files/svm.pdf}

Egami, Naoki and Kosuke Imai. 2015. ``Causal Interaction in High-Dimension.'' Working paper. \url{http://imai.princeton.edu/research/files/int.pdf}
}

\author{Naoki Egami, Marc Ratkovic and Kosuke Imai.}

%\seealso{\link{FindIt}}

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

