% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConditionalEffect.R
\name{ConditionalEffect}
\alias{ConditionalEffect}
\title{Estimating the Conditional Effects with the CausalANOVA.}
\usage{
ConditionalEffect(
  object,
  treat.fac = NULL,
  cond.fac = NULL,
  base.ind = 1,
  round = 3,
  inference = NULL,
  verbose = TRUE
)
}
\arguments{
\item{object}{The output from \code{CausalANOAV} function.}

\item{treat.fac}{The name of factor acting as the main treatment variable.}

\item{cond.fac}{The name of factor acting as the conditioning (moderating)
variable.}

\item{base.ind}{An indicator for the baseline of the treatment factor.
Default is 1.}

\item{round}{Digits to round estimates. Default is 3.}

\item{inference}{(optional). This argument is mainly for internal use. It
indicates whether \code{CausalANOVA} has done inference or not.}

\item{verbose}{Whether it prints the progress.}
}
\value{
\item{CondtionalEffects}{The summary of estimated conditional
effects.} \item{...}{Arguments for the internal use.}
}
\description{
\code{ConditionalEffect} estimates a variety of conditional effects using
the ouput from \code{CausalANOVA}.
}
\details{
See Details in \code{CausalANOVA}.
}
\examples{

data(Carlson)
## Specify the order of each factor
Carlson$newRecordF<- factor(Carlson$newRecordF,ordered=TRUE,
                            levels=c("YesLC", "YesDis","YesMP",
                                     "noLC","noDis","noMP","noBusi"))
Carlson$promise <- factor(Carlson$promise,ordered=TRUE,levels=c("jobs","clinic","education"))
Carlson$coeth_voting <- factor(Carlson$coeth_voting,ordered=FALSE,levels=c("0","1"))
Carlson$relevantdegree <- factor(Carlson$relevantdegree,ordered=FALSE,levels=c("0","1"))

## ####################################### 
## Without Screening and Collapsing
## ####################################### 
#################### AMEs and two-way AMIEs ####################
fit2 <- CausalANOVA(formula=won ~ newRecordF + promise + coeth_voting + relevantdegree,
                    int2.formula = ~ newRecordF:coeth_voting,
                    data=Carlson, pair.id=Carlson$contestresp,diff=TRUE,
                    cluster=Carlson$respcodeS, nway=2)
summary(fit2)
plot(fit2, type="ConditionalEffect", fac.name=c("newRecordF","coeth_voting"))
ConditionalEffect(fit2, treat.fac="newRecordF", cond.fac="coeth_voting")
}
\references{
Egami, Naoki and Kosuke Imai. 2019. Causal Interaction in
Factorial Experiments: Application to Conjoint Analysis, Journal of the American Statistical Association.
\url{http://imai.fas.harvard.edu/research/files/int.pdf}

Lim, M. and Hastie, T. 2015. Learning interactions via hierarchical
group-lasso regularization. Journal of Computational and Graphical
Statistics 24, 3, 627--654.

Post, J. B. and Bondell, H. D. 2013. ``Factor selection and structural
identification in the interaction anova model.'' Biometrics 69, 1, 70--79.
}
\seealso{
\link{CausalANOVA}.
}
\author{
Naoki Egami and Kosuke Imai.
}
