\name{maketwoway}
\alias{maketwoway}
\title{Constructing a named matrix of all two-way interactions.}
\description{
  \code{maketwoway} takes a set of pre-treatment covariates and returns all sensible two-way interactions.
}
\usage{
maketwoway(X, wts = 1, center = TRUE, deletion = TRUE, threshold = 0.99999,
frame.meanPre, predict = FALSE)
}
\arguments{
  \item{X}{A matrix of pre-treatment covariates.}
  \item{wts}{An optional set of extrapolation weights.}
  \item{center}{Whether the columns should be centered before interacting.}
  \item{deletion}{
  Whether to automatically drop pairwise linearly dependent columns and columns of all 0's or 1's.
}
  \item{threshold}{
  Threshold to drop correlated columns.
}
  \item{frame.meanPre}{
  The argument to use when \code{maketwoway} is called within \code{predict}.
}
  \item{predict}{
  The argument to use when \code{maketwoway} is called within \code{predict}.
}
}
\details{Useful for constructing a matrix for \code{FindIt}.  Takes a matrix and returns a matrix of two-way interactions.  A named matrix with an intercept, all two-way interactions and square terms are returned.  If \code{center=TRUE}, interactions are formed by (1) centering, (2) scaling to sd 1, and (3) interacting.  
Non-sensible interactions, such as interacting dummy variables for male and female, as well as squared terms for dummy variables, are not returned.	A  rescaling vector is returned, so that the columns of the standardized matrix can be returned to their original scale.
	}
\value{
  \item{X}{A matrix with all linear terms, squared terms, and two-way interactions from the original matrix.}
  \item{scale.X}{A vector of weights for reutrning X to its original scale.}
}


\author{Naoki Egami, Marc Ratkovic and Kosuke Imai.}

%\seealso{\link{summary.FindIt}}

\examples{
#See the help page for FindIt() for an example.
	
}