\name{makellway}
\alias{makeallway}
\title{Constructing named all possible interactions from a given set of factors.}
\description{
  \code{makeallway} takes a data frame of factors, each denoting a treatment level, and returns a matrix of 0's and 1's for all possible interactions.The current version can produce up to four way interactions.}
\usage{
makeallway(X, threshold = 0.999999, deletion = TRUE, make.reference = TRUE,
sparse.use = FALSE, nway)
}
\arguments{
  \item{X}{A matrix of data frame of treatments, each column interpreted as a factor.}
  \item{threshold}{
  Threshold to drop correlated columns.
}
  \item{deletion}{
  Whether to automatically drop pairwise linearly dependent columns and columns of all 0's or 1's.
}
  \item{make.reference}{
  Whether to make a reference matrix to check which columns are dropped when \code{makeallway=TRUE}.
}
  \item{sparse.use}{
  Whether to use a sparse matrix or not.
}
  \item{nway}{The order of interaction terms. All possible interactions up to this value are computed.}
}
\details{Useful for construcing a matrix for \code{FindIt} with \code{type="multiple"}.  Takes a data frame of factors and returns a design matrix with all possible interactions.  A named matrix of 0's and 1's, for every possible interaction, is returned.  Pairwise linearly dependent columns and columns of all 0's or 1's are eliminated.  
	}
\value{
  \item{FinalData}{A named matrix of 0's and 1's.}
  \item{reference}{A matrix that keeps track of all columns and shows which columns are dropped.}
}


\author{Naoki Egami, Marc Ratkovic and Kosuke Imai.}

%\seealso{\link{summary.FindIt}}

\examples{
#See the help page for FindIt() for an example.
	
}