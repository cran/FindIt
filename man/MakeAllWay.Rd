\name{makellway}
\alias{makeallway}
\title{Constructing named all possible interactions from a given set of factors.}
\description{
  \code{makeallway} takes a data frame of factors, each denoting a treatment level, and returns a matrix of 0's and 1's for all possible interactions.
}
\usage{
	  makeallway(X, wts = 1)
}
\arguments{
  \item{X}{A matrix of data frame of treatments, each column interpreted as a factor.}
  \item{wts}{An optional set of extrapolation weights.}
}
\details{Useful for construcing a matrix for \code{FindIt} with \code{type="multiple"}.  Takes a data frame of factors and returns a design matrix with all possible interactions.  A named matrix of 0's and 1's, for every possible interaction, is returned.  Pairwise linearly dependent columns and columns of all 0's or 1's are eliminated.  
	}
\value{
  \item{X}{A named matrix of 0's and 1's.}
}


\author{Marc Ratkovic and Kosuke Imai.}

%\seealso{\link{summary.FindIt}}

\examples{
#See the help page for FindIt() for an example.
	
}