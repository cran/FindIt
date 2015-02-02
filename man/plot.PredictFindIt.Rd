\name{plot.PredictFindIt}
\alias{plot.PredictFindIt}
\title{
Plot estimated treatment effects or predicted outcomes for each treatment combination.
}
\description{
Plot estimated treatment effects when \code{treat.type="single"} and predicted outcomes for each treatment combination when \code{treat.type="multiple"}.
}
\usage{ 
\method{plot}{PredictFindIt}(x, labels = "index", ...)
}
\arguments{
  \item{x}{
  Output from \code{predict.FindIt}.
}
  \item{labels}{
  Arguments passed to plot.
}
  \item{\dots}{
  Further arguments passed to or from other methods.
}
}
\details{
	Plot estimated treatment effects when \code{treat.type="single"} and predicted outcomes for each treatment combination when \code{treat.type="multiple"}.
}
\value{
   \item{plot}{Plot estimated treatment effects when \code{treat.type="single"} and predicted outcomes for each treatment combination when \code{treat.type="multiple"}.}
}

\author{Naoki Egami, Marc Ratkovic and Kosuke Imai.}

%\seealso{\link{FindIt}}

\examples{
## See the help page for FindIt() for an example.
	
}
