\name{safe.scale}
\alias{safe.scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale the phenotypes and genotypes in a Dataset
}
\description{
Scale the phenotypes and genotypes in a Dataset, safely taking into account the possibility a genotype may have zero variance.
}
\usage{
safe.scale(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{A numeric matrix with no missing values}
%%     ~~Describe \code{x} here~~
}
\details{
Scales the columns of mat by substracting the column mean and dividing by the column standard deviation. If the standard deviation is zero the column is set to zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A matrix with the same dimensions as mat in which each column has been scaled.
}
\references{
Mott et al Genetics 2020
}
\author{
Richard Mott
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
