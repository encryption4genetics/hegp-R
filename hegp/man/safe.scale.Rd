\name{safe.scale}
\alias{safe.scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale genotypes in a dosage matrix by subtracting the column means and dividing by the column standard deviations
}
\description{
Scale genotypes safely taking into account the possibility a genotype may have zero variance.
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
Mott et al Genetics 2020 https://doi.org/10.1534/genetics.120.303153
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

