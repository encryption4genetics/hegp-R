\name{rustiefel}
\alias{rustiefel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rustiefel
}
\description{
Simulate a random orthogonal matrix of dimensions \eqn{m*R} using the Steiefel manifold
}
\usage{
rustiefel(m, R=m, shuffle=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{ the number of rows of the simulated matrix}
  \item{R}{ the number of columns }
  \item{shuffle}{ optionally shuffle the data using the computer's entropy for added security}
%%     ~~Describe \code{x} here~~
}
\details{
Function adapted from R package rsteifel.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ..
Returns an orthogonal matrix of dimension \eqn{m*R}
}
\references{
Mott et al Genetics 2020 https://doi.org/10.1534/genetics.120.303153
R package rstiefel https://cran.r-project.org/web/packages/rstiefel/index.html
}
\author{
Richard Mott
}
\note{
Function adapted from R package rsteifel https://cran.r-project.org/web/packages/rstiefel/index.html.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
