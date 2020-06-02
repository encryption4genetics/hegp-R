\name{qnorm.D}
\alias{qnorm.D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{qnorm.D
%%  ~~function to do ... ~~
}
\description{
Replace a phenotype and each vector of genotype dosages by their Normal quantiles
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qnorm.D(D,digits=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{ A Dataset object }
  \item{digits}{ Optionally truncate the digits of the quantiles, if digits>0}
%%     ~~Describe \code{x} here~~
}
\details{
The phenotype \code{D$y} and each column of the genotype matrix \code{D$genos} are replaced by a permutation of the standard Normal quantiles, to improve the security of the encryption. If digits>0 then in addition only the first few decimal digits of each quantil are kept. If digits=NA then no truncation is performed.

This function should be used after encryption by \code{\link{encrypt.D}} to add an additional level of protection to a Dataset.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
A Dataset object with quantile-normalised phenotype and genotypes. Other elements of the input Dataset are copied verbatim.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{build.D}}
\code{\link{encrypt.D}}
}
\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

