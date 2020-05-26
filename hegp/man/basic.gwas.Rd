\name{basic.gwas}
\alias{basic.gwas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{basic.gwas
%%  ~~function to do ... ~~
}
\description{
Perform a standard genome wide association analysis. Used to check that plaintext and ciphertext data produce the same gwas results.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
basic.gwas( D, mc.cores=10 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{A Dataset}
  \item{mc.cores}{Number of cores for parellelisation}
%%     ~~Describe \code{x} here~~
}
\details{
Each vector of SNP dosages is tested for assocation with the phenotype by simple linear regression.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A dataframe containing the logP of each tested SNP joined to the columns of D$map}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

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
