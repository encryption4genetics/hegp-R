\name{build.D}
\alias{build.D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
build.D
%%  ~~function to do ... ~~
}
\description{
Create a Dataset object from its constituent components
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
build.D( y, dosages, cov=NULL, map=NULL, kinship=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Numeric phenotype vector }
  \item{dosages}{ Matrix of genotype dosages}
  \item{cov}{Optional matrix of covariates}
  \item{map}{Optional data frame of information about genotypes. If supplied, the i'th row of map refers to the i'th column of the genotype dosages. }
  \item{kinship}{Optional switch to generate a geneetic relationship matrix from the genotype dosages}

%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the components
y=y.s, geno=geno, cov=cov, map=map, maf=af
\item{y}{vector of phenotypes, scaled to have zero mean and variance equal to one}
\item{geno}{matrix of genotype dosages, each column (SNP) scaled to have zero mean and variance equal to one}
\item{cov}{matrix of covariates. If the input covariate matrix is NULL this is a vector of ones}
\item{map}{optional dataframe of information about SNPs, e.g. chromosome and base-pair coordinate}
\item{af}{allele frequencies of the SNPs, computed from the genotype dosages}
\item{kinship}{optional genetic relationship matrix}
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
No missing values are allowed. The dimensions of the phenotypes and genotypes are made compatible by matching the rownames of the genotypes with the names of the phenotypes. If a genetic relationship matrix is calculated it uses the function make.kinship from the library mixed.model.
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
