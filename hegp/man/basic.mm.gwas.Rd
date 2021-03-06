\name{basic.mm.gwas}
\alias{basic.mm.gwas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixed Model GWAS
}
\description{
Perform a mixed-model GWAS to check that plaintext and ciphertext produce the same results.  }
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

\usage{
basic.mm.gwas(D, mc.cores=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{A Dataset as generated by  \code{\link{build.D}} or \code{\link{encrypt.D}}}
  \item{mc.cores}{Number of cores over which to parallelize computation}
%%     ~~Describe \code{x} here~~
}
\details{
A standard mixed model is fitted to the data, using a SNP-based genetic reelationship matrix. The phenotype and genotype are then transformed and each transformed SNP is tested for association with the transformed phenotype. Uses the function mixed.model.gwas from the package mixed.model.gwas
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A dataframe containing the logP of each tested SNP joined to the columns of D$map
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
 \code{\link{build.D}}
 \code{\link{encrypt.D}}
 \code{\link{basic.gwas}}
}
\examples{

}
