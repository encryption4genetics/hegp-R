\name{build.D}
\alias{build.D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
build.D
%%  ~~function to do ... ~~
}
\description{
Create a Dataset object from its constituent components, namely the phenotype, genotype dosages, covariates, and optionally the physical map  and kinship matrix. Once created, such a dataset can be manipulated by other functions such as \code{\link{encrypt.D}}, \code{link{basic.gwas}}, \code{\link{basic.mm.gwas}}
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
  \item{kinship}{Optional switch to generate a genetic relationship matrix from the genotype dosages}

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
Mott et al Genetics 2020 https://doi.org/10.1534/genetics.120.303153
}
\author{
Richard Mott
}
\note{
No missing values are allowed. The dimensions of the phenotypes and genotypes are made compatible by matching the rownames of the genotypes with the names of the phenotypes. If a genetic relationship matrix is calculated it uses the function \code{make.kinship} from the library mixed.model.gwas.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{encrypt.D}}
\code{\link{basic.gwas}}
\code{\link{basic.mm.gwas}}
}
\examples{
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

