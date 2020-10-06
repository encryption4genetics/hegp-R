\name{make.encrypter}
\alias{make.encrypter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create encryption keys for a Dataset
}
\description{
Sample a series of orthogonal matrices suitable for encrypting a given Dataset object.
}
\usage{
make.encrypter( D, blocksize=0, shuffle=FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{ A Dataset object, created by \code{\link{build.D}} }
  \item{blocksize}{Optional size of encryption blocks. Each block of individuals is encrypted separately. If blocksize is zero then a single encrypter is generated.}
  \item{shuffle}{Permute the random numbers as an additional security measure }
%%     ~~Describe \code{x} here~~
}
\details{
Create random orthogonal encryption keys for the dataset D created by a call to \code{\link{build.D}}. Each encryption key is a random orthogonal matrix generated from the Steifel manifold. If the dataset contains \eqn{N} individuals then if \eqn{blocksize>0}, \eqn{N/blocksize + 1} keys are generated. Most keys are of dimension \eqn{blocksize*blocksize} with the final key with smaller dimension to make the sum of the dimensions of the keys equal to \eqn{N}.
The Dataset can then be encrypted using the orthogonal keys by \code{\link{encrypt.D}}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with elements
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{blocks }{The number of blocks}
 \item{block}{ a list of encryption keys, each a matrix}
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
}
\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
