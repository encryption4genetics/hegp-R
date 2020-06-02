\name{.HEGP}
\alias{HEGP}

\title{Homomorphic Encryption of Genotypes and Phenotypes}
\description{Overview of the use of orthogonal matrix keeys for homomorphic encryption of quantitative genetic data}
\details{
This package implements a method of homomorphic encryption of genotypes and phenotypes suitable for quantitative genetic analysis. Full details are in the paper (Mott et al Genetics, 2020). Briefly, the encryption key is a randomly generated orthogonal matrix that is multiplied into the plaintext phenotype, genotypes and optionally any covariates to produce a ciphertext that closely resembles samples from a Normal distribution. The othogonal transformation leaves the log-likelihood of the data unchanged and hence all inferences about parameter estimates, heritability and p-values of association are unchanged by the transformation.

Specifically, suppose we have a standardised vector \eqn{y} of \eqn{n} phenotypes and an \eqn{n \times p}
matrix \eqn{G} of standardised SNP dosages.
Suppose further the phenotypic variance-covariance matrix \eqn{V} is modelled as

\deqn{V=K\sigma^2_g+I\sigma^2_e}

where \eqn{K=GG^T/p} is the additive genetic relationship matrix (sometimes called a kinship matrix). Then the mixed model
\deqn{y =  X\beta + e}
describes the relationship betwen the phenotype and a set of fixed effects, represented by the \eqn{n \times k} design matrix \eqn{X}, one of which will be a SNP (a column of \eqn{G}), others may be covariates such as sex. The residual vector \eqn{e} has the same variance matrix \eqn{V} as \eqn{y}. The variance matrix \eqn{V} can be decomposed into its matrix square root via its eigen decomposition:
\deqn{V = E^T \Lambda E = (E^T \Lambda^{0.5} E)(E^T \Lambda^{0.5} E) = A^2 }
where \eqn{E} is the orthogonal matrix of eigenvectors and \eqn{\Lambda} a diagonal matrix of positive eigenvalues. Then the mixed model transformation
\deqn{A^{-1}y =  A^{-1}X\beta + A^{-1}e}
converts the mixed model to ordinary least squares in which the variance matrix of the error is the identity matrix.

Now suppose \eqn{P} is an \eqn{n \times n} orthogonal matrix (ie so \eqn{PP^T = I}. Then if we replace the plaintext \deqn{D=[y,G,X,V]} by the ciphertext \deqn{D(P) = [Py,PG,PX,PVP^T]}  the resulting ciphertext mixed model has the same log-likelihood, and the maximum likelihood estimators of the fixed effects \eqn{\beta} are unchanged, as are the p-values of genetic association and the estimates of the variance components \eqn{\sigma^2_g, \sigma^2_e},  the heritability \eqn{h^2=\frac{\sigma^2_g}{\sigma^2_g+ \sigma^2_e}} and the linkage disequilibrium between SNPs. Thus the linear mixed model is invariant under orthogonal transformation.

Furthermore, if we sample the orthogonal transformation at random from the Stiefel Manifold then the transformed ciphertext \eqn{D(P)} closely resemble a sample of Gaussian deviates, and can be thought of as a homomorphic encryption of the plaintext \eqn{D}, that preserves all the essential characteristics of the data while obscuring the genotypes of the individuals. The Stiefel Manifold samples orthogonal matrices as follows:

(i) Simulate an \eqn{n\times n} matrix \eqn{M} whose entries are all iid N(0,1).

(ii) Compute the eigen-decomposition of the symmetric matrix \eqn{M^TM=Q^TSQ} where \eqn{Q} is \eqn{n\times n} orthogonal and \eqn{S} is diagonal with positive entries.

(iii) Return the orthogonal matrix \eqn{ P=MQ^TS^{-0.5}Q} where \eqn{S^{-0.5}} is the diagonal matrix whose elements are the reciprocals of the square roots of the eigenvalues.

The HEGP package implements this scheme. It encapsulates the phenotype, genotypes, covariates and any ancilliary files such as the genetic relationship matrix and the map of SNP coordinates as a special dataset \eqn{D}, implemented as a list of conformant R objects using the function \code{build.D()}. These Datasets can then be manipulated by the other functions in the package. For instance orthogonal encryption is acocmplished by \code{encrypt.D}. Functions to sample random orthogonal matrices are provided in \code{\link{make.encrypter}}, \code{\link{rsteifel}} (based on the Rsteifel package), and functions to perform mixed model GWAS \code{basic.mm.gwas} and simple GWAS \code{basic.gwas} are also provided. There is also a function \code{\link{qnorm.D}} to add an addition level of security by replacing each column of encrypted data by quantiles from a Gaussian distribution.
}
\author{Richard Mott}
\references{
Mott et al Genetics 2020 https://doi.org/10.1534/genetics.120.303153
}
