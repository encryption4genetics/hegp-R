
# hegp.R
# Author: Richard Mott, UCL Genetics Institute.
# (C) 2020 Richard Mott


# Generate a set of orthogonal matrices of max dimension <=bloocksize suitable for encrypting genotype dosage and phenotype datasets of the same dimension as the one in D. The dimension of the last matrix is reduced in order to ensure the sum of the dimensions exactly equals the number of individuals in D. If blocksize<=0 then a single block to encrypt the dataset is generated

make.encrypter <- function( D, blocksize=0, shuffle=TRUE ) {
  N = nrow(D$geno)
  if ( blocksize <= 0  ) blocksize = N;

  start = 1
  end = min( blocksize, N)
  encrypter = list(N=N)
  block = list()
  b = 1
  df = NULL
  while ( end <= N ) {
    if ( end > N-100 ) end = N
    bsize = end-start+1
    block[[b]] = rustiefel( bsize, bsize, shuffle )
    df = rbind( df, c(start, end, bsize ))
    start = end+1
    if ( end == N ) {
      end = N+1
    } else {
      end = min( start + blocksize-1, N )
    }
    b = b+1
  }
  df = data.frame(df)
  names(df) = c( "start", "end", "size" )
  encrypter$blocks = df
  encrypter$block = block
  return ( encrypter )
}

# Encypt a dataset D using the encrypter

encrypt.D <- function( D, encrypter, invert=FALSE, kinship=FALSE ) {

  blocks = encrypter$blocks
  block = encrypter$block
  encrypted = list( y=rep( 0, length=length(D$y)), geno = matrix( 0, nrow=nrow(D$geno), ncol=ncol(D$geno)), cov=matrix(0, nrow=nrow(D$cov), ncol=ncol(D$cov)))
  colnames(encrypted$geno) = colnames(D$geno)
  colnames(encrypted$cov) = colnames(encrypted$cov)
  for( i in 1:nrow(blocks) ) {
    P = encrypter$block[[i]]
    if ( invert ) P = t(P)
    idx = blocks$start[i]:blocks$end[i]
    encrypted$y[idx] = P %*% D$y[idx]

    encrypted$geno[idx,] = P %*% D$geno[idx,]
    encrypted$cov[idx,] = P %*% D$cov[idx,]
    encrypted$map = D$map
  }
  names(encrypted$y) = names(D$y)
  rownames(encrypted$geno) = rownames(D$geno)
  rownames(encrypted$cov) = rownames(D$cov)
  if ( kinship ) {
    encrypted$kinship = make.kinship(encrypted$geno)
    colnames(encrypted$kinship) = rownames(encrypted$geno)
    rownames(encrypted$kinship) = rownames(encrypted$geno)
  }
  return( encrypted )
}



# create a dataset comprising a genotype dosage matrix, phenotype vector and covariate matrix

build.D <- function( y, dosages, cov=NULL, map=NULL, kinship=FALSE ) {

  y = y[!is.na(y)]
  ids = intersect( names(y) , rownames(dosages))
  y = y[match(ids, names(y), nomatch=0)]
  dosages = dosages[match(ids, rownames(dosages), nomatch=0),]
  N = length(y)

  if ( !is.numeric(y)) {
    warning( "y is not numeric")
    return(NULL)
  }

  if ( is.null(cov)) {
    cov = matrix(1, nrow=N, ncol=1)
    colnames(cov) = c("intercept")
  }

  if ( sum(is.na(y)) + sum(is.na(dosages)) + sum(is.na(cov)) > 0 ) {
    warning( "data contain missing values" )
    return(NULL)
  }

  if ( N != nrow(dosages ) | ( !is.null(cov) & N != nrow(cov) )) {
    warning( "dimensions incompatible")
    return(NULL)
  }

  y.s = scale(y)
  names(y.s) = names(y)
  af = apply(dosages, 2, mean)
  af = ifelse( af < 0.5, af, 1-af)
  geno = safe.scale(dosages)
  cov = safe.scale(cov)

  D = list( y=y.s, geno=geno, cov=cov, map=map, maf=af )
  if ( kinship ) {
    K = make.kinship(D$geno)
    D$kinship = K
  }
  return(D)
}

safe.scale <- function( mat) {
  # scale the columns of a matrix delaing safely with situation where variance of column is 0
  # if a column is not numeric then it is converted to c column of 0's
  # row and column names are preserved.

  m = apply(mat,2, function(x) {
    if ( is.numeric(x)) {
      s = sd(x)
      if ( s > 0.0) {
        x = (x-mean(x))/s
      } else {
        x = rep(0,len=length(x))
      }
    }
    else {
      x = rep(0,len=length(x))
    }
    })
  colnames(m) = colnames(mat)
  rownames(m) = rownames(mat)
  return(m)
}
# Perform a simple quantitative trait GWAS without a mixed model, regressing phenotype on genotype dosage.

basic.gwas <- function( D, mc.cores=10 ) {

  to = floor((1:mc.cores)*ncol(D$geno)/mc.cores)
  to[mc.cores] = ncol(D$geno)
  start = c(1, to[1:(mc.cores-1)]+1)
  lp = mclapply( 1:mc.cores, function( segment ) {
    #	   cat("segment", segment, start[segment],to[segment],"\n")
    r = as.numeric(cor( D$y, D$geno[,start[segment]:to[segment]]))
    r2 = r*r
    n = length(D$y)
    t = r * sqrt( (n-2)/(1-r2) )
    logP = -pt( abs(t), df=n-2, lower.tail=FALSE, log.p=TRUE)/log(10)
    return(logP) })
  logP = unlist(lp)
  return( data.frame(cbind(D$map, logP)))
}

# A simple GWAS with a mixed model...

basic.mm.gwas <- function( D, mc.cores=10) {
  y = D$y
  if( is.null(names(y))) names(y) = as.character(1:length(y))
  genotypes = D$geno
  rownames(genotypes) = names(y)
  mm = mixed.model.gwas( y, genotypes, kinship=D$kinship, nperm=0 )
  return( data.frame(cbind(D$map, logP=-log10(mm$pval[1,]))))
}

# Simulate a random orthogonal matrix of dimensions m*R using the Steifel Manifold (function adapted from R package rsteifel)
rustiefel <- function (m, R=m, shuffle=TRUE)
{
  rn = rnorm(m * R)
  if ( shuffle ) { # permute the random numbers for extra security
   f = file("/dev/random", "rb")
   data = readBin(f, integer(), size=4, n = m*R, endian="little")
   r = order(data)
   rn = rn[r]
  }

  X <- matrix(rn, m, R)
  tmp <- eigen(t(X) %*% X)
  X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow = R)) %*% t(tmp$vec))
}

# Quantile normalise a dataset

qnorm.D <- function( D, digits=NA) {
  n = length(D$y)
  q = qnorm((1:n)/(n+1))
  if ( !is.na(digits)) q= signif(q,digits=digits)
  Dq = list()
  Dq$y[order(D$y)] = q
  Dq$geno = apply( D$geno, 2, function( x, q) { z=rep(0,length(q)); z[order(x)]=q; z}, q)
  Dq$cov = apply( D$cov, 2, function( x, q){ z=rep(0,length(q)); z[order(x)]=q; z}, q)
  Dq$map = D$map
  return(Dq)
}


