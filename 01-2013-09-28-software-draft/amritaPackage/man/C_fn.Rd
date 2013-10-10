\docType{methods}
\name{C_fn}
\alias{C_fn}
\title{Denominator term for two default statistics: Multi-locus Burden and Linear Kernel.}
\usage{
  C_fn(map_object, p_hat, r_hat)
}
\arguments{
  \item{map_object}{is the user input map file of all the
  markers where markers are rows and 4 columns: marker
  name, chromosome, base pair, user-specified weights.
  Internal weights as function of estimated minor alelle
  frequency will be used if user has not specified
  weights.}

  \item{p_hat}{is the estimated minor allele frequency per
  marker}

  \item{r_hat}{is the estimated inter-marker correlation
  coefficient matrix}
}
\description{
  This function returns the value of $c_s$ term that is
  part of the denominator for Burden and Kernel statistics.
}
\examples{
data(example_data)
 genotype = geno_object[,2:ncol(geno_object)]
 p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(C_fn(map_object,p_hat,r_hat))
}
\author{
  Ray and Gong
}

