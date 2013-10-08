\docType{methods}
\name{C_fn}
\alias{C_fn}
\alias{kernel_fn}
\alias{kernel_statistic_fn}
\title{Denominator term for two default statistics: Multi-locus Burden and Linear Kernel.}
\usage{
  C_fn(map_object, p_hat, r_hat)

  kernel_fn(genotype, ped_object, Psi, p_hat, r_hat,
    map_object)

  kernel_statistic_fn(genotype, ped_object, Psi, p_hat,
    r_hat, map_object)
}
\arguments{
  \item{map_object}{is the user input map file of all the
  markers}

  \item{p_hat}{is the estimated minor allele frequency per
  marker}

  \item{r_hat}{is the estimated inter-marker correlation
  coefficient matrix}

  \item{genotype}{is the genotype matrix}

  \item{ped_object}{is the user input pedigree data}

  \item{Psi}{is the matrix of twice kinship coefficients
  between a pair of individuals}

  \item{p_hat}{is the vector of estimated minor allele
  frequency per marker}

  \item{r_hat}{is the matrix of estimated inter-marker
  correlation coefficients}

  \item{map_object}{is the user input mapfile of the
  markers}

  \item{genotype}{is the genotype matrix}

  \item{ped_object}{is the user input pedigree data}

  \item{Psi}{is the matrix of twice kinship coefficients
  between a pair of individuals}

  \item{p_hat}{is the vector of estimated minor allele
  frequency per marker}

  \item{r_hat}{is the matrix of estimated inter-marker
  correlation coefficients}

  \item{map_object}{is the user input mapfile of the
  markers}
}
\description{
  This function returns the value of c_s term that is part
  of the denominator for Burden and Kernel statistics.

  This function returns the linear Kernel statistic (Schaid
  et al.)

  This function returns the multi-locus burden statistic
}
\section{Kernel}{
  Linear
}
\examples{
data(example_data)
 genotype = geno_object[,2:ncol(geno_object)]
 p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(C_fn(map_object,p_hat,r_hat))
data(example_data)
 genotype = geno_object[,2:ncol(geno_object)]
Psi = 2*kinship_fn(ped_object)
p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(kernel_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
data(example_data)
genotype = geno_object[,2:ncol(geno_object)]
Psi = 2*kinship(ped_object)
p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(burden_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
}
\author{
  Ray and Gong

  Ray and Gong

  Amrita and Gail (emails)
}
\references{
  Schaid
}

