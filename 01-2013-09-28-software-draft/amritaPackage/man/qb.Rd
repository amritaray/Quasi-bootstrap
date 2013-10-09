\docType{methods}
\name{C_fn}
\alias{C_fn}
\alias{kernel_statistic_fn}
\alias{mb_statistic_fn}
\alias{p_hat_fn}
\title{Denominator term for two default statistics: Multi-locus Burden and Linear Kernel.}
\usage{
  C_fn(map_object, p_hat, r_hat)

  kernel_statistic_fn(genotype, ped_object, Psi, p_hat,
    r_hat, map_object)

  kernel_statistic_fn(genotype, ped_object, Psi, p_hat,
    r_hat, map_object)

  mb_statistic_fn(genotype, ped_object, Psi, p_hat, r_hat,
    map_object)

  p_hat_fn(genotype, epsilon = 1e-04)
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

  \item{genotype}{is the genotype matrix with individuals
  as rows and columns with marker genotypes as number of
  minor alleles.}

  \item{ped_object}{is the user input pedigree data, where
  rows are individuals, and 6 columns as pedigree id,
  individual id, father id, mother id, gender, and
  affection status}

  \item{Psi}{is the matrix of twice kinship coefficients
  between a pair of individuals}

  \item{p_hat}{is the vector of estimated minor allele
  frequency per marker}

  \item{r_hat}{is the matrix of estimated inter-marker
  correlation coefficients}

  \item{map_object}{is the user input mapfile of the
  markers}

  \item{genotype}{is the genotype matrix with individuals
  as rows and columns with marker genotypes as number of
  minor alleles.}

  \item{ped_object}{is the user input pedigree data, where
  rows are individuals, and 6 columns as pedigree id,
  individual id, father id, mother id, gender, and
  affection status}

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

  \item{Genotype}{is the user input genotype data with rows
  as individuals and columns as markers with number of
  minor alleles.}

  \item{epsilon}{is a small quantity, if the estimate is
  less or equal to 0 the function returns epsilon; if the
  estimate is greater or equal to 1 the function returns
  1-epsilon.}
}
\description{
  This function returns the value of $c_s$ term that is
  part of the denominator for Burden and Kernel statistics.

  This function returns the multi-locus burden statistic

  Computes linear kernel statistic and degrees of freedom

  This function returns the multi-locus burden statistic

  This function returns the estimate of minor allele
  frequency for each marker.
}
\details{
  This function returns the linear Kernel statistic (Schaid
  et al.)
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
Psi = 2*kinship(ped_object)
p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(burden_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
data(example_data)
 genotype = geno_object[,2:ncol(geno_object)]
Psi = 2*kinship_fn(ped_object)
p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(kernel_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
data(example_data)
genotype = geno_object[,2:ncol(geno_object)]
Psi = 2*kinship_fn(ped_object)
p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(mb_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
data(example_data)
 genotype=geno_object[,2:ncol(geno_object)]
 p_hat_fn(genotype,epsilon)
}
\author{
  Ray and Gong

  Amrita and Gail (emails)

  Ray and Gong

  Ray and Gail
}
\references{
  Schaid (2013)- Ask Alice

  Madsen and Browning (2009) "A Groupwise Association Test
  for Rare Mutations Using a Weighted Sum Statistic" PLoS
  Genet 5(2): e1000384
}

