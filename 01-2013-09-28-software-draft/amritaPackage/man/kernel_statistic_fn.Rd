\docType{methods}
\name{kernel_statistic_fn}
\alias{kernel_statistic_fn}
\title{Linear Kernel statistic}
\usage{
  kernel_statistic_fn(genotype, ped_object, Psi, p_hat,
    r_hat, map_object)
}
\arguments{
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
}
\description{
  Computes linear kernel statistic and degrees of freedom
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
Psi = 2*kinship_fn(ped_object)
p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(kernel_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
}
\author{
  Ray and Gong
}
\references{
  Schaid (2013)- Ask Alice
}

