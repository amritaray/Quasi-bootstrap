\docType{methods}
\name{mb_statistic_fn}
\alias{mb_statistic_fn}
\title{Madsen-Browning statistic}
\usage{
  mb_statistic_fn(genotype, ped_object, Psi, p_hat, r_hat,
    map_object)
}
\arguments{
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
  This function returns the multi-locus burden statistic
}
\examples{
data(example_data)
genotype = geno_object[,2:ncol(geno_object)]
Psi = 2*kinship_fn(ped_object)
p_hat = p_hat_fn(genotype)
r_hat=r_hat_fn(genotype)
print(mb_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
}
\author{
  Ray and Gail
}
\references{
  Madsen and Browning (2009) "A Groupwise Association Test
  for Rare Mutations Using a Weighted Sum Statistic" PLoS
  Genet 5(2): e1000384
}

