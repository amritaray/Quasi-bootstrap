\docType{methods}
\name{p_hat_fn}
\alias{p_hat_fn}
\title{Minor allele frequency estimate}
\usage{
  p_hat_fn(genotype, epsilon = 1e-04)
}
\arguments{
  \item{Genotype}{is the user input genotype data with rows
  as individuals and columns as markers with number of
  minor alleles.}

  \item{epsilon}{is a small quantity, if the estimate is
  less or equal to 0 the function returns epsilon; if the
  estimate is greater or equal to 1 the function returns
  1-epsilon.}
}
\description{
  This function returns the estimate of minor allele
  frequency for each marker.
}
\examples{
data(example_data)
 genotype=geno_object[,2:ncol(geno_object)]
 p_hat_fn(genotype,epsilon)
}

