\name{r_hat_fn}
\alias{r_hat_fn}
\title{Marker correlation}
\usage{
  r_hat_fn(genotype, epsilon = 1e-04)
}
\arguments{
  \item{Genotype}{is the user input genotype matrix.}

  \item{epsilon}{is a small quantity that is added or or
  subtracted from genotype depending on the number of minor
  alleles per marker. This adjustment #'is done so the
  genotypic variance at a marker is non-zero.}
}
\description{
  This function returns estimate of inter-marker
  correlation matrix.
}

