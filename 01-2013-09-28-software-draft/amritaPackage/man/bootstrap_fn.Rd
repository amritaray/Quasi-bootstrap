\name{bootstrap_fn}
\alias{bootstrap_fn}
\title{Quasi-Bootstrap}
\usage{
  bootstrap_fn(N_bootstrap_reps, genotype, ped_object,
    test_statistic_fns, ...)
}
\arguments{
  \item{N_bootstrap_reps}{is the number of bootstrap
  replications}

  \item{genotype}{is the genotype matrix}

  \item{ped_object}{is the user input pedigree data}

  \item{test_statistic_fns}{is a list of test statistics.
  This includes the default list of three statistics
  (Burden, Kernel and Madsen-Browning), and any user
  specified statistic.}

  \item{map_object}{is the user input mapfile of the
  markers}
}
\description{
  Compute quasi-bootstrap pvalues for association
  statistics
}
\details{
  This function returns a list of the observed statistics,
  number of bootstrap replications and the quasi-bootstrap
  pvalues. The quasi-bootstrap method can be applied to any
  genetic data with design (case control, pedigree) to
  compute association statistics and corresponding pvalues.
  The idea is to bootstrap from the decorrelated genotype
  matrix to circumvent the problem that subjects' genotypes
  at any marker may be correlated.
}
\examples{
data(example_data)
 genotype = geno_object[,2:ncol(geno_object)]
 test_statistic_fns = list(
     burden = burden_statistic_fn,
     kernel = kernel_statistic_fn,
     mb = mb_statistic_fn)
 print(bootstrap_fn(100, genotype, ped_object, test_statistic_fns, map_object))
}
\author{
  Ray and Gong
}

