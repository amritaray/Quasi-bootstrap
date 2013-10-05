\name{C_fn}
\alias{C_fn}
\title{Denominator term for two default statistics: Multi-locus Burden and Linear Kernel.
This function returns the value of c_s term that is part of the denominator for Burden and Kernel statistics.}
\usage{
  C_fn(map_object, p_hat, r_hat)
}
\arguments{
  \item{map_object}{is the user input map file of all the
  markers}

  \item{p_hat}{is the estimated minor allele frequency per
  marker}

  \item{r_hat}{is the estimated inter-marker correlation
  coefficient matrix}
}
\description{
  Denominator term for two default statistics: Multi-locus
  Burden and Linear Kernel. This function returns the value
  of c_s term that is part of the denominator for Burden
  and Kernel statistics.
}
\examples{
x <- c(2,4,3,2,4)
 doughnut(x)
 ## Add lables
 doughnut(x, labels=LETTERS[1:5])
}

