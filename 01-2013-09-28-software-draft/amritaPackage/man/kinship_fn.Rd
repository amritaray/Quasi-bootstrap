\docType{methods}
\name{kinship_fn}
\alias{kinship_fn}
\title{Kinship matrix}
\usage{
  kinship_fn(ped_object)
}
\arguments{
  \item{ped_object}{is the user input pedigree file. This
  file has individuals as rows, and 6 columns: family id,
  individual id, father id, mother id, gender, and
  affection status (0 = unaffected, 1 = affected, NA =
  missing).}
}
\description{
  This function returns a matrix of the kinship
  coefficients of a pair of individuals
}
\examples{
data(example_data)
 kinship_object=kinship_fn(ped_object)
head(kinship_object)
}

