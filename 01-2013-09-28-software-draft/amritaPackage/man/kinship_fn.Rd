\name{kinship_fn}
\alias{kinship_fn}
\title{Kinship matrix}
\usage{
  kinship_fn(ped_object)
}
\arguments{
  \item{ped_object}{is the user input pedigree file.}
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

