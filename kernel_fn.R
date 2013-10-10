#' Multi-locus Burden statistic
#'  
#'  This function returns the multi-locus burden statistic
#' 
#' @param genotype is the genotype matrix with individuals as 
#' rows and columns with marker genotypes as number of minor alleles.
#' @param ped_object is the user input pedigree data, where rows are individuals, 
#' and 6 columns as pedigree id, individual id, father id, mother id, gender, and 
#' affection status  
#' @param Psi is the matrix of twice kinship coefficients between a pair of individuals
#' @param p_hat is the vector of estimated minor allele frequency per marker
#' @param r_hat is the matrix of estimated inter-marker correlation coefficients
#' @param map_object is the user input mapfile of the markers
#' 
#' @author  Amrita and Gail (emails)
#' @docType methods
#' @rdname qb
#' 
#' @examples
#'  data(example_data)
#' genotype = geno_object[,2:ncol(geno_object)]
#' Psi = 2*kinship(ped_object)
#' p_hat = p_hat_fn(genotype)
#' r_hat=r_hat_fn(genotype)
#' print(burden_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
#' @export

kernel_statistic_fn = function(
  genotype, ped_object, Psi, p_hat, r_hat, map_object
){
  kernel_fn(
    genotype, ped_object, Psi, p_hat, r_hat, map_object
  )$chi_square
}
