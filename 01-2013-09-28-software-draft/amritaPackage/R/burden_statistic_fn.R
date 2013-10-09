#' Multi-locus Burden statistic
#'
#' This function returns the multi-locus burden statistic
#' 
#' @param genotype is the genotype matrix
#' @param ped_object is the user input pedigree data
#' @param Psi is the matrix of twice kinship coefficients between a pair of individuals
#' @param p_hat is the vector of estimated minor allele frequency per marker
#' @param r_hat is the matrix of estimated inter-marker correlation coefficients
#' @param map_object is the user input mapfile of the markers 
#'          
#' @author Ray and Gong
#' @export
#' @examples
#'  data(example_data)
#'  genotype = geno_object[,2:ncol(geno_object)]
#' Psi = 2*kinship_fn(ped_object)
#' p_hat = p_hat_fn(genotype)
#' r_hat=r_hat_fn(genotype)
#' print(burden_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
#' @section Burden: pvalue

burden_statistic_fn = function(
  genotype, ped_object, Psi, p_hat, r_hat, map_object
){
  case_control = ped_object$case_control
  uuu = case_control - mean(case_control)
  weight = weight_fn(map_object, p_hat)
  CCC = C_fn(map_object, p_hat, r_hat)
  cs = sum(CCC) 
  
  xi = weight %*% t(genotype)
  
  chi_square = ( uuu %*% t(xi) )^2 / ( cs * t(uuu) %*% Psi %*% uuu )
  attributes(chi_square) = NULL
  #burden
  list(
    statistic = chi_square, 
    degrees_freedom = 1)
}
