#' Denominator term for two default statistics: Multi-locus Burden and Linear Kernel.
#' 
#' This function returns the value of $c_s$ term that is part of the denominator for 
#' Burden and Kernel statistics. 
#' 
#' @param map_object is the user input map file of all the markers where markers are rows 
#' and 4 columns: marker name, chromosome, base pair, user-specified weights. Internal 
#' weights as function of estimated minor alelle frequency will be used if user has not 
#' specified weights. 
#' @param p_hat is the estimated minor allele frequency per marker
#' @param r_hat is the estimated inter-marker correlation coefficient matrix
#' 
#' @author Ray and Gong
#' @docType methods
#' @export
#' @examples
#' data(example_data)
#'  genotype = geno_object[,2:ncol(geno_object)]
#'  p_hat = p_hat_fn(genotype)
#' r_hat=r_hat_fn(genotype)
#' print(C_fn(map_object,p_hat,r_hat))

C_fn = function(map_object, p_hat, r_hat){
  weight = weight_fn(map_object, p_hat)
  w_matrix = weight %*% t(weight) 
  sss = sqrt(p_hat * (1 - p_hat))
  s_matrix =  sss %*% t(sss)
  CCC = 2 * w_matrix * r_hat * s_matrix
  CCC
}
