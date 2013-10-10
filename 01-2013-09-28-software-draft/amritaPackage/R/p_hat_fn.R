#' Minor allele frequency estimate 
#' 
#' This function returns the estimate of minor allele frequency for each marker.
#' 
#' @param Genotype is the user input genotype data with rows as individuals and 
#' columns as markers with number of minor alleles.
#' @param epsilon is a small quantity, if the estimate is less or equal to 0 the function 
#' returns epsilon; if the estimate is greater or equal to 1 the function 
#' returns 1-epsilon.
#' 
#' @docType methods
#' @examples
#'  data(example_data)
#'  genotype=geno_object[,2:ncol(geno_object)] 
#'  p_hat_fn(genotype,epsilon)
#' @export

p_hat_fn = function(genotype, epsilon = 0.0001){
  p_hat = colMeans(genotype)/2
  if( any(p_hat <= 0) ){
    print("p_hat <= 0")
    p_hat[p_hat <= 0] = epsilon
    
  }
  if(any(p_hat >= 1)){
    print("p_hat >= 1")
    p_hat[p_hat >= 1] =  1 - epsilon  
  }
  p_hat
}
