#' Weights 
#' 
#' This function assigns weights per marker. Weight = user specified weight in the map file, 
#'  else function of estimated minor allele frequency $\hat{p}(1-\hat{p})^{-0.5}$
#' 
#' @param map_object Data frame of marker information that user inputs. 
#' This file has markers as rows and columns as name, chromosome, base pairs, 
#' user-specified weights. Internal weights as function of sample minor allele 
#' frequency will be used if user does not specify weights.
#' @param p_hat Estimate of minor allele frequency from the input genotype file.    
#' 
#' @docType methods
#' @examples
#'  data(example_data)
#' genotype = geno_object[,2:ncol(geno_object)]
#' p_hat = p_hat_fn(genotype)
#' print(weight_fn(map_object, p_hat))
#' @export

weight_fn = function(map_object, p_hat){
  weight = if("weight" %in% names(map_object)){
    map_object$weight
  }else{
    ( p_hat * (1-p_hat) )^(-0.5)
  }
  weight
}
