#' Marker correlation
#' 
#' This function returns estimate of inter-marker correlation matrix.
#' 
#' @param Genotype is the user input genotype matrix.
#' @param epsilon is a small quantity that is added or or subtracted from genotype depending on the number of minor alleles per marker. This adjustment #'is done so the genotypic variance at a marker is non-zero. 
#' 
#' @docType methods
#' @examples
#'  data(example_data)
#'  genotype=geno_object[,2:ncol(geno_object)] 
#'  print(r_hat_fn(genotype,epsilon))
#' @export

r_hat_fn = function(genotype, epsilon = 0.0001){
  NNN = nrow(genotype)
  MMM = ncol(genotype)
  genotype_fixed = do.call(cbind, lapply(1:MMM, function(mmm){
    column = genotype[, mmm]
    if(length(unique(column)) == 1){
      print("r_hat found column with only one value")
      ggg = unique(column)
      nnn = sample(seq(1, NNN), 1)
      increment = if(ggg < 0.5) epsilon else -epsilon
      column[nnn] = column[nnn] + increment
    }
    column
    }))
  cor(genotype_fixed)
  }
