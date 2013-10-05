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
