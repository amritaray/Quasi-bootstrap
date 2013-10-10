#' Quasi-Bootstrap
#'
#' Compute quasi-bootstrap pvalues for association statistics
#' 
#' This function returns a list of the observed statistics, 
#' number of bootstrap replications and the quasi-bootstrap pvalues. The quasi-bootstrap 
#' method can be applied to any genetic data with design (case control, pedigree) to compute 
#' association statistics and corresponding pvalues. The idea is to bootstrap from the 
#' decorrelated genotype matrix to circumvent the problem that subjects' genotypes at any 
#' marker may be correlated.
#' 
#' @param N_bootstrap_reps is the number of bootstrap replications
#' @param genotype is the genotype matrix
#' @param ped_object is the user input pedigree data
#' @param test_statistic_fns is a list of test statistics. This includes the 
#'          default list of three statistics (Burden, Kernel and Madsen-Browning),  
#'          and any user specified statistic.
#' @param map_object is the user input mapfile of the markers 
#'          
#' @author Ray and Gong
#' @export
#' @docType methods
#' @examples
#'  data(example_data)
#'  genotype = geno_object[,2:ncol(geno_object)]
#'  test_statistic_fns = list(
#'      burden = burden_statistic_fn,
#'      kernel = kernel_statistic_fn,
#'      mb = mb_statistic_fn)
#'  print(bootstrap_fn(100, genotype, ped_object, test_statistic_fns, map_object))

bootstrap_fn = function(
  N_bootstrap_reps, genotype, ped_object, test_statistic_fns, ...
){
  Psi = 2 * kinship_fn(ped_object)                                        
  keep_q = !is.na(ped_object$case_control)
  genotype = genotype[keep_q, ]
  ped_object = ped_object[keep_q, ]  
  Psi = Psi[keep_q, keep_q]
  
  p_hat = p_hat_fn(genotype) 
  r_hat = r_hat_fn(genotype)
  observed = data.frame(t(sapply(names(test_statistic_fns), function(name){
    test_statistic_fn = test_statistic_fns[[name]]
    test_statistic_fn(genotype, ped_object, Psi, p_hat, r_hat, ...)$statistic
  })))
  
  LLL = t(chol(Psi))
  L_inverse = solve(LLL)
  genotype_centered = genotype - rep(1, nrow(genotype)) %*% t(2 * p_hat)  
  AAA = L_inverse %*% as.matrix(genotype_centered)                        
  
  bootstrap_reps = do.call(rbind, lapply(
    1:N_bootstrap_reps, 
    function(dummy){ 
      NNN = nrow(AAA)
      indices = sample(1:NNN, size=NNN, replace=TRUE)
      A_b = AAA[indices,]
      genotype_b = LLL %*% A_b
      data.frame(t(sapply(names(test_statistic_fns), function(name){
        test_statistic_fn = test_statistic_fns[[name]]
        test_statistic_fn(genotype_b, ped_object, Psi, p_hat, r_hat, ...)$statistic
      })))   
   }))
  p_value = data.frame(t(mapply(function(ooo, bbb){ 
    mean(ooo <= bbb)
  }, observed, bootstrap_reps)))
  list(observed = observed,
       bootstrap_reps = bootstrap_reps,
       p_value = p_value)                     
}
