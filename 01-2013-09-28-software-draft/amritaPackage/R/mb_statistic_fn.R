#' @export
mb_statistic_fn = function(
  genotype, ped_object, Psi, p_hat, r_hat, map_object
){
  case_control = ped_object$case_control
  uuu = case_control - mean(case_control)
  weight = weight_fn(map_object, p_hat)
  CCC = C_fn(map_object, p_hat, r_hat)
  cs = sum(CCC) 
  
  xi = weight %*% t(genotype)

  rank_xi = rank(xi)

  term_one = sum(rank_xi[case_control == 1])
  
  N_cases = sum(case_control == 1)
  N_control = sum(case_control == 0)
  NNN = length(case_control)
  term_two = N_cases * (NNN + 1) / 2
  term_deno = N_cases * N_control * (NNN + 1)/12
  
  (term_one - term_two)/sqrt(term_deno)

}
