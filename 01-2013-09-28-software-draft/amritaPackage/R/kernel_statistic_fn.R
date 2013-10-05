#' @export
kernel_statistic_fn = function(
  genotype, ped_object, Psi, p_hat, r_hat, map_object
){
  kernel_fn(
    genotype, ped_object, Psi, p_hat, r_hat, map_object
  )$chi_square
}
