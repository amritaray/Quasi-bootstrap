#' @export
C_fn = function(map_object, p_hat, r_hat){
  weight = weight_fn(map_object, p_hat)
  w_matrix = weight %*% t(weight) 
  sss = sqrt(p_hat * (1 - p_hat))
  s_matrix =  sss %*% t(sss)
  CCC = 2 * w_matrix * r_hat * s_matrix
  CCC
}
