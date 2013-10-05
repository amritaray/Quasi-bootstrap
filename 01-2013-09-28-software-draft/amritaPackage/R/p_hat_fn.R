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
