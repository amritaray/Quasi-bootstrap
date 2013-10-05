#' @export
weight_fn = function(map_object, p_hat){
  weight = if("weight" %in% names(map_object)){
    map_object$weight
  }else{
    ( p_hat * (1-p_hat) )^(-0.5)
  }
  weight
}
