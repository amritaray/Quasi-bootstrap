#' Kinship matrix
#'
#' This function returns a matrix of the kinship coefficients of a pair of individuals
#' 
#' @param ped_object is the user input pedigree file. This file has individuals as rows, 
#' and 6 columns: family id, individual id, father id, mother id, gender, and affection 
#' status (0 = unaffected, 1 = affected, NA = missing).
#'
#' @docType methods
#' @examples
#'  data(example_data)
#'  kinship_object=kinship_fn(ped_object)
#' head(kinship_object)
#' @export

kinship_fn = function(ped_object){
  kinship(with(data.frame(
    id=ped_object$id,
    dadid=ped_object$father_id,
    momid=ped_object$mother_id,
    sex=ped_object$sex,
    famid=ped_object$ped_id), pedigree(id,dadid,momid,sex)))
}
