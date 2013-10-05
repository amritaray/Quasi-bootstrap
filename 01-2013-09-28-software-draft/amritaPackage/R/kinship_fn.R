#' @export
kinship_fn = function(ped_object){
  kinship(with(data.frame(
    id=ped_object$id,
    dadid=ped_object$father_id,
    momid=ped_object$mother_id,
    sex=ped_object$sex,
    famid=ped_object$ped_id), pedigree(id,dadid,momid,sex)))
}
