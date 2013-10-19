#' Quadratic Kernel statistic
#'
#' Computes linear kernel statistic and degrees of freedom
#' 
#' This function returns the linear Kernel statistic (Schaid et al.)
#' 
#' @param genotype is the genotype matrix with individuals as 
#'  rows and columns with marker genotypes as number of minor alleles.
#' @param ped_object is the user input pedigree data, where rows are individuals, 
#'  and 6 columns as pedigree id, individual id, father id, mother id, gender, and 
#'  affection status  
#' @param Psi is the matrix of twice kinship coefficients between a pair of individuals
#' @param p_hat is the vector of estimated minor allele frequency per marker
#' @param r_hat is the matrix of estimated inter-marker correlation coefficients
#' @param map_object is the user input mapfile of the markers 
#'          
#' @author Ray and Gong
#' @references Schaid (2013)- Ask Alice
#' @docType methods
#' @export
#' 
#' @examples
#'  data(example_data)
#'  genotype = geno_object[,2:ncol(geno_object)]
#' Psi = 2*kinship_fn(ped_object)
#' p_hat = p_hat_fn(genotype)
#' r_hat=r_hat_fn(genotype)
#' print(kernel_statistic_fn(genotype,ped_object, Psi, p_hat, r_hat, map_object))
#' @section Kernel: Linear

quadratic_kernel_statistic_fn = function(
  genotype, ped_object, Psi, p_hat, r_hat, map_object
){
  case_control = ped_object$case_control
  uuu = case_control - mean(case_control)
  weight = weight_fn(map_object, p_hat)
  CCC = C_fn(map_object, p_hat, r_hat)
  cs = sum(CCC)
  
  genotype = as.matrix(genotype)
  MMM = ncol(genotype)
  WWW = diag(weight, nrow = MMM, ncol = MMM)
  
  kkk = genotype %*% WWW %*% t(genotype)
  quadk = matrix(nrow=nrow(kkk),ncol=nrow(kkk),0)
   for(i in 1:nrow(kkk)){
         for(j in 1:nrow(kkk)){
               quadk[i,j]=(1+kkk[i,j])^2
               }
         }
  
  quad_kernel = t(uuu)%*%quadk%*%uuu
  
  attributes(quad_kernel) = NULL
  
  list(
       statistic = quad_kernel 
       )
}
