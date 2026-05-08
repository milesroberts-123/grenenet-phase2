#' Get upper triangle only of a matrix
#'
#' @param mat a matrix of numbers
#'
#' @return matrix
#' @export
#'
#' @examples
get_upper_tri <- function(mat){
  mat[lower.tri(mat)]<- NA
  diag(mat) <- NA
  return(mat)
}