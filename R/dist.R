#' Get pairwise distance
#'
#' @param codesA data.frame
#' @param codesB data.frame
#' @export
get_pair_dist <- function(codesA, codesB){
  if(!all(c(all(codesA$is_valid), all(codesB$is_valid), all(codesA$size %in% c(13, 9), all(codesB$size %in% c(13, 9)))))){
    stop("Error: either some of your codes aren't valid, or not size 9, 13")
  }
  vec_A <- get_cass_vec(codesA$code)
  vec_B <- get_cass_vec(codesB$code)
  return(retrieve_dist_pair(vec_A, vec_B))
}
