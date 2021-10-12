#' Computes Cohen's h (Cohen, 1988) for the difference in two proportions:
#'  $h = 2arcsin(\sqrt{p1}) - 2arcsin(\sqrt{p2})$ 
#' 
#' @param p1 The first proportion.
#' @param y The second proportion.
#' @return The Cohen's h value.
#' @examples
#' cohens_h(0.7, 0.75)
#' cohens_h(0.3, 0.4)

cohens_h <- function(p1, p2) {
  #  Computes Cohen's h (Cohen, 1988) for the difference in two proportions:
  #  $h = 2arcsin(\sqrt{p1}) - 2arcsin(\sqrt{p2})$ 
  
  # Args:
  #  p1: the first proportion 
  #  p2: the second proportion
  
  # Returns: 
  #  Cohen's h  
  
  #trial
  
  h <-  2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
  return(h)
}