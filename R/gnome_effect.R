#' Function to calculate variance increase due to AC covariance
#'
#' @param a Additive genetic variance
#' @param c Shared environment variance
#' @param e Non-shared environmental variance
#' @param g Variance due to cultural transmission
#' @param b Variance due to sibling interaction
#'
#' @return Proportion of variance increase
#' @export
#'
#' @examples
#' gnome_effect(a = sqrt(.4), c = sqrt(.3), e = sqrt(.3), g = 0, b = sqrt(.05))
#' gnome_effect(a = sqrt(.4), c = sqrt(.3), e = sqrt(.3), g = sqrt(.05), b = 0)
#' gnome_effect(a = sqrt(.4), c = sqrt(.3), e = sqrt(.3), g = sqrt(.05), b = sqrt(.05))
gnome_effect <- function(a, c, e, g, b) {

  # Phenotypic covariances for MZ and DZ
  smz = 2 * (g + a/2 + b/2)**2 + (a/2 + b/2) * a + (a/2 + b/2) * b + c**2 + e**2
  sdz = 2 * (g + a/2 + b/2)**2 + a**2/2 + b**2/2 + c**2 + e**2

  # Variance Increase
  effect_mz <- smz - 1
  effect_dz <- sdz - 1

  effect = list(mz = effect_mz, dz = effect_dz)
  return(effect)
}
