#' Effect Size Function
#'
#' Function to calculate proportion of phenotypic variance increase due to AC covariance
#'
#' @param a Additive genetic path coefficient
#' @param c Shared environment path coefficient
#' @param e Non-shared environmental path coefficient
#' @param g Path coefficient cultural transmission
#' @param b Path coefficient sibling interaction
#' @param varA Additive genetic variance
#' @param varC Shared environment variance
#' @param varE Non-shared environmental variance
#'
#' @return Proportion of phenotypic variance increase due to AC covariance
#' @examples
#' gnome_effect(a = sqrt(.4), c = sqrt(.3), e = sqrt(.3), g = 0, b = sqrt(.05))
#' gnome_effect(a = sqrt(.4), c = sqrt(.3), e = sqrt(.3), g = sqrt(.05), b = 0)
#' gnome_effect(a = sqrt(.4), c = sqrt(.3), e = sqrt(.3), g = sqrt(.05), b = sqrt(.05))
#' @export
gnome_effect <- function(a, c, e, g, b, varA = 1, varC = 1, varE = 1) {

  # Phenotypic covariances for MZ and DZ
  smz = 2 * (g + a/2 + b/2)**2*varA + (a*varA/2 + b*varA/2) * a + (a*varA/2 + b*varA/2) * b + c**2*varC + e**2
  sdz = 2 * (g + a/2 + b/2)**2*varA + a**2*varA/2 + b**2*varA/2 + c**2*varC + e**2*varE

  # Variance Increase
  varPh_noAC <- a**2*varA + c**2*varC + e**2*varE
  effect_mz <- (smz - varPh_noAC) / varPh_noAC
  effect_dz <- (sdz - varPh_noAC) / varPh_noAC

  effect = list(mz = effect_mz, dz = effect_dz)
  return(effect)
}
