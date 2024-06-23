#' Power Function
#'
#' Function to calculate power from the non-centrality parameter
#'
#' @importFrom stats pchisq qchisq
#' @param alpha Type II error
#' @param df Degrees of freedom
#' @param ncp Non-centrality parameter
#'
#' @return Power
#' @export
#'
#' @examples
#' gnome_power(df = 1, ncp = 10)
#' gnome_power(alpha = .01, df = 1, ncp = 5)

gnome_power <- function(alpha = .05, df, ncp) {
  critical_chi2 <- qchisq(alpha, df, lower.tail = F)
  if (abs(ncp) < .0001) { ncp = 0 }
  power <- pchisq(critical_chi2, df, ncp, lower.tail = F)
  power
}
