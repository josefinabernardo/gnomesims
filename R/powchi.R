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
#' powchi(df = 1, ncp = 10)
#' powchi(alpha = .01, df = 1, ncp = 5)
powchi = function(alpha = .05, df, ncp) {
  crit = qchisq(alpha, df, lower.tail = F)
  if (abs(ncp) < .0001) { ncp = 0 }
  power = pchisq(crit, df, ncp, lower.tail = F)
  power
}
