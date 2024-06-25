#' Power Results
#'
#' Power results from OpenMx for 9 settings
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{nmz}{Number of MZ twin pairs}
#'   \item{ndz}{Number of DZ twin pairs}
#'   \item{a}{Path coefficient additive genetic variance}
#'   \item{c}{Path coefficient shared environment variance}
#'   \item{e}{Path coefficient unique environment variance}
#'   \item{CT}{Input parameter cultural transmission}
#'   \item{SI}{Input parameters sibling interaction}
#'   \item{x}{Input parameter phenotypic sibling interaction}
#'   \item{PGS}{Proportion of A explained by the polygenic score}
#'   \item{A}{Proportion of A left unexplained by the polygenic score}
#'   \item{CT(m1) MZDZ}{CT for CT-only model, MZ&DZ sample}
#'   \item{SI(m2) MZDZ}{SI for SI-only model, MZ&DZ sample}
#'   \item{CT(m3) MZDZ}{CT for combined model, MZ&DZ sample}
#'   \item{SI(m3) MZDZ}{SI for combined model, MZ&DZ sample}
#'   \item{CT(m1) DZ}{CT for CT-only model, DZ-only sample}
#'   \item{SI(m2) DZ}{SI for SI-only model, DZ-only sample}
#'   \item{CT(m3) DZ}{CT for combined model, DZ-only sample}
#'   \item{SI(m3) DZ}{SI for combined model, DZ-only sample}
#'   \item{Smz}{Percentage of overall variance due to AC covariance for MZ twins}
#'   \item{Sdz}{Percentage of overall variance due to AC covariance for DZ twins}
#' }
"gnome_power_data"
