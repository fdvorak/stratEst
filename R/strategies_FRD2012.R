#' strategies.FRD2012
#'
#' List of eleven prisoner's dilemma strategies (Fudenberg, Rand, and Dreber 2012).
#'
#' c("ALLC","TFT","TF2T","TF3T","T2FT","T2F2T","GRIM","GRIM2","GRIM3","ALLD","DTFT")
#' \describe{
#'  The prisoner's dilemma strategies are:
#'  \item{ALLC}{Strategy which always cooperates.}
#'  \item{TFT}{Strategy which cooperates unless the partner defected in the last round.}
#'  \item{TF2T}{Strategy which cooperates unless the partner defected in the last two rounds.}
#'  \item{TF3T}{Strategy which cooperates unless the partner defected in the last three rounds.}
#'  \item{T2FT}{Strategy which cooperates unless the partner defected in either of the last two rounds.}
#'  \item{T2F2T}{Strategy which cooperates unless the partner defected for two consecutive rounds of the last three rounds.}
#'  \item{GRIM}{Strategy which cooperates until one player defects, then GRIM defects forever.}
#'  \item{GRIM2}{Strategy which cooperates until two consecutive rounds occur in which one player defected, then GRIM2 defects forever.}
#'  \item{GRIM3}{Strategy which cooperates until three consecutive rounds occur in which one player defected, then GRIM3 defects forever.}
#'  \item{ALLD}{Strategy which always defects.}
#'  \item{DTFT}{Strategy which starts with defection, then plays according to TFT.}
#' }
#'
#' @format Each strategy is encoded as a data.frame object. The rows of the data frame represent the states of the automaton. The first row is the start state of the automaton. Each data.frame object contains the following variables:
#' \describe{
#'   \item{\code{output.d}}{Probability to defect.}
#'   \item{\code{output.c}}{Probability to cooperate.}
#'   \item{\code{tremble}}{Probability of a tremble.}
#'   \item{\code{input.cc}}{State transition for the input cc.}
#'   \item{\code{input.cd}}{State transition for the input cd.}
#'   \item{\code{input.dc}}{State transition for the input dc.}
#'   \item{\code{input.dd}}{State transition for the input dd.}
#' }
#' @usage data(strategies.FRD2012)
#' @examples
#' strategies <- strategies.FRD2012[c("ALLC","ALLD","TFT","GRIM","PTFT")]
#' @references
#' Fudenberg, D., Rand, D. G., and Dreber, A. (2012): Slow to Anger and Fast to Forgive: Cooperation in an Uncertain World, \emph{American Economic Review}, 102, 720-749.
#'
"strategies.FRD2012"
