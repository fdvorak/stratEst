#' strategies.DF2011
#'
#' List of six prisoner's dilemma strategies (Dal Bo and Frechette 2011).
#'
#' \describe{
#'  The prisoner's dilemma strategies are:
#'  \item{ALLD}{Strategy which always defects.}
#'  \item{ALLC}{Strategy which always cooperates.}
#'  \item{GRIM}{Strategy which cooperates until one player defects, then GRIM defects forever.}
#'  \item{TFT}{Strategy which cooperates unless the partner defected in the last round.}
#'  \item{WSLS}{Strategy which cooperates if both players chose the same action last round, otherwise WSLS defects.Also known as PTFT.}
#'  \item{T2}{Strategy which cooperates until either player defects, then it defects twice and returns to cooperation (regardless of the actions during the punishment phase).}
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
#' @usage data(strategies.DF2011)
#' @examples
#' strategies <- strategies.DF2011[c("ALLD","ALLC","TFT","GRIM")]
#' @references
#' Dal Bo, P. and G. R. Frechette (2011): The evolution of cooperation in infinitely repeated games: Experimental evidence, \emph{American Economic Review}, 101, 411-429.
"strategies.DF2011"
