#' SGRIM
#'
#' Matrix represention of the semi grim strategy (Breitmoser 2015). The strategy plays C if both players played C, and D if both players played D. If one player played D and the other C, play C with some probability.
#'
#' @format A matrix with 1 row and 6 columns:
#' \describe{
#'   \item{Rows}{Each row corresponds to one states of the automaton.}{}
#'   \item{Column 1}{Enumerates the states of the automaton.}
#'   \item{Column 2}{Probability to play C given the current state of the automaton. }
#'   \item{Column 3}{State transition if the history of play in the last round was CC (input is 1).}
#'   \item{Column 4}{State transition if the history of play in the last round was CD (input is 2).}
#'   \item{Column 5}{State transition if the history of play in the last round was DC (input is 3).}
#'   \item{Column 6}{State transition if the history of play in the last round was DD (input is 4).}
#' }
#' @usage data(SGRIM)
#' @examples
#' strategies <- rbind(SGRIM,GRIM,ALLD,ALLC,TFT)
#' @references
#' Breitmoser, Y. (2015): Cooperation, but no reciprocity: Individual strategies in the repeated prisoner's dilemma, \emph{American Economic Review}, 105, 2882-2910.
"SGRIM"
