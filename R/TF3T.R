#' TF3T
#'
#' Matrix representation of the prisoner's dilemma strategy which plays C unless partner played D in all of the last 3 rounds.
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
#' @usage data(TF3T)
#' @examples
#' strategies <- rbind(TF3T,GRIM,ALLD,ALLC,TFT)
"TF3T"
