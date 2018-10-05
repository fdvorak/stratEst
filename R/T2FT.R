#' T2FT
#'
#' Matrix representation of the prisoner's dilemma strategy which plays C unless partner played D in either of the last 2 rounds (2 rounds of punishment if partner plays D).
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
#' @usage data(T2FT)
#' @examples
#' strategies <- rbind(T2FT,GRIM,ALLD,ALLC,TFT)
"T2FT"
