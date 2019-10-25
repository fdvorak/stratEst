#' ALLC
#'
#' Prisoner's dilemma strategy which always cooperates as a data frame object.
#'
#' @format A data frame with 1 row and 6 columns:
#' \describe{
#'   \item{Rows}{Each row corresponds to one state of the ALLC automaton.}{}
#'   \item{Column \code{name}}{Contains the names of the rows of the automaton.}
#'   \item{Column \code{state}}{Enumerates the states of the automaton.}
#'   \item{Column \code{r1}}{Probability to play C (output is 1) if the automaton is in this state .}
#'   \item{Column \code{trembles}}{Probability of a tremble in this state .}
#'   \item{Column \code{t1}}{State transition if the action profile in the last round is CC (input is 1).}
#'   \item{Column \code{t2}}{State transition if the action profile in the last round is CD (input is 2).}
#'   \item{Column \code{t3}}{State transition if the action profile in the last round is DC (input is 3).}
#'   \item{Column \code{t4}}{State transition if the action profile in the last round is DD (input is 4).}
#' }
#' @usage data(ALLC)
#' @examples
#' strategies <- rbind(ALLC,ALLD,TFT,GRIM,PTFT)
"ALLC"
