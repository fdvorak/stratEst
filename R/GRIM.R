#' GRIM
#'
#' Matrix representation of the prisoner's dilemma strategy which plays C until either player plays D, then it plays D forever.
#'
#' @format A data frame with 8 variables:
#' \describe{
#'   \item{\code{state}}{Indicates the state of the automaton.}
#'   \item{\code{r0}}{Probability to play D (output is 0).}
#'   \item{\code{r1}}{Probability to play C (output is 1).}
#'   \item{\code{tremble}}{Probability of a tremble.}
#'   \item{\code{t1}}{State transition if the action profile is CC (input is 1).}
#'   \item{\code{t2}}{State transition if the action profile is CD (input is 2).}
#'   \item{\code{t3}}{State transition if the action profile is DC (input is 3).}
#'   \item{\code{t4}}{State transition if the action profile is DD (input is 4).}
#' }
#' @usage data(GRIM)
#' @examples
#' strategies <- rbind(GRIM,ALLD,ALLC,TFT,PTFT)
"GRIM"
