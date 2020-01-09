#' M1BF
#'
#' Matrix representation of the prisoner's dilemma strategy which plays if both players played C, and D if both players played D. If the own action was C and the other player played D, play C with some probability. If the own action was D and the other player played C, play C with some (potentially different) probability.
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
#' @usage data(M1BF)
#' @examples
#' strategies <- rbind(M1BF,GRIM,ALLD,ALLC,TFT)
"M1BF"
