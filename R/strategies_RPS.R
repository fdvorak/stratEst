#' strategies.RPS
#'
#' List of seven rock-paper-scissors strategies.
#'
#' \describe{
#'  The rock-paper-scissors strategies are:
#'  \item{rock}{Strategy which plays rock.}
#'  \item{paper}{Strategy which plays paper.}
#'  \item{scissors}{Strategy which plays scissors.}
#'  \item{nash}{Nash equilibrium strategy which plays every action with probability one-third.}
#'  \item{mixed}{Strategy which plays every action witha certain probability.}
#'  \item{anti}{Strategy which randomizes in the first period and then plays rock after scissors, paper after rock, and scissors after paper.}
#'  \item{cyclic}{Strategy which randomizes in the first period and then plays rock after paper, paper after scissors, and scissors after rock.}
#' }
#'
#' @format Each strategy is encoded as a stratEst.strategy object. The rows of the data frame represent the states of the strategy. The first row is the start state of the strategy. Each stratEst.strategy object contains the following variables:
#' \describe{
#'   \item{\code{output.rock}}{Probability to play rock.}
#'   \item{\code{output.paper}}{Probability to play paper.}
#'   \item{\code{output.scissors}}{Probability to play scissors.}
#'   \item{\code{tremble}}{Probability of a tremble.}
#'   \item{\code{input.rock}}{State transition for the input rock.}
#'   \item{\code{input.paper}}{State transition for the input paper.}
#'   \item{\code{input.scissors}}{State transition for the input scissors.}
#' }
#' @usage data(strategies.RPS)
#' @examples
#' strategies <- strategies.RPS[c("nash","mixed","cyclic")]
"strategies.RPS"
