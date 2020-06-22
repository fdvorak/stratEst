#' Data of the rock-paper-scissors game from Wang, Xu, and Zhou (2014)
#'
#' Experimental data of 72 participants playing 300 periods of the rock-paper-scissors game in matching groups of six.
#'
#' @format A data.frame with 21.600 rows and 6 variables:
#' \describe{
#'   \item{id}{Variable which identifies the participant.}
#'   \item{game}{The identifier of the game.}
#'   \item{period}{The period within the game.}
#'   \item{choice}{A factor with three levels which indicates if the player chooses rock, paper or scissors.}
#'   \item{other_choice}{A factor with three levels which indicates if the other player chooses rock, paper or scissors.}
#'   \item{result}{A factor with three levels which indicates if the result for the player.}
#' }
#' @usage data(WXZ2014)
#' @references
#' Wang, Z., Xu, B., and Zhou, H. (2014): Social cycling and conditional responses in the Rock-Paper-Scissors game. \emph{Scientific Reports} 4, 5830.
#'
"WXZ2014"
