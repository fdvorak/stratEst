#' Data of the rock-paper-scissors game from roshambo.me
#'
#' The dataset of 10000 observations of rock-paper-scissors as a data.frame.
#'
#' @format A data frame with 10000 rows and 6 variables:
#' \describe{
#'   \item{id}{Variable which identifies the player of a unique game.}
#'   \item{game}{The identifier of a unique game.}
#'   \item{period}{The period within the game.}
#'   \item{choice}{A factor with three levels which indicates if the player chooses rock, paper or scissors.}
#'   \item{other_choice}{A factor with three levels which indicates if the other player chooses rock, paper or scissors.}
#'   \item{result}{A factor with three levels which indicates if the result for the player.}
#' }
#' @usage data(RPS)
#' @source \url{https://roshambo.me/}
"RPS"
