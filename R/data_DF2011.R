#' Data of Dal Bo and Frechette (2011)
#'
#' The dataset of Dal Bo and Frechette (DF2011, 2011) as stratEst.data frame.
#'
#' @format A data frame with 7358 rows and 8 variables:
#' \describe{
#'   \item{treatment}{A treatment identifier of the experiment.}
#'   \item{id}{Variable which identifies a participant.}
#'   \item{game}{The supergame number.}
#'   \item{period}{The period of the supergame.}
#'   \item{choice}{A factor with two levels which is indicates if the participant cooperates (c) or defects (d) in the current period.}
#'   \item{other_choice}{A factor with two levels which indicates if the other participant cooperates (c) or defects (d) in the current period.}
#'   \item{input}{A factor with four levels which is indicates the action profile in the previous round. The first letter indicates the action of the participant, the second letter the action of the partner in the previous round. In the first round of a game the input is NA.}
#'   \item{output}{A factor with two levels which indicates if the other participant cooperates (c) or defects (d) in the current period.}
#' }
#' @usage data(data.DF2011)
#' @source \url{https://www.aeaweb.org/articles?id=10.1257/aer.101.1.411}
#' @references
#' Dal Bo, P. and G. R. Frechette (2011): The evolution of cooperation in infinitely repeated games: Experimental evidence, \emph{American Economic Review}, 101, 411-429.
#'
"data.DF2011"
