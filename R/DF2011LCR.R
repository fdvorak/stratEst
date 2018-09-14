#' Data of Dal Bo and Frechette (2011)
#'
#' A dataset to set up to perform latent class regression on the repeated prisonner's dilemma experiment by Dal Bo and Frechette.
#'
#' @format A data frame with 37042 rows and 8 variables:
#' \describe{
#'   \item{supergame}{The supergame number.}
#'   \item{period}{Period of the supergame.}
#'   \item{coop}{Dummy which is one if the participant cooperated in the current round.}
#'   \item{date}{Date of the session.}
#'   \item{r}{The stage game parameter of treatment.}
#'   \item{delta}{Discount factor of the treatment.}
#'   \item{group}{Group id of two matched participants.}
#'   \item{id}{Variable which identifies a unique participant-supergame combination.}
#' }
#' @usage data(DF2011LCR)
#' @source \url{https://www.aeaweb.org/articles?id=10.1257/aer.101.1.411}
#' @references
#' Dal Bo, P. and G. R. Frechette (2011): The evolution of cooperation in innitely repeated games: Experimental evidence, \emph{American Economic Review}, 101, 411-429.
#'
"DF2011LCR"
