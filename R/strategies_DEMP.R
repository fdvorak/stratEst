#' strategies.DEMP
#'
#' List of strategies for indefinitely repeated games of strategic substitutes and complements (Dvorak, Embrey, Mengel, Peeters, 2020). Four action games with choices a, b, c, and d.
#'
#' \describe{
#'  The memory-one strategies are chracterized by a vector of four responses to the four inputs a, b, c, and d which represent the partner's action in the previous period:
#'  \item{UC}{Unconditional cooperation. Response vector aaaa.}
#'  \item{CC1}{Conditional cooperation variant 2. Response vector abcd.}
#'  \item{CC2}{Conditional cooperation variant 1. Response vector abcc.}
#'  \item{NR}{Nash reversion. Response vector accc.}
#'  \item{PNR}{Partial Nash reversion. Response vector bbcc.}
#'  \item{N}{Nash play.  Response vector cccc.}
#'  \item{P}{Punishment strategy.  Response vector dddd.}
#'  \item{M1}{Myopic strategy variant 1. Response vector dccc.}
#'  \item{M2}{Myopic strategy variant 2. Response vector bccc.}
#'
#'  Strategy variants which start with action a in period one.
#'  \item{UC.a}{Unconditional cooperation. Response vector a-aaaa.}
#'  \item{CC1.a}{Conditional cooperation variant 2. Response vector a-abcd.}
#'  \item{CC2.a}{Conditional cooperation variant 1. Response vector a-abcc.}
#'  \item{NR.a}{Nash reversion. Response vector a-accc.}
#'  \item{PNR.a}{Partial Nash reversion. Response vector a-bbcc.}
#'  \item{N.a}{Nash play.  Response vector a-cccc.}
#'  \item{P.a}{Punishment strategy.  Response vector a-dddd.}
#'  \item{M1.a}{Myopic strategy variant 1. Response vector a-dccc.}
#'  \item{M2.a}{Myopic strategy variant 2. Response vector a-bccc.}
#'
#'  Strategy variants which start with action b in period one.
#'  \item{UC.b}{Unconditional cooperation. Response vector b-aaaa.}
#'  \item{CC1.b}{Conditional cooperation variant 2. Response vector b-abcd.}
#'  \item{CC2.b}{Conditional cooperation variant 1. Response vector b-abcc.}
#'  \item{NR.b}{Nash reversion. Response vector b-accc.}
#'  \item{PNR.b}{Partial Nash reversion. Response vector b-bbcc.}
#'  \item{N.b}{Nash play.  Response vector b-cccc.}
#'  \item{P.b}{Punishment strategy.  Response vector b-dddd.}
#'  \item{M1.b}{Myopic strategy variant 1. Response vector b-dccc.}
#'  \item{M2.b}{Myopic strategy variant 2. Response vector b-bccc.}
#'
#'  Strategy variants which start with action c in period one.
#'  \item{UC.c}{Unconditional cooperation. Response vector c-aaaa.}
#'  \item{CC1.c}{Conditional cooperation variant 2. Response vector c-abcd.}
#'  \item{CC2.c}{Conditional cooperation variant 1. Response vector c-abcc.}
#'  \item{NR.c}{Nash reversion. Response vector c-accc.}
#'  \item{PNR.c}{Partial Nash reversion. Response vector c-bbcc.}
#'  \item{N.c}{Nash play.  Response vector c-cccc.}
#'  \item{P.c}{Punishment strategy.  Response vector c-dddd.}
#'  \item{M1.c}{Myopic strategy variant 1. Response vector c-dccc.}
#'  \item{M2.c}{Myopic strategy variant 2. Response vector c-bccc.}
#'
#'  Strategy variants which start with action d in period one.
#'  \item{UC.d}{Unconditional cooperation. Response vector d-aaaa.}
#'  \item{CC1.d}{Conditional cooperation variant 2. Response vector d-abcd.}
#'  \item{CC2.d}{Conditional cooperation variant 1. Response vector d-abcc.}
#'  \item{NR.d}{Nash reversion. Response vector d-accc.}
#'  \item{PNR.d}{Partial Nash reversion. Response vector d-bbcc.}
#'  \item{N.d}{Nash play.  Response vector d-cccc.}
#'  \item{P.d}{Punishment strategy.  Response vector d-dddd.}
#'  \item{M1.d}{Myopic strategy variant 1. Response vector d-dccc.}
#'  \item{M2.d}{Myopic strategy variant 2. Response vector d-bccc.}
#'
#' }
#'
#' @format Each strategy is encoded as a data.frame object. The rows of the data frame represent the states of the automaton. The first row is the start state of the automaton. Each data.frame object contains the following variables:
#' \describe{
#'   \item{\code{output.a}}{Probability for action a.}
#'   \item{\code{output.b}}{Probability for action b.}
#'   \item{\code{output.c}}{Probability for action c.}
#'   \item{\code{output.d}}{Probability for action d.}
#'   \item{\code{tremble}}{Probability of a tremble.}
#'   \item{\code{input.a}}{State transition for the input a.}
#'   \item{\code{input.b}}{State transition for the input b.}
#'   \item{\code{input.c}}{State transition for the input c.}
#'   \item{\code{input.d}}{State transition for the input d.}
#' }
#' @usage data(strategies.DEMP)
#' @examples
#' strategies <- strategies.DEMP[c("UC","CC1","CC2","NR","N")]
#' @references
#' Dvorak, F., Embrey, M., Mengel, F., Peeters, R. (2020): Eliciting strategies in indefinitely repeated games of strategic substitutes and complements, \emph{Working Paper}.
"strategies.DEMP"
