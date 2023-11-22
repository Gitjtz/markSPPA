

#' @title Multitype \emph{k}th Nearest Neighbour Distance Function (i-to-any)
#'
#' @description
#' For a multitype point pattern, estimate the distribution of the distance from a point of type
#' \emph{i} to the \emph{k}th nearest other point of any type.
#'
#' @usage
#' NDdot(X, i, r = NULL, breaks = NULL, ...,
#'    k=1,correction = c("km", "rs", "han"))
#'
#' @param X A "\verb{ppp}" object. Same as in \verb{Gdot} from the package \verb{spatstat}
#' @param i The type (mark value) of the points in \verb{X }from which distances are measured. Same as in \verb{Gdot} from the package \verb{spatstat}
#' @param r Same as in \verb{Gdot} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Gdot} from the package \verb{spatstat}
#' @param ... Same as in \verb{Gdot} from the package \verb{spatstat}
#' @param k Integer, indicating the \emph{k}th nearest neighbour.
#' @param correction Same as in \verb{Gdot} from the package \verb{spatstat}
#'
#' @details
#' This function is the same as \verb{Gdot} in \verb{spatatat} when \verb{k=1}.
#'
#' @return A object of "\verb{fv}".
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked as.ppp
#' @importFrom spatstat.explore rebadge.as.dotfun
#'
#' @export
#'
#' @seealso \code{\link{NDcross}},\code{\link{NDmulti}},
#' \code{\link{ND_NDdot}},\code{\link{ND_NDdot2}},
#' \code{\link{triNDcross}}, \code{\link{NDest}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(NDdot(amacrine,"off",k=2))
#' }
#'
NDdot<-function (X, i, r = NULL, breaks = NULL, ...,
                 k=1,correction = c("km", "rs", "han")){
  X <- as.ppp(X)
  if (!is.marked(X))
    stop(paste("point pattern has no", sQuote("marks")))
  stopifnot(is.multitype(X))
  marx <- marks(X, dfok = FALSE)
  if (missing(i))
    i <- levels(marx)[1]
  I <- (marx == i)
  if (sum(I) == 0)
    stop("No points are of type i")
  J <- rep.int(TRUE, X$n)
  result <- NDmulti(X, I, J, r, breaks, k, disjoint = FALSE, ...,
                   correction = correction)
  result <- rebadge.as.dotfun(result, "ND", NULL, i)
  attr(result,"yexp")<-substitute({ND[list(i~ dot)]^{k}}(r),
                                  list(k=k,i=i,dot=~symbol("\xb7")))
  return(result)
}
