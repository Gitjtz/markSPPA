#' @title Multitype \emph{k}th Nearest Neighbour Distance Function (i-to-j)
#'
#' @description
#' For a multitype point pattern, estimate the distribution of the distance from a point of type
#' \emph{i} to the \emph{k}th nearest point of type \emph{j}.
#'
#' @usage
#' NDcross(X, i, j, r = NULL, breaks = NULL, ...,
#'    k=1, correction = c("rs", "km", "han"))
#'
#' @param X A "\verb{ppp}" object. Same as in \verb{Gcross} from the package \verb{spatstat}
#' @param i The type (mark value) of the points in \verb{X }from which distances are measured. Same as in \verb{Gcross} from the package \verb{spatstat}
#' @param j The type (mark value) of the points in \verb{X }from which distances are measured. Same as in \verb{Gcross} from the package \verb{spatstat}
#' @param r Same as in \verb{Gcross} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Gcross} from the package \verb{spatstat}
#' @param ... Same as in \verb{Gcross} from the package \verb{spatstat}
#' @param k Integer, indicating the \emph{k}th nearest neighbour.
#' @param correction Same as in \verb{Gcross} from the package \verb{spatstat}
#'
#' @details
#' This function is the same as \verb{Gcross} in \verb{spatatat} package when \verb{k=1}.
#'
#' Use the function \code{\link{rpoints.r1fn.mkdf}} to simulate a CSR of points corresponding to \verb{yes}
#'
#' @return A object of "\verb{fv}".
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv rebadge.as.crossfun
#'
#' @export
#' @references
#' R package \verb{spatstat}
#'
#' @seealso \code{\link{NDmulti}}, \code{\link{NDest}},  \code{\link{ND_NDcross}},
#' \code{\link{ND_NDcross2}},\code{\link{NDdot}},
#' \code{\link{triNDcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(NDcross(amacrine,"off","off",k=2))
#' }
#'
NDcross<-function (X, i, j, r = NULL, breaks = NULL, ...,
                   k=1, correction = c("rs", "km", "han")){
  X <- as.ppp(X)
  if (!is.marked(X, dfok = FALSE))
    stop(paste("point pattern has no", sQuote("marks")))
  stopifnot(is.multitype(X))
  marx <- marks(X, dfok = FALSE)
  if (missing(i))
    i <- levels(marx)[1]
  if (missing(j))
    j <- levels(marx)[2]
  I <- (marx == i)
  if (sum(I) == 0)
    stop("No points are of type i")
  if (i == j) {
    result <- NDest(X[I], r = r, breaks = breaks, k=k,...)
  }
  else {
    J <- (marx == j)
    if (sum(J) == 0)
      stop("No points are of type j")
    result <- NDmulti(X, I, J, r = r, breaks = breaks, disjoint = FALSE,
                     ...,k=k, correction = correction)
  }
  result <- rebadge.as.crossfun(result, "ND", NULL, i, j)
  attr(result,"yexp")<-substitute({ND[list(i, j)]^{k}}(r),
                               list(k=k,i=i,j=j))
  return(result)
}

