#'
#' @title Permutation of a mark variable
#'
#' @description
#' Permuting a mark variable. In the case of more than one mark variable,
#' only the given variable is permuted.
#'
#' @usage perm.onemkvar(ppp.marks,V.perm,...)
#'
#' @param ppp.marks A marked "\verb{ppp}" object with more than one mark.
#' @param V.perm The variable for permutation
#' @param ... further arguments passed to or from other methods.
#'
#' @return A "\verb{ppp}" object.
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#'
#' @export
#'
#' @seealso \code{\link{perm.r1r2.mkdf}}, \code{\link{perm.r1fn.mkdf}},
#' \code{\link{rpoints.r1fn.mkdf}}, \code{\link{rpoints.r1f2.1mk}}, \code{\link{simul.shift.r1f2.1mk}},
#'  \code{\link{simul.shift.r1fn.mkdf}}, \code{\link{simul.jitter.r1f2.1mk}}, \code{\link{simul.jitter.r1fn.mkdf}},
#'   \code{\link{mkcorr.d8.t6}},\code{\link{mkcorr.d8}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' perm.onemkvar(ppp.ash.infect.dbh.c3,V.perm="DBH")
#' }
#'
perm.onemkvar<-function(ppp.marks,V.perm,...){

  stopifnot(is.ppp(ppp.marks))
  stopifnot(V.perm %in% names(marks(ppp.marks)))

  mk<-marks(ppp.marks)
  mk[V.perm]<-sample(mk[,V.perm],npoints(ppp.marks))
  marks(ppp.marks)<-mk
  return(ppp.marks)
}
