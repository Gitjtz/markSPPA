#'
#' @title \verb{Kdot} difference estimation, tidy data
#'
#' @description
#' Estimates difference between two dot-type \emph{K} functions (\verb{Kdot}).
#'
#' @usage
#' K_Kdot(ppp.marks,yes="Yes",no="No",...,
#'     r=NULL, breaks=NULL, correction,
#'     ratio=FALSE)
#'
#' @param ppp.marks A "\verb{ppp}" object with a dichotomous mark(for example, "Yes" and "No").
#' @param yes The value of successful category in the dichotomous mark, for example, "Yes".
#' If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category in the dichotomous mark, for example, "No".
#' If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in \verb{Kdot()} from the package \verb{spatstat}
#' @param r Same as in \verb{Kdot()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Kdot()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Kdot()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{Kdot()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference \emph{K} functions is estimated as (Raventós, et al., 2010).
#' \deqn{K_{l,l+m}(r)-K_{m,l+m}(r)}
#' where \eqn{K()} is the Ripley \emph{K} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests, for example, the host density-dependent infection of a forest disease.
#'
#'
#' Use the function \verb{rlabel()} from the package \verb{spatstat} to shuffle
#' the binary mark for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' @seealso \code{\link{K_Kdot2}} which does the same job,
#' \code{\link{g_gdot}}, \code{\link{L_Ldot}}, \code{\link{J_Jdot}}, \code{\link{ND_NDdot}},
#' \code{\link{triK_K}}, \code{\link{triKcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(K_Kdot(ppp.ash.c3))
#' nsim=19
#' ev.k_kdot.ash.c3<-envelope(ppp.ash.c3,K_Kdot,nsim=nsim,
#'     funargs=list(yes="Yes",no="No"),simulate=rlabel)
#' plot(ev.k_kdot.ash.c3)}
#'

K_Kdot<-function(ppp.marks,yes="Yes",no="No",...,
                 r=NULL, breaks=NULL, correction,
                 ratio=FALSE){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]
  if(missing(no)) no<-levels(marks(ppp.marks))[1]
  stopifnot(all(c(yes, no) %in% levels(marks(ppp.marks))))
  stopifnot(yes!=no)

  k.al<-Kdot(ppp.marks,yes,...,
             r=r, breaks=breaks, correction=correction,
             ratio=ratio)
  k.am<-Kdot(ppp.marks,no,...,
             r=r, breaks=breaks, correction=correction,
             ratio=ratio)
  k_k<-eval.fv(k.al-k.am)

  attr(k_k,"fname")<-c("K", paste("list(",yes,"-", no,",", "~symbol(\"\\267\")",")"))
  labl.name<-names(k_k)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(k_k,"labl")<-labls

  return(k_k)
}
