#'
#' @title \verb{Jdot} difference estimation, tidy data
#'
#' @description
#' Estimates difference between two dot-type \emph{J} functions (\verb{Jdot}).
#'
#' @usage
#' J_Jdot(ppp.marks,yes="Yes",no="No",...,
#'     eps=NULL, r=NULL, breaks=NULL, correction=NULL)
#'
#' @param ppp.marks A "\verb{ppp}" object with a dichotomous mark(for example, "Yes" and "No").
#' @param yes The value of successful category in the dichotomous mark, for example, "Yes".
#' If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category in the dichotomous mark, for example, "No".
#' If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in \verb{Jdot()} from the package \verb{spatstat}
#' @param r Same as in \verb{Jdot()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Jdot()} from the package \verb{spatstat}
#' @param eps Same as in \verb{Jdot()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Jdot()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference \emph{J} functions is estimated as (Raventós, et al., 2010).
#' \deqn{J_{l,l+m}(r)-J_{m,l+m}(r)}
#' where \eqn{J()} is \emph{J} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests, for example, the host density-dependent infection of a forest disease.
#'
#'
#' Use the function \verb{rlabel} from the package \verb{spatstat} to shuffle
#' the binary mark for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Jdot eval.fv
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' @seealso \code{\link{J_Jdot2}} which does the same job,
#' \code{\link{K_Kdot}}, \code{\link{L_Ldot}}, \code{\link{g_gdot}},\code{\link{ND_NDdot}},
#' \code{\link{triJcross}}, \code{\link{triJ_J}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(J_Jdot(ppp.ash.c3))
#' nsim=19
#' ev.j_jdot.ash.c3<-envelope(ppp.ash.c3,J_Jdot,nsim=nsim,
#'     funargs=list(yes="Yes",no="No"),simulate=rlabel)
#' plot(ev.j_jdot.ash.c3)
#' }

J_Jdot<-function(ppp.marks,yes="Yes",no="No",...,
                 eps=NULL, r=NULL, breaks=NULL, correction=NULL){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]
  if(missing(no)) no<-levels(marks(ppp.marks))[1]
  stopifnot(all(c(yes, no) %in% levels(marks(ppp.marks))))
  stopifnot(yes!=no)

  k.al<-Jdot(ppp.marks,yes,...,
             eps=eps, r=r, breaks=breaks, correction=correction)
  k.am<-Jdot(ppp.marks,no,...,
             eps=eps, r=r, breaks=breaks, correction=correction)
  j_j<-eval.fv(k.al-k.am)

  attr(j_j,"fname")<-c("J", paste("list(",yes,"-", no,",", "~symbol(\"\\267\")",")"))
  labl.name<-names(j_j)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(j_j,"labl")<-labls

  return(j_j)
}
