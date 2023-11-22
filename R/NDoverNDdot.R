#'
#' @title \verb{NDdot} ratio estimation, tidy data
#'
#' @description
#' Estimates ratio between two dot-type \verb{ND} functions (\verb{NDdot}).
#'
#' @usage
#' NDoverNDdot(ppp.marks,yes="Yes",no="No",...,
#'     r = NULL, breaks = NULL,
#'     k=1,correction = c("km", "rs", "han"))
#'
#' @param ppp.marks A "\verb{ppp}" object with a dichotomous mark(for example, "Yes" and "No").
#' @param yes The value of successful category in the dichotomous mark, for example, "Yes".
#' If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category in the dichotomous mark, for example, "No".
#' If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in \verb{NDdot()} from the package \verb{spatstat}
#' @param r Same as in \verb{NDdot()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{NDdot()} from the package \verb{spatstat}
#' @param k Integer, indicating the \emph{k}th nearest neighbor.
#' @param correction Same as in \verb{NDdot()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference \emph{ND} functions is estimated as
#' \deqn{\frac{ND_{l,l+m}(r)}{ND_{m,l+m}(r)}}
#' where \eqn{ND()} is the kth nearest neighbor distance function, \emph{l} and \emph{m} are
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
#' @seealso \code{\link{NDoverNDdot2}} which does the same job,
#'  \code{\link{govergdot}}, \code{\link{LoverLdot}}, \code{\link{JoverJdot}},\code{\link{KoverKdot}},
#'  \code{\link{triNDcross}}, \code{\link{NDdot}} .
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(NDoverNDdot(ppp.ash.c3))
#' nsim=19
#' ev.ndovernddot.ash.c3<-envelope(ppp.ash.c3,NDoverNDdot,nsim=nsim,
#'     funargs=list(yes="Yes",no="No"),simulate=rlabel)
#' plot(ev.ndovernddot.ash.c3)
#' }
#'

NDoverNDdot<-function(ppp.marks,yes="Yes",no="No",...,
                   r = NULL, breaks = NULL,
                   k=1,correction = c("km", "rs", "han")){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]
  if(missing(no)) no<-levels(marks(ppp.marks))[1]
  stopifnot(all(c(yes, no) %in% levels(marks(ppp.marks))))
  stopifnot(yes!=no)

  k.al<-NDdot(ppp.marks,yes,...,
              r=r, breaks=breaks, correction=correction,
              k=k)
  k.am<-NDdot(ppp.marks,no,...,
              r=r, breaks=breaks, correction=correction,
              k=k)
  ndovernd<-eval.fv(k.al/k.am)

  attr(ndovernd,"fname")<-c("ND", paste("list(",yes,"-", no,",", "~symbol(\"\\267\")",")"))
  labl.name<-names(ndovernd)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(ndovernd,"labl")<-labls

  return(ndovernd)
}
