#'
#' @title \verb{Kdot} ratio estimation, untidy data
#'
#' @description
#' Estimates the ratio of two dot-type \emph{K} functions (\verb{Kdot} ).
#'
#' @usage
#' KoverKdot2(ppp.marks,V.species="Species", host="Ash",
#'     V.bistate="Infected",yes="Yes",no="No",...,
#'     r=NULL, breaks=NULL, correction,
#'     ratio=FALSE)
#'
#' @param ppp.marks A "\verb{ppp}" object with two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark(for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in \verb{Kdot()} from the package \verb{spatstat}
#' @param r Same as in \verb{Kdot()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Kdot()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Kdot()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{Kdot()} from the package \verb{spatstat}
#'
#' @details
#' Ratio of two \verb{Kdot} fuctions is estimated as (Wiegand, 2018).
#' \deqn{\frac{K_{l,l+m}(r)}{K_{m,l+m}(r)}}
#' where \eqn{K()} is the Ripley \emph{K} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests the strength of density on mortality.
#'
#' Use the function \code{\link{perm.r1fn.mkdf}} to permute \verb{V.bistate} mark of the
#' focal species \verb{host} for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#'
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' @seealso \code{\link{KoverKdot}} which does the same job,
#' \code{\link{govergdot2}}, \code{\link{LoverLdot2}}, \code{\link{JoverJdot2}}, \code{\link{NDoverNDdot2}}
#' \code{\link{triK_K}}, \code{\link{triKcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(KoverKdot2(ppp.marks.c3))
#' nsim=24999
#' ev.koverkdot2.ash.c3<-envelope(ppp.marks.c3,KoverKdot2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.koverkdot2.ash.c3)
#' }
#'
KoverKdot2<-function(ppp.marks,V.species="Species", host="Ash",
                     V.bistate="Infected",yes="Yes",no="No",...,
                     r=NULL, breaks=NULL, correction,
                     ratio=FALSE){

  stopifnot(is.ppp(ppp.marks) && is.marked(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))
  stopifnot(host %in% marks(ppp.marks)[,V.species] && length(levels(marks(ppp.marks)[1,V.bistate]))==2)
  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]

  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  k.al<-Kdot(ppp.lm,yes,
             r=r, breaks=breaks, correction=correction,
             ratio=ratio)
  k.am<-Kdot(ppp.lm,no,
             r=r, breaks=breaks, correction=correction,
             ratio=ratio)
  k_k<-eval.fv(k.al/k.am)

  attr(k_k,"yexp")<-substitute(K[list(yes, host)](r)/K[list(no, host)](r),
                               list(yes=yes,no=no,host=host))
  attr(k_k,"fname")<-c("K", paste("list(",yes,"/", no,",", host,")"))
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
