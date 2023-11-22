#'
#' @title \verb{Ldot} ratio estimation, untidy data
#'
#' @description
#' Estimates the ratio of two dot-type \emph{L} functions (\verb{Ldot}).
#'
#' @usage
#' LoverLdot2(ppp.marks,V.species="Species", host="Ash",
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
#' @param correction Same as in \verb{Ldot()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{Kdot()} from the package \verb{spatstat}
#'
#' @details
#' Ratio of two \verb{Ldot} fuctions is estimated as (Wiegand, 2018).
#' \deqn{\frac{L_{l,l+m}(r)}{L_{m,l+m}(r)}}
#' where \eqn{L()} is the \emph{L} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests the strength of density on mortality.
#'
#' Use the function \code{\link{perm.r1fn.mkdf}} to permute \verb{V.bistate} mark of the
#' focal species \verb{host} for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot Ldot Lcross eval.fv
#'
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' @seealso \code{\link{LoverLdot}} which does the same job,
#' \code{\link{govergdot2}}, \code{\link{JoverJdot2}}, \\code{\link{NDoverNDdot2}}, \code{\link{KoverKdot2}},
#'  \code{\link{triL_L}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(LoverLdot2(ppp.marks.c3))
#' nsim=2499
#' ev.loverldot2.ash.c3<-envelope(ppp.marks.c3,LoverLdot2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.loverldot2.ash.c3)
#' }
#'
LoverLdot2<-function(ppp.marks,V.species="Species", host="Ash",
                     V.bistate="Infected",yes="Yes",no="No",...,
                     r=NULL, breaks=NULL, correction,
                     ratio=FALSE){

  stopifnot(is.ppp(ppp.marks) && is.marked(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))
  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]
  stopifnot(host %in% marks(ppp.marks)[,V.species] && length(levels(marks(ppp.marks)[1,V.bistate]))==2)

  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  k.al<-Ldot(ppp.lm,yes,
             r=r, breaks=breaks, correction=correction,
             ratio=ratio)
  k.am<-Ldot(ppp.lm,no,
             r=r, breaks=breaks, correction=correction,
             ratio=ratio)
  l_l<-eval.fv(k.al/k.am)

  attr(l_l,"yexp")<-substitute(L[list(yes, host)](r)/L[list(no, host)](r),
                               list(yes=yes,no=no,host=~ symbol("\xb7")))
  attr(l_l,"fname")<-c("L", paste("list(",yes,"/", no,",", "~symbol(\"\\267\")",")"))
  labl.name<-names(l_l)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(l_l,"labl")<-labls


  return(l_l)
}
