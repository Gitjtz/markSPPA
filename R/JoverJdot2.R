#'
#' @title \verb{Jdot} ratio estimation, untidy data
#'
#' @description
#' Estimates the ratio of two dot-type \emph{J} functions (\verb{Jdot}).
#'
#' @usage
#' JoverJdot2(ppp.marks,V.species="Species", host="Ash",
#'     V.bistate="Infected",yes="Yes",no="No",...,
#'     eps=NULL, r=NULL, breaks=NULL, correction=NULL)
#'
#' @param ppp.marks A "\verb{ppp}" object with two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark(for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in \verb{Jdot()} from the package \verb{spatstat}
#' @param r Same as in \verb{Jdot()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Jdot()} from the package \verb{spatstat}
#' @param eps Same as in \verb{Jdot()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Jdot()} from the package \verb{spatstat}
#'
#' @details
#' Ratio of two \verb{Jdot} fuctions is estimated as (Wiegand, 2018).
#' \deqn{\frac{J_{l,l+m}(r)}{J_{m,l+m}(r)}}
#' where \eqn{J()} is the \emph{J} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests the strength of density on mortality.
#'
#' Use the function \code{\link{perm.r1fn.mkdf}} to permute \verb{V.bistate} mark of the
#' focal species \verb{host} for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Jdot eval.fv
#'
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' @seealso \code{\link{JoverJdot}} which does the same job,
#' \code{\link{govergdot2}}, \code{\link{KoverKdot2}}, \code{\link{LoverLdot2}},
#' \code{\link{NDoverNDdot2}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(JoverJdot2(ppp.marks.c3))
#' nsim=19
#' ev.jovjdot2.ash.c3<-envelope(ppp.marks.c3,JoverJdot2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.jovjdot2.ash.c3)
#' }
#'
JoverJdot2<-function(ppp.marks,V.species="Species", host="Ash",
                     V.bistate="Infected",yes="Yes",no="No",...,
                     eps=NULL, r=NULL, breaks=NULL, correction=NULL){

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

  k.al<-Jdot(ppp.lm,yes,
             eps=eps, r=r, breaks=breaks, correction=correction)
  k.am<-Jdot(ppp.lm,no,
             eps=eps, r=r, breaks=breaks, correction=correction)
  jovj<-eval.fv(k.al/k.am)

  attr(jovj,"fname")<-c("J", paste("list(",yes,"/", no,",", "~symbol(\"\\267\")",")"))
  labl.name<-names(jovj)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(jovj,"labl")<-labls


  return(jovj)
}
