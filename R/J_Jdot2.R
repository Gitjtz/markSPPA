#'
#' @title \verb{Jdot} difference estimation, untidy data
#'
#' @description
#' Estimates difference between two dot-type \emph{J} functions (\verb{Jdot}).
#'
#' @usage
#' J_Jdot2(ppp.marks,V.species="Species", host="Ash", V.bistate="Infected",
#'     yes="Yes",no="No",...,
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
#' Bivariate difference \emph{J} function is estimated as (Raventós, et al., 2010).
#' \deqn{J_{l,l+m}(r)-J_{m,l+m}(r)}
#' where \eqn{J()} is the \emph{J} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests, for example, host density-dependent infection of a forest disease.
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
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' @seealso \code{\link{J_Jdot}} which does the same job,
#' \code{\link{K_Kdot2}}, \code{\link{L_Ldot2}}, \code{\link{g_gdot2}},\code{\link{ND_NDdot2}},
#' \code{\link{triJcross}}, \code{\link{triJ_J}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(J_Jdot2(ppp.marks.c3))
#' nsim=2499
#' ev.j_jdot2.ash.c3<-envelope(ppp.marks.c3,J_Jdot2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.j_jdot2.ash.c3)
#' }
#'
J_Jdot2<-function(ppp.marks,V.species="Species", host="Ash",
                  V.bistate="Infected",yes="Yes",no="No",...,
                  eps=NULL, r=NULL, breaks=NULL, correction=NULL){

  stopifnot(is.ppp(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))
  stopifnot(host %in% marks(ppp.marks)[,V.species] && length(levels(marks(ppp.marks)[1,V.bistate]))==2)
  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]

  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  k.al<-Jdot(ppp.lm,yes,...,
             eps=eps, r=r, breaks=breaks, correction=correction)
  k.am<-Jdot(ppp.lm,no,...,
             eps=eps, r=r, breaks=breaks, correction=correction)
  j_j<-eval.fv(k.al-k.am)

  attr(j_j,"yexp")<-substitute(J[list(yes, host)](r)-J[list(no, host)](r),
                               list(yes=yes,no=no,host=host))
  attr(j_j,"fname")<-c("J", paste("list(",yes,"-", no,",", host,")"))
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
