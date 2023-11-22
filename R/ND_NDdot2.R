#'
#' @title \verb{NDdot} difference estimation, untidy data
#'
#' @description
#' Estimates difference between two dot-type \emph{k}th neareast neighbor distance functions (\verb{NDdot}).
#'
#' @usage
#' ND_NDdot2(ppp.marks,V.species="Species", host="Ash", V.bistate="Infected",
#'     yes="Yes",no="No",...,
#'     r = NULL, breaks = NULL,
#'     k=1,correction = c("km", "rs", "han"))
#'
#' @param ppp.marks A "\verb{ppp}" object with two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark(for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in \verb{NDdot()} from the package \verb{spatstat}
#' @param r Same as in \verb{NDdot()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{NDdot()} from the package \verb{spatstat}
#' @param k Integer, indicating the \emph{k}th nearest neighbor.
#' @param correction Same as in \verb{NDdot()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference \emph{ND} function is estimated as
#' \deqn{ND_{l,l+m}(r)-ND_{m,l+m}(r)}
#' where \eqn{ND()} is the kth nearest neighbor distance function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests, for example, host density-dependent infection of a forest disease.
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
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' @seealso \code{\link{ND_NDdot}} which does the same job,
#' \code{\link{g_gdot2}}, \code{\link{K_Kdot2}}, \code{\link{J_Jdot2}}, \code{\link{L_Ldot2}}
#' \code{\link{triND_ND}}, \code{\link{triNDcross}}.
#'
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(ND_NDdot2(ppp.marks.c3))
#' nsim=2499
#' ev.nd_nddot2.ash.c3<-envelope(ppp.marks.c3,ND_NDdot2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.nd_nddot2.ash.c3)
#' }
#'
ND_NDdot2<-function(ppp.marks,V.species="Species", host="Ash",
                  V.bistate="Infected",yes="Yes",no="No",...,
                  r = NULL, breaks = NULL,
                  k=1,correction = c("km", "rs", "han")){

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

  # host<-paste0(yes,"+",no)

  k.al<-NDdot(ppp.lm,yes,...,
              r=r, breaks=breaks, k=k, correction=correction)
  k.am<-NDdot(ppp.lm,no,...,
              r=r, breaks=breaks, k=k, correction=correction)
  nd_nd<-eval.fv(k.al-k.am)

  attr(nd_nd,"fname")<-c("ND", paste("list(",yes,"-", no,",", "~symbol(\"\\267\")",")"))
  labl.name<-names(nd_nd)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(nd_nd,"labl")<-labls

  return(nd_nd)
}
