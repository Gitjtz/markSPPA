#'
#' @title Trivariate \verb{Kcross}
#'
#' @description
#' Estimates trivariate cross-type \emph{K} function (\verb{Kcross}).
#'
#' @usage
#' triKcross(ppp.marks,V.species="Species",host="Ash", nonhost="Cork",
#'    V.bistate="Infected",yes="Yes",...,
#'    r=NULL, breaks=NULL, correction,
#'    ratio=FALSE)
#'
#' @param ppp.marks A "\verb{ppp}" object with at least two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark (for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned. If \verb{missing}, the 1st level of \verb{V.species} will be assigned.
#' @param nonhost One of levels of \verb{V.species}, for example, a nonhost tree species. If \verb{missing}, the 2nd level of \verb{V.species} will be assigned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes A value indicates the state of the host tree of concerned, for example, "Yes" or "No". If \verb{missing}, the second level of \verb{V.bistate } will be assigned.
#' @param ... Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{Kcross()} from the package \verb{spatstat}
#'
#' @details
#' The trivariate cross type  \eqn{K_{l,h}(r)},where h represents a heterospecies,
#' was extended from bivariate cross type \eqn{K_{l,m}(r)}
#'
#' The mark is a data frame. The mark variable (for example, \verb{Species})
#' contains a focal tree species and other tree species.
#' Another mark variable \verb{V.bistate} is only connected to the focal species.
#'
#' This function tests, for example, the infected host trees aggregated/segregated with a nonhost species.
#'
#' Use the function \code{\link{rpoints.r1fn.mkdf}} to simulate a CSR pattern of the
#' antecedent pattern for envelope simulation, or use the function \code{\link{perm.r1fn.mkdf}} to
#' permute \verb{V.bistate} mark of the focal species.
#'
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
#' @seealso \code{\link{triKcross}},\code{\link{K_Kcross}}, \code{\link{K_Kcross2}},
#' \code{\link{g_gcross2}},\code{\link{K_Kdot2}},\code{\link{KoverKdot2}},
#' \code{\link{triK_K}},\code{\link{trigcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(triKcross(ppp.marks.c3))
#' nsim=2499
#' ev.triKcross.cork.c3<-envelope(ppp.marks.c3,triKcross,nsim=nsim,
#'     funargs=list(V.species = "Species",host="Ash",nonhost="Cork",V.bistate="Infected",yes="Yes"),
#'     simulate = expression(rpoints.r1fn.mkdf(ppp.marks.c3,V.species="Species",sp="Cork")))
#' plot(ev.triKcross.cork.c3)
#' }

triKcross<-function(ppp.marks,V.species="Species",host="Ash", nonhost="Cork",
                    V.bistate="Infected",yes="Yes",...,
                    r=NULL, breaks=NULL, correction,
                    ratio=FALSE){

  stopifnot(is.ppp(ppp.marks) && is.marked(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))

  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]
  if(missing(nonhost)) nonhost<-levels(marks(ppp.marks)[1,V.species])[2]
  stopifnot(all(c(host,nonhost) %in% marks(ppp.marks)[,V.species]))
  stopifnot(nlevels(marks(ppp.marks)[1,V.bistate])==2)

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]

  ppp.al<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                      marks(ppp.marks)[V.species]==host &
                      marks(ppp.marks)[V.bistate]==yes,]

  marks(ppp.al)<-factor(marks(ppp.al)[,V.species])

  k.al<-Kcross(ppp.al,host,nonhost,...,
               r=r, breaks=breaks, correction=correction,
               ratio=ratio)

  host.yes<-paste0(host,".",yes)
  attr(k.al,"yexp")<-substitute(K[list(yes, nonhost)](r),
                                list(yes=host.yes, nonhost=nonhost))
  attr(k.al,"fname")<-c("K", paste("list(",host.yes,",", nonhost,")"))

  return(k.al)

}
