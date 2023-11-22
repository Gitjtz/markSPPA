#'
#' @title \verb{pcfdot} difference estimation, untidy data
#'
#' @description
#' Estimates difference between two dot-type pair correlation functions (\verb{pcfdot}).
#' Same as g_gdot.
#'
#' @usage g_gdot2(ppp.marks,V.species="Species", host="Ash",
#'     V.bistate="Infected",yes="Yes",no="No",...,
#'     r = NULL,
#'     kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
#'     correction = c("isotropic", "Ripley", "translate"),
#'     divisor = c("r", "d"),ratio = FALSE)
#'
#' @param ppp.marks A "\verb{ppp}" object with two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark(for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param r Same as in \verb{pcfdot()}from the package \verb{spatstat}
#' @param kernel Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param bw Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param stoyan Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param correction Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param divisor Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{pcfdot()} from the package \verb{spatstat}
#'
#'
#' @details
#' Bivariate difference pcf (or \verb{g} function) is estimated as (Raventós, et al., 2010).
#' \deqn{g_{l,l+m}(r)-g_{m,l+m}(r)}
#' where \eqn{g()} is the  pair-correlation function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' Use the function \code{\link{perm.r1fn.mkdf}} to permute \verb{V.bistate} mark of the
#' focal species \verb{host} for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross eval.fv
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
#' @seealso \code{\link{g_gdot}} which does the same job,
#' \code{\link{K_Kdot2}}, \code{\link{L_Ldot2}}, \code{\link{J_Jdot2}}, \code{\link{ND_NDdot2}},
#' \code{\link{trig_g}}, \code{\link{trigcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(g_gdot2(ppp.marks.c3),ylim = c(-3,0))
#' nsim=2499
#' ev.g_gdot2.ash.c3<-envelope(ppp.marks.c3,g_gdot2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.g_gdot2.ash.c3,ylim=c(-4,4))
#' }

g_gdot2<-function(ppp.marks,V.species="Species", host="Ash",
                  V.bistate="Infected",yes="Yes",no="No",...,
                  r = NULL,
                  kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
                  correction = c("isotropic", "Ripley", "translate"),
                  divisor = c("r", "d"),
                  ratio = FALSE){

  stopifnot(is.ppp(ppp.marks) && is.marked(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))
  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]
  stopifnot(host %in% marks(ppp.marks)[,V.species])

  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  pcf.al<-pcfdot(ppp.lm,yes,...,
                 r = r,
                 kernel = kernel, bw = bw, stoyan = stoyan,
                 correction = correction,
                 divisor = divisor,
                 ratio = ratio)
  pcf.am<-pcfdot(ppp.lm,no,...,
                 r = r,
                 kernel = kernel, bw = bw, stoyan = stoyan,
                 correction = correction,
                 divisor = divisor,
                 ratio = ratio)
  g_g<-eval.fv(pcf.al-pcf.am)


  attr(g_g,"yexp")<-substitute(g[list(yes, host)](r)-g[list(no, host)](r),
                               list(yes=yes,no=no,host=host))
  attr(g_g,"fname")<-c("g", paste("list(",yes,"-", no,",", host,")"))
  labl.name<-names(g_g)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(g_g,"labl")<-labls

  return(g_g)

}
