#'
#' @title \verb{pcfdot} difference estimation, tidy data
#'
#' @description
#' Estimates difference between two dot-type pair correlation functions (\verb{pcfdot}).
#' Same as g_gdot2.
#'
#' @usage
#' g_gdot(ppp.marks,yes="Yes",no="No",...,
#'     r = NULL,
#'     kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
#'     correction = c("isotropic", "Ripley", "translate"),
#'     divisor = c("r", "d"),ratio = FALSE)
#'
#' @param ppp.marks A "\verb{ppp}" object with a dichotomous mark(for example, "Yes" and "No").
#' @param yes The value of successful category in the dichotomous mark, for example, "Yes".
#' If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category in the dichotomous mark, for example, "No".
#' If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param r Same as in \verb{pcfdot()}from the package \verb{spatstat}
#' @param kernel Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param bw Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param stoyan Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param correction Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param divisor Same as in \verb{pcfdot()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{pcfdot()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference pcfdot (or \emph{g} function) is estimated as (Raventós, et al., 2010).
#' \deqn{g_{l,l+m}(r)-g_{m,l+m}(r)}
#' where \eqn{g()} is the  pair-correlation function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' Use the function \verb{rlabel()} from the package \verb{spatstat} to shuffle
#' the binary mark for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity
#' @importFrom spatstat.explore pcfcross pcfdot eval.fv
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
#' @seealso \code{\link{g_gdot2}} which does the same job,
#' \code{\link{K_Kdot}}, \code{\link{L_Ldot}}, \code{\link{J_Jdot}}, \code{\link{ND_NDdot}},
#' \code{\link{trig_g}}, \code{\link{trigcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(g_gdot(ppp.ash.c3),ylim=c(-3,1))
#' nsim=2499
#' ev.g_gdot.ash.c3<-envelope(ppp.ash.c3,g_gdot,nsim=nsim,
#'     funargs=list(yes="Yes",no="No"),simulate=rlabel)
#' plot(ev.g_gdot.ash.c3,ylim=c(-4,4))
#' }
#'
g_gdot<-function(ppp.marks,yes="Yes",no="No",...,
                 r = NULL,
                 kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
                 correction = c("isotropic", "Ripley", "translate"),
                 divisor = c("r", "d"),
                 ratio = FALSE){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]
  if(missing(no)) no<-levels(marks(ppp.marks))[1]
  stopifnot(all(c(yes, no) %in% levels(marks(ppp.marks))))
  stopifnot(yes!=no)

  pcf.al<-pcfdot(ppp.marks,yes,...,
                 r = r,
                 kernel = kernel, bw = bw, stoyan = stoyan,
                 correction = correction,
                 divisor = divisor,
                 ratio = ratio)
  pcf.am<-pcfdot(ppp.marks,no,...,
                 r = r,
                 kernel = kernel, bw = bw, stoyan = stoyan,
                 correction = correction,
                 divisor = divisor,
                 ratio = ratio)
  g_g<-eval.fv(pcf.al-pcf.am)

  attr(g_g,"fname")<-c("g", paste("list(",yes,"-", no,",", "~symbol(\"\\267\")",")"))
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
