#'
#' @title  \verb{pcfcross} difference estimation, untidy data
#'
#' @description
#' Estimates difference between two corss-type pair correlation functions (\verb{pcfcross}).
#'
#' @usage
#' g_gcross2(ppp.marks,V.species="Species", host="Ash",
#'     V.bistate="Infected",yes="Yes",no="No",t.type="yn-nn",...,
#'     r = NULL, kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
#'     correction = c("isotropic", "Ripley", "translate"),
#'     divisor = c("r", "d"),ratio = FALSE)
#'
#' @param ppp.marks A ppp object with two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark(for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param t.type Character, indicating what type of test to be estimated. The value can be "yn-nn", "ny-nn", "yy-nn", and "ny-yn",see \verb{Details}
#' @param ... Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{pcfcross()}from the package \verb{spatstat}
#' @param kernel Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param bw Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param stoyan Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param correction Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param divisor Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{pcfcross()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference pcf (or \emph{g} function) is estimated as (Raventos, et al., 2010).
#'
#' if t.type=="yn-nn", \deqn{g_{m,l}(r)-g_{l,l}(r)}
#' if t.type=="ny-nn", \deqn{g_{l,m}(r)-g_{l,l}(r)}
#' if t.type=="yy-nn", \deqn{g_{m,m}(r)-g_{l,l}(r)}
#' if t.type=="ny-yn", \deqn{g_{l,l}(r)-g_{m,l}(r)}
#' where \eqn{g()} is the pair-correlation function (pcf), \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' Note that there is a symmetry in \emph{i} and \emph{j}, so that \eqn{g_{ij}(r) = g_{ji}(r)},
#' However, estimates of these functions may not be exactly symmetric in i and j,
#' because of the effect of edge corrections.
#'
#' This function tests, for example, host density-dependent infection of a forest disease.
#'
#' Performs the same job as \code{\link{g_gcross}}, but for untidy data,
#' i.e., the mark is a data frame. A mark variable (for example, \verb{Species})
#' contains a focal tree species and other tree species.
#' Another mark variable \verb{V.bistate} is only connected to the focal species.
#'
#' Use the function \code{\link{perm.r1fn.mkdf}} to permute \verb{V.bistate} mark of the
#' focal species \verb{host} for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks is.ppp is.marked
#' @importFrom spatstat.explore pcfcross eval.fv
#'
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Raventos, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' Baddeley A.,Rubak E.,Turner R. (2016). Spatial Point Patterns: Methodology and Applications with R.
#' Boca Raton, FL: CRC Press.
#'
#' @seealso \code{\link{g_gcross}} which does the same job,
#' \code{\link{L_Lcross2}}, \code{\link{J_Jcross2}},\code{\link{ND_NDcross2}},\code{\link{trigcross}},
#' \code{\link{trig_g}}.
#'
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(g_gcross2(ppp.marks.c3,t.type="ny-yn"),ylim = c(-0.2,0.5))
#' nsim=2499
#' ev.g_gcross2.ash.c3<-envelope(ppp.marks.c3,g_gcross2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.g_gcross2.ash.c3,ylim=c(-2,2))
#' }
#'
g_gcross2<-function(ppp.marks,V.species="Species", host="Ash",
                    V.bistate="Infected",yes="Yes",no="No",t.type="yn-nn",...,
                    r = NULL,
                    kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
                    correction = c("isotropic", "Ripley", "translate"),
                    divisor = c("r", "d"),
                    ratio = FALSE){

  stopifnot(is.ppp(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))
  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]
  stopifnot(host %in% marks(ppp.marks)[,V.species])

  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  switch(t.type,"yn-nn"={
    pcf.yn<-pcfcross(ppp.lm,yes,no,...,
                     r = r,
                     kernel = kernel, bw = bw, stoyan = stoyan,
                     correction = correction,
                     divisor = divisor,
                     ratio = ratio)
    pcf.nn<-pcfcross(ppp.lm,no,no,...,
                     r = r,
                     kernel = kernel, bw = bw, stoyan = stoyan,
                     correction = correction,
                     divisor = divisor,
                     ratio = ratio)
    g_g<-eval.fv(pcf.yn-pcf.nn)
    attr(g_g,"fname")<-c("g", paste("list(",yes,"-", no,",", no,")"))

  },"ny-nn"={
    pcf.ny<-pcfcross(ppp.lm,no,yes,...,
                     r = r,
                     kernel = kernel, bw = bw, stoyan = stoyan,
                     correction = correction,
                     divisor = divisor,
                     ratio = ratio)
    pcf.nn<-pcfcross(ppp.lm,no,no,...,
                     r = r,
                     kernel = kernel, bw = bw, stoyan = stoyan,
                     correction = correction,
                     divisor = divisor,
                     ratio = ratio)
    g_g<-eval.fv(pcf.ny-pcf.nn)
    attr(g_g,"fname")<-c("g", paste("list(",no,",", yes,"-", no,")"))

  }, "yy-nn"={
    pcf.yy<-pcfcross(ppp.lm,yes,yes,...,
                     r = r,
                     kernel = kernel, bw = bw, stoyan = stoyan,
                     correction = correction,
                     divisor = divisor,
                     ratio = ratio)
    pcf.nn<-pcfcross(ppp.lm,no,no,...,
                     r = r,
                     kernel = kernel, bw = bw, stoyan = stoyan,
                     correction = correction,
                     divisor = divisor,
                     ratio = ratio)
    g_g<-eval.fv(pcf.yy-pcf.nn)
    attr(g_g,"fname")<-c("g", paste("list(",yes,"-", yes,",", no,"-", no,")"))

  },"ny-yn"={
    pcf.ny<-pcfcross(ppp.lm,no,yes,...,
                     r = r,
                     kernel = kernel, bw = bw, stoyan = stoyan,
                     correction = correction,
                     divisor = divisor,
                     ratio = ratio)
    pcf.yn<-pcfcross(ppp.lm,yes,no,...,
                     r = r,
                     kernel = kernel, bw = bw, stoyan = stoyan,
                     correction = correction,
                     divisor = divisor,
                     ratio = ratio)
    g_g<-eval.fv(pcf.ny-pcf.yn)
    attr(g_g,"fname")<-c("g", paste("list(",no,"-", yes,",", yes,"-", no,")"))

  },stop("Unsupported test type!")
  )
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
