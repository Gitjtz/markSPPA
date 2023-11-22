#'
#' @title Conspecific trivariate random labeling, untidy data
#'
#' @description
#' Trivariate random labeling function for one variate with a qualitative mark.
#'
#' @usage
#' contrl2(ppp.marks, V.species="Species",host="Ash",
#'     V.bistate="Infected",yes="Yes",...,
#'     r = NULL, kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
#'     correction = c("isotropic", "Ripley", "translate"),
#'     divisor = c("r", "d"),ratio = FALSE)
#'
#' @param ppp.marks A "\verb{ppp}" object with more than one marks. One is a multiple category mark (for example, species in a forest) and another is a dichotomous mark(for example, "Yes" and "No" to indicates whether the host species is infected).
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of \verb{V.species}, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes A value indicates the state of the host tree of concerned, for example, "Yes" or "No".
#' If \verb{missing}, the second level of \verb{V.bistate } will be assigned.
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
#' For example, this function tests the influence of the spatial pattern of
#' conspecific host trees on infection (Raventos et al. 2010).
#' \deqn{p_{l,l+m}(r)=\frac{p_lg_{l,l+m}(r)}{g_{l+m,l+m}(r)}}
#' where \eqn{p_l} is the proportions of type \emph{l} points in the pattern,
#' \eqn{g()} is the  pair-correlation function, \emph{m} is another type points.
#'
#' Performs the same job as \code{\link{contrl}}, but for untidy data, i.e., the mark is a data frame. The mark variable (for example, \verb{Species})
#'  contains a focal tree species and other
#' tree species. Another mark variable \verb{V.bistate} is only connected to the focal species.
#'
#' The goal is to test whether the spatial pattern of the focal tree affects one of the state ("yes") of
#' the mark \verb{V.bistate}.
#'
#' Use the function \code{\link{perm.r1fn.mkdf}} to permutate the labels of the
#' focal species for envelope simulation, or use the function \code{\link{rpoints.r1fn.mkdf}}
#' to generate simulations of Complete Spatial Randomness (CSR).
#'
#' @return An object of class "\verb{fv}".
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp
#' @importFrom spatstat.explore pcfcross eval.fv pcfdot pcf
#'
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall CRC, Boca Raton.
#'
#' Raventos, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' @seealso \code{\link{contrl}}, \code{\link{heterotrl}},\code{\link{trig_g}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' nsim<-2499
#' plot(contrl2(ppp.marks.c3))
#' ev.contrl2.ash.c3<-envelope(ppp.marks.c3,contrl2,nsim=nsim,
#'     funargs=list(V.species="Species",
#'     host="Ash",V.bistate="Infected",yes="Yes"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     sp="Ash")))
#' plot(ev.contrl2.ash.c3,ylim=c(0,0.8))
#' }
#'
contrl2<-function(ppp.marks,V.species="Species",
                  host="Ash",V.bistate="Infected",yes="Yes",...,
                  r = NULL,
                  kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
                  correction = c("isotropic", "Ripley", "translate"),
                  divisor = c("r", "d"),
                  ratio = FALSE){

  stopifnot(is.ppp(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]))
  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]
  stopifnot(host %in% marks(ppp.marks)[,V.species] && yes %in% marks(ppp.marks)[,V.bistate])
  stopifnot(nlevels(marks(ppp.marks)[1,V.bistate])==2)

  ppp.a<-subset(ppp.marks,marks(ppp.marks)[V.species]==host)
  marks(ppp.a)<-marks(ppp.a)[,V.bistate]

  intens<-intensity(ppp.a)[yes]/sum(intensity(ppp.a))
  pcf.al<-pcfdot(ppp.a,yes,...,
                 r = r,
                 kernel = kernel, bw = bw, stoyan = stoyan,
                 correction = correction,
                 divisor = divisor,
                 ratio = ratio)
  pcf.alm<-pcf(ppp.a,...)

  trl<-eval.fv(intens*pcf.al/pcf.alm)

  host.yes<-paste0(host,".",yes)

  attr(trl,"yexp")<-substitute(p[list(yes, host)](r),list(yes=host.yes,host=host))
  attr(trl,"fname")<-c("p", paste("list(",host.yes,",",host,")"))
  labl.name<-names(trl)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(trl,"labl")<-labls

  return(trl)
}

