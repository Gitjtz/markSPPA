#'
#' @title Heterospecific trivariate random labeling
#'
#' @description
#' Trivariate random labeling function for data with two qualitative marks.
#'
#' @usage
#' heterotrl(ppp.marks,V.species="Species",
#'     host="Ash",nonhost="Cork",
#'     V.bistate="Infected",yes="Yes",...,
#'     r = NULL,
#'     kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
#'     correction = c("isotropic", "Ripley", "translate"),
#'     divisor = c("r", "d"),ratio = FALSE)
#'
#' @param ppp.marks A "ppp" object with at least two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark(for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A level of \verb{V.species}, indicating the focal tree species. That is the host tree of a disease which is concerned. If missing, the 1st level of \verb{V.species}
#' @param nonhost A level of V.species, indicating one of nonhost tree species. If missing, the 2nd level of \verb{V.species}
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes A level indicates the state of the host tree of concerned, for example, "Yes" or "No". If missing, the 2nd level of \verb{V.bistate}
#' If \verb{missing}, the second level of \verb{V.bistate } will be assigned.
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
#' This function tests the influence of non host tree species on the infection of host species (Raventós et al. 2010).
#' \deqn{p_{l,h}(r)=\frac{p_lg_{l,h}(r)}{g_{l+m,h}(r)}}
#' where \eqn{p_l} is the proportions of type \emph{l} points in the the host pattern,
#' \eqn{g()} is the pair-correlation function, \emph{m} is another host type points,
#' \emph{h} represents the points of non host species.
#'
#' Trivariate random labeling explores the effect of an antecedent focal
#' pattern \emph{h} on the process that distributes a qualitative
#' mark (type \emph{l} and type \emph{m}) .
#' That is, to explore whether the qualitative mark of the second pattern
#' (of a focal species) depends on the distance from a point of the antecedent
#' pattern (of heterospecies).
#'
#' In trivariate case, randomization was only occurred on the antecedent pattern
#' (nonhost). Use the function \code{\link{rpoints.r1fn.mkdf}} to carry out
#' simulation of a random point pattern of a nonhost species. However,the function
#' \code{\link{perm.r1fn.mkdf}} can be used to carry out a permutation of marks of
#' a host species.
#'
#' @return A "\verb{fv}" object.
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot eval.fv
#'
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial segregation hypothesis: a test with nine-year survivorship data in a Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' @seealso \code{\link{contrl}}, \code{\link{contrl2}},\code{\link{trig_g}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(heterotrl(ppp.marks.c3))
#' nsim<-2499
#' ev.heterotrl.cork.c3.csr<-envelope(ppp.marks.c3,heterotrl,nsim=nsim,
#'     funargs=list(V.species="Species",host="Ash",nonhost="Cork",
#'     V.bistate="Infected",yes="Yes"),
#'     simulate=expression(rpoints.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     sp="Cork")))
#' plot(ev.heterotrl.cork.c3.csr)
#' ev.heterotrl.cork.c3.shuffle<-envelope(ppp.marks.c3,heterotrl,nsim=nsim,
#'     funargs=list(V.species="Species",host="Ash",nonhost="Cork",
#'     V.bistate="Infected",yes="Yes"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,
#'     V.species="Species",sp="Ash")))
#' plot(ev.heterotrl.cork.c3.shuffle)
#' }
#'
heterotrl<-function(ppp.marks,V.species="Species",host="Ash",nonhost="Cork",
                    V.bistate="Infected",yes="Yes",...,
                    r = NULL,
                    kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
                    correction = c("isotropic", "Ripley", "translate"),
                    divisor = c("r", "d"),
                    ratio = FALSE){

  stopifnot(is.ppp(ppp.marks) && is.marked(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))
  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]
  if(missing(nonhost)) nonhost<-levels(marks(ppp.marks)[1,V.species])[2]
  stopifnot(all(c(host,nonhost) %in% marks(ppp.marks)[,V.species]) && yes %in% marks(ppp.marks)[,V.bistate])


  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  intens<-intensity(ppp.lm)[yes]/sum(intensity(ppp.lm))
  ppp.al<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                      marks(ppp.marks)[V.species]==host & marks(ppp.marks)[V.bistate]==yes,]
  marks(ppp.al)<-factor(marks(ppp.al)[,V.species])

  ppp.alm<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                       marks(ppp.marks)[V.species]==host,]
  marks(ppp.alm)<-factor(marks(ppp.alm)[,V.species])

  pcf.al<-pcfcross(ppp.al,host,nonhost,...,
                   r = r,
                   kernel = kernel, bw = bw, stoyan = stoyan,
                   correction = correction,
                   divisor = divisor,
                   ratio = ratio)
  pcf.alm<-pcfcross(ppp.alm,host,nonhost,...,
                    r = r,
                    kernel = kernel, bw = bw, stoyan = stoyan,
                    correction = correction,
                    divisor = divisor,
                    ratio = ratio)
  trl<-eval.fv(intens*pcf.al/pcf.alm)

  host.yes<-paste0(host,".",yes)
  attr(trl,"yexp")<-substitute(p[list(host, nonhost)](r),
                               list(host=host.yes,nonhost=nonhost))
  attr(trl,"fname")<-c("p", paste("list(",host.yes,",",nonhost,")"))

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
