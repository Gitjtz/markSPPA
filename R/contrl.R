#'
#' @title Conspecific trivariate random labeling, tidy data
#'
#' @description
#' Trivariate random labeling function for one variate with a binary qualitative mark.
#'
#' @usage
#' contrl(ppp.marks,yes="Yes",no="No",...,
#'     r = NULL,
#'     kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
#'     correction = c("isotropic", "Ripley", "translate"),
#'     divisor = c("r", "d"),ratio = FALSE)
#'
#' @param ppp.marks A multitype "\verb{ppp}" object, of which the binary levels are
#' yes and no.
#' @param yes The value of successful category in the dichotomous mark,
#' for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category in the dichotomous mark,
#' for example, "No". If \verb{missing}, the first level of the mark will be assigned.
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
#' conspecific host tree on infection (Raventos et al. 2010).
#' \deqn{p_{l,l+m}(r)=p_l\frac{g_{l,l+m}(r)}{g_{l+m,l+m}(r)}}
#' where \eqn{p_l} is the proportions of type \emph{l} points in the pattern,
#' \eqn{g()} is the  pair-correlation function, \emph{m} is another type points.
#'
#' Performs the same job as \code{\link{contrl2}}, but for tidy data, i.e, with only one mark.
#'
#' For envelope simulation, let \verb{simulate=NULL} in the
#' function \verb{envelope} from the package \verb{spatstat} to generate simulations of
#' Complete Spatial Randomness (CSR). Alternatively,
#' use \verb{rlabel} from the package \verb{spatstat} to permutate labels.
#'
#' @return An object of class "\verb{fv}".
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity marks<- is.ppp
#' @importFrom spatstat.explore pcfcross eval.fv pcfdot pcf
#'
#' @export
#'
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall CRC, Boca Raton.
#'
#' Raventos, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' @seealso \code{\link{contrl2}}, \code{\link{heterotrl}},\code{\link{trig_g}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' nsim<-2499
#' plot(contrl(ppp.ash.c3,yes="Yes",no="No"))
#' ev.contrl.ash.c3<-envelope(ppp.ash.c3,contrl,nsim=nsim,
#'     funargs=list(yes="Yes",no="No"),
#'     simulate=rlabel)
#' plot(ev.contrl.ash.c3)
#' }
#'
contrl<-function(ppp.marks,yes="Yes",no="No",...,
                 r = NULL,
                 kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
                 correction = c("isotropic", "Ripley", "translate"),
                 divisor = c("r", "d"),
                 ratio = FALSE) {

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]
  if(missing(no)) no<-levels(marks(ppp.marks))[1]
  stopifnot(all(c(yes, no) %in% levels(marks(ppp.marks))))
  stopifnot(yes!=no)

  intens<-intensity(ppp.marks)[yes]/sum(intensity(ppp.marks))
  pcf.al<-pcfdot(ppp.marks,yes,...,
                 r = r,
                 kernel = kernel, bw = bw, stoyan = stoyan,
                 correction = correction,
                 divisor = divisor,
                 ratio = ratio)
  pcf.alm<-pcf(ppp.marks,...)

  trl<-eval.fv(intens*pcf.al/pcf.alm)

  attr(trl,"yexp")<-substitute(p[list(yes, host)](r),list(yes=yes,host=~ symbol("\xb7")))
  attr(trl,"fname")<-c("p", paste("list(",yes,",","~symbol(\"\\267\")",")"))

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
