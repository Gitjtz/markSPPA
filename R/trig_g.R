#'
#' @title Trivariate type \verb{pcfcross} difference estimation
#'
#' @usage
#' trig_g(ppp.marks,V.species="Species",host="Ash", nonhost="Cork",
#'     V.bistate="Infected",yes="Yes",no="No",...,
#'     r = NULL,
#'     kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
#'     correction = c("isotropic", "Ripley", "translate"),
#'     divisor = c("r", "d"),
#'     ratio = FALSE)
#'
#' @description
#' Estimates difference between two trivariate type pcfcross functions
#'
#' @param ppp.marks A marked "\verb{ppp}" object with two qualitative marks, one is multiple and the other is binary.
#' @param V.species The multiple qualitative mark variable
#' @param host One of levels of \verb{V.species}, for example, a host tree species. If \verb{missing}, the 1st level of \verb{V.species} will be assigned.
#' @param nonhost One of levels of \verb{V.species}, for example, a nonhost tree species. If \verb{missing}, the 2nd level of \verb{V.species} will be assigned.
#' @param V.bistate A binary mark only connected to \verb{host}.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{pcfcross()}from the package \verb{spatstat}
#' @param kernel Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param bw Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param stoyan Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param correction Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param divisor Same as in \verb{pcfcross()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{pcfcross()} from the package \verb{spatstat}
#'
#'
#' @details
#' Trivariate difference K functions is estimated as (Raventós, et al., 2010).
#' \deqn{g_{l,a}(r)-g_{m,a}(r)}
#' where \eqn{g()} is the pair-correlation function, \emph{l} and \emph{m} are
#' two levels of a binary mark of the focal pattern, \emph{a} represents an antecedent pattern.
#'
#' For example, this function tests whether the infected host trees have more neighboring nonhost trees than uninfected host trees.
#'
#' Use the function \code{\link{rpoints.r1fn.mkdf}} to simulate a CSR pattern of the
#' antecedent pattern for envelope simulation, or use the function \code{\link{perm.r1fn.mkdf}} to
#' permute \verb{V.bistate} mark of the focal species.
#'
#' @return A "\verb{fv}" object
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
#' @seealso \code{\link{triK_K}},\code{\link{g_gcross}},\code{\link{g_gdot}},
#' \code{\link{govergdot}}, \code{\link{trigcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(trig_g(ppp.marks.c3))
#' nsim=2499
#' ev.trig_g.cork.c3<-envelope(ppp.marks.c3,trig_g,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash",nonhost="Cork"),
#'     simulate = expression(rpoints.r1fn.mkdf(ppp.marks.c3,V.species="Species",sp="Cork")))
#' plot(ev.trig_g.cork.c3)
#' }
#'
trig_g<-function(ppp.marks,V.species="Species",host="Ash", nonhost="Cork",
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
  if(missing(nonhost)) nonhost<-levels(marks(ppp.marks)[1,V.species])[2]
  stopifnot(all(c(host,nonhost) %in% marks(ppp.marks)[,V.species]))
  stopifnot(length(levels(marks(ppp.marks)[1,V.bistate]))==2)

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  ppp.al<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                      marks(ppp.marks)[V.species]==host &
                      marks(ppp.marks)[V.bistate]==yes,]


  marks(ppp.al)<-factor(marks(ppp.al)[,V.species])

  ppp.am<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                      marks(ppp.marks)[V.species]==host &
                      marks(ppp.marks)[V.bistate]==no,]

  marks(ppp.am)<-factor(marks(ppp.am)[,V.species])

  pcf.al<-pcfcross(ppp.al,host,nonhost,...,
                   r = r,
                   kernel = kernel, bw = bw, stoyan = stoyan,
                   correction = correction,
                   divisor = divisor,
                   ratio = ratio)
  pcf.am<-pcfcross(ppp.am,host,nonhost,...,
                   r = r,
                   kernel = kernel, bw = bw, stoyan = stoyan,
                   correction = correction,
                   divisor = divisor,
                   ratio = ratio)

  g_g<-eval.fv(pcf.al-pcf.am)

  attr(g_g,"yexp")<-substitute(g[list(yes, nonhost)](r)-g[list(no, nonhost)](r),
                               list(yes=yes, no=no, nonhost=nonhost))
  attr(g_g,"fname")<-c("g", paste("list(",yes,"-", no,",", nonhost,")"))

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
