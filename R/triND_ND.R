#'
#' @title Trivariate type \verb{Kcross} difference estimation
#'
#' @usage
#' triND_ND(ppp.marks,V.species="Species",host="Ash", nonhost="Cork",
#'     V.bistate="Infected",yes="Yes",no="No",...,
#'     r = NULL, breaks = NULL,
#'     k=1, correction = c("rs", "km", "han"))
#'
#' @description
#' Estimates difference between two trivariate type Kcross functions
#'
#' @param ppp.marks A marked "\verb{ppp}" object with two qualitative marks, one is multiple and the other is binary.
#' @param V.species The multiple qualitative mark variable
#' @param host One of levels of \verb{V.species}, for example, a host tree species. If \verb{missing}, the 1st level of \verb{V.species} will be assigned.
#' @param nonhost One of levels of \verb{V.species}, for example, a nonhost tree species. If \verb{missing}, the 2nd level of \verb{V.species} will be assigned.
#' @param V.bistate A binary mark only connected to \verb{host}.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param ... Same as in as in \verb{Kcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param k Integer, indicating the \emph{k}th nearest neighbour.
#' @param correction Same as in \verb{Kcross()} from the package \verb{spatstat}
#'
#' @details
#' Trivariate difference K functions is estimated as (Raventós, et al., 2010).
#' \deqn{ND_{l,a}(r)-ND_{m,a}(r)}
#' where \eqn{ND()} is the \verb{NDcross} function, \emph{l} and \emph{m} are
#' two levels of a binary mark of the focal pattern, \emph{a} represents an antecedent pattern.
#'
#' For example, this function tests whether the infected host trees have more neighboring nonhost trees than uninfected host trees.
#'
#' Use the function \code{\link{rpoints.r1fn.mkdf}} to simulate a CSR pattern of the
#' antecedent pattern for envelope simulation, or use the function \code{\link{perm.r1fn.mkdf}} to
#' permute \verb{V.bistate} mark of the focal species.
#'
#' @return An object of "\verb{fv}".
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
#' @seealso \code{\link{triNDcross}},\code{\link{NDcross}},\code{\link{triK_K}},
#' \code{\link{triL_L}},\code{\link{trig_g}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(triND_ND(ppp.marks.c3))
#' nsim=2499
#' ev.trind_nd.cork.c3<-envelope(ppp.marks.c3,triND_ND,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash",nonhost="Cork"),
#'     simulate = expression(rpoints.r1fn.mkdf(ppp.marks.c3,V.species="Species",sp="Cork")))
#' plot(ev.trind_nd.cork.c3)
#' }
#'
triND_ND<-function(ppp.marks,V.species="Species",host="Ash", nonhost="Cork",
                 V.bistate="Infected",yes="Yes",no="No",...,
                 r = NULL, breaks = NULL,
                 k=1, correction = c("rs", "km", "han")){

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

  ppp.al<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                      marks(ppp.marks)[V.species]==host &
                      marks(ppp.marks)[V.bistate]==yes,]
  marks(ppp.al)<-factor(marks(ppp.al)[,V.species])

  ppp.am<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                      marks(ppp.marks)[V.species]==host &
                      marks(ppp.marks)[V.bistate]==no,]
  marks(ppp.am)<-factor(marks(ppp.am)[,V.species])

  k.al<-NDcross(ppp.al,host,nonhost,...,
                r=r, breaks=breaks, k=k, correction=correction)
  k.am<-NDcross(ppp.am,host,nonhost,...,
                r=r, breaks=breaks,k=k, correction=correction)

  nd_nd<-eval.fv(k.al-k.am)

  attr(nd_nd,"yexp")<-substitute({ND[list(yes, nonhost)]^{k}}(r)-{ND[list(no, nonhost)]^{k}}(r),
                                 list(k=k,yes=yes, no=no, nonhost=nonhost))
  attr(nd_nd,"fname")<-c("ND", paste("list(",yes,"-", no,",", nonhost,")"))
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
