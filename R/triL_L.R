#'
#' @title Trivariate type \verb{Lcross} difference estimation
#'
#' @usage
#' triL_L(ppp.marks,V.species="Species",host="Ash", nonhost="Cork",
#'     V.bistate="Infected",yes="Yes",no="No",...,
#'     r=NULL, breaks=NULL, correction,
#'     ratio=FALSE)
#'
#' @description
#' Estimates difference between two trivariate type \emph{L} functions
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
#' @param correction Same as in \verb{Lcross()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{Kcross()} from the package \verb{spatstat}
#'
#' @details
#' Trivariate difference \emph{L} functions is estimated as
#' \deqn{L_{l,a}(r)-L_{m,a}(r)}
#' where \eqn{L()} is the \emph{L}  function, \emph{l} and \emph{m} are
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
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot Lcross Ldot eval.fv
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
#' @seealso \code{\link{trig_g}},\code{\link{K_Kcross}},\code{\link{K_Kdot}},
#' \code{\link{KoverKdot}},\code{\link{triK_K}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(triL_L(ppp.marks.c3))
#' nsim=2499
#' ev.tril_l.cork.c3<-envelope(ppp.marks.c3,triL_L,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash",nonhost="Cork"),
#'     simulate = expression(rpoints.r1fn.mkdf(ppp.marks.c3,V.species="Species",sp="Cork")))
#' plot(ev.tril_l.cork.c3)
#' }
#'
triL_L<-function(ppp.marks,V.species="Species",host="Ash", nonhost="Cork",
                 V.bistate="Infected",yes="Yes",no="No",...,
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
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  ppp.al<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                      marks(ppp.marks)[V.species]==host &
                      marks(ppp.marks)[V.bistate]==yes,]
  marks(ppp.al)<-factor(marks(ppp.al)[,V.species])

  ppp.am<-ppp.marks[marks(ppp.marks)[V.species]==nonhost |
                      marks(ppp.marks)[V.species]==host &
                      marks(ppp.marks)[V.bistate]==no,]
  marks(ppp.am)<-factor(marks(ppp.am)[,V.species])

  k.al<-Lcross(ppp.al,host,nonhost,...,
               r=r, breaks=breaks, correction=correction,
               ratio=ratio)
  k.am<-Lcross(ppp.am,host,nonhost,...,
               r=r, breaks=breaks, correction=correction,
               ratio=ratio)

  l_l<-eval.fv(k.al-k.am)

  attr(l_l,"yexp")<-substitute(L[list(yes, nonhost)](r)-L[list(no, nonhost)](r),list(yes=yes, no=no, nonhost=nonhost))
  attr(l_l,"fname")<-c("L", paste("list(",yes,"-", no,",", nonhost,")"))
  labl.name<-names(l_l)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(l_l,"labl")<-labls

  return(l_l)

}
