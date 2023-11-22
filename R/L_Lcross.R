#'
#' @title \verb{Lcross} difference estimation, tidy data
#'
#' @description
#' Estimates difference between two cross-type \emph{L} functions (\verb{Lcross}).
#'
#' @usage
#' L_Lcross(ppp.marks,yes="Yes",no="No",t.type="yn-nn",...,
#'     r=NULL, breaks=NULL, correction,
#'     ratio=FALSE)
#'
#' @param ppp.marks A "\verb{ppp}" object with a dichotomous mark (for example, "Yes" and "No").
#' @param yes The value of successful category in the dichotomous mark, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category in the dichotomous mark, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param t.type Character, indicating what type of test to be estimated. The value can be "yn-nn", "ny-nn", "yy-nn", and "ny-yn",see \verb{Details}
#' @param ... Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Lcross()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{Kcross()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference \emph{L} functions is estimated as (Raventós, et al., 2010).
#'
#' if t.type=="yn-nn", \deqn{L_{m,l}(r)-L_{l,l}(r)}
#' if t.type=="ny-nn", \deqn{L_{l,m}(r)-L_{l,l}(r)}
#' if t.type=="yy-nn", \deqn{L_{m,m}(r)-L_{l,l}(r)}
#' if t.type=="ny-yn", \deqn{L_{l,l}(r)-L_{m,l}(r)}
#' where \eqn{L()} is the \emph{L}  function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests, for example, the host density-dependent infection of a
#' forest disease.
#'
#' Performs the same job as \code{\link{L_Lcross2}}, but used for tidy data,
#' i.e, with only one mark.
#'
#' Use the function \verb{rlabel()} from the package \verb{spatstat} to permute
#' the binary mark for envelope simulation.
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
#' @seealso \code{\link{L_Lcross2}} which does the same job,
#' \code{\link{g_gcross}}, \code{\link{K_Kcross}}, \code{\link{J_Jcross}},\code{\link{ND_NDcross}},
#' \code{\link{triKcross}}, \code{\link{triL_L}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(L_Lcross(ppp.ash.c3,t.type="yy-nn"))
#' nsim=2499
#' ev.l_lcross.ash.c3<-envelope(ppp.ash.c3,L_Lcross,nsim=nsim,
#'     funargs=list(yes="Yes",no="No"),simulate=rlabel)
#' plot(ev.l_lcross.ash.c3)
#' }
#'

L_Lcross<-function(ppp.marks,yes="Yes",no="No",t.type="yn-nn",...,
                   r=NULL, breaks=NULL, correction,
                   ratio=FALSE){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]
  if(missing(no)) no<-levels(marks(ppp.marks))[1]
  stopifnot(all(c(yes, no ) %in% levels(marks(ppp.marks))))
  stopifnot(yes!=no)

  switch(t.type,"yn-nn"={
    l.yn<-Kcross(ppp.marks,yes,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l.nn<-Kcross(ppp.marks,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l_l<-eval.fv(l.yn-l.nn)
    attr(l_l,"fname")<-c("L", paste("list(",yes,"-", no,",", no,")"))

  },"ny-nn"={
    l.ny<-Kcross(ppp.marks,no,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l.nn<-Kcross(ppp.marks,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l_l<-eval.fv(l.ny-l.nn)
    attr(l_l,"fname")<-c("L", paste("list(",no,",", yes,"-", no,")"))

  }, "yy-nn"={
    l.yy<-Kcross(ppp.marks,yes,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l.nn<-Kcross(ppp.marks,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l_l<-eval.fv(l.yy-l.nn)
    attr(l_l,"fname")<-c("L", paste("list(",yes,"-", yes,",", no,"-", no,")"))

  },"ny-yn"={
    l.ny<-Kcross(ppp.marks,no,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l.yn<-Kcross(ppp.marks,yes,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l_l<-eval.fv(l.ny-l.yn)
    attr(l_l,"fname")<-c("L", paste("list(",no,"-", yes,",", yes,"-", no,")"))

  },stop("Unsupported test type!")
  )

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
