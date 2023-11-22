#'
#' @title \verb{Jcross} difference estimation, tidy data
#'
#' @description
#' Estimates difference between two cross-type \emph{J} functions (\verb{Jcross}).
#'
#' @usage
#' J_Jcross(ppp.marks,yes="Yes",no="No",t.type="yn-nn",...,
#'     eps=NULL, r=NULL, breaks=NULL, correction=NULL)
#'
#' @param ppp.marks A "\verb{ppp}" object with a dichotomous mark (for example, "Yes" and "No").
#' @param yes The value of successful category in the dichotomous mark,
#' for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category in the dichotomous mark,
#' for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param t.type Character, indicating what type of test to be estimated. The value can be "yn-nn", "ny-nn", "yy-nn", and "ny-yn",see \verb{Details}
#' @param ... Same as in \verb{Jcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{Jcross()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Jcross()} from the package \verb{spatstat}
#' @param eps Same as in \verb{Jcross()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Jcross()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference \emph{J} functions is estimated as (Raventós, et al., 2010).
#'
#' if t.type=="yn-nn", \deqn{J_{m,l}(r)-J_{l,l}(r)}
#' if t.type=="ny-nn", \deqn{J_{l,m}(r)-J_{l,l}(r)}
#' if t.type=="yy-nn", \deqn{J_{m,m}(r)-J_{l,l}(r)}
#' if t.type=="ny-yn", \deqn{J_{l,l}(r)-J_{m,l}(r)}
#' where \eqn{J()} is the \emph{J} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests, for example, the host density-dependent infection of a
#' forest disease.
#'
#' Performs the same job as \code{\link{J_Jcross2}}, but used for tidy data,
#' i.e, with only one mark.
#'
#' Use the function \verb{rlabel()} from the package \verb{spatstat} to permute
#' the binary mark for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Jcross Kdot eval.fv
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' @seealso \code{\link{J_Jcross2}} which does the same job,
#' \code{\link{K_Kcross}}, \code{\link{L_Lcross}}, \code{\link{g_gcross}},\code{\link{ND_NDcross}},
#' \code{\link{triJcross}}, \code{\link{triJ_J}}.
#'
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(J_Jcross(ppp.ash.c3,t.type="ny-nn"))
#' plot(J_Jcross(ppp.ash.c3,t.type="ny-nn",correction="Ripley"))
#' nsim=19
#' ev.j_jcross.ash.c3<-envelope(ppp.ash.c3,J_Jcross,nsim=nsim,
#'     funargs=list(yes="Yes",no="No"),simulate=rlabel)
#' plot(ev.j_jcross.ash.c3)
#' }
#'

J_Jcross<-function(ppp.marks,yes="Yes",no="No",t.type="yn-nn",...,
                   eps=NULL, r=NULL, breaks=NULL, correction=NULL){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]
  if(missing(no)) no<-levels(marks(ppp.marks))[1]
  stopifnot(all(c(yes, no ) %in% levels(marks(ppp.marks))))
  stopifnot(yes!=no)

  switch(t.type,"yn-nn"={
    k.yn<-Jcross(ppp.marks,yes,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    k.nn<-Jcross(ppp.marks,no,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    j_j<-eval.fv(k.yn-k.nn)
    attr(j_j,"fname")<-c("K", paste("list(",yes,"-", no,",", no,")"))

  },"ny-nn"={
    k.ny<-Jcross(ppp.marks,no,yes,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    k.nn<-Jcross(ppp.marks,no,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    j_j<-eval.fv(k.ny-k.nn)
    attr(j_j,"fname")<-c("K", paste("list(",no,",", yes,"-", no,")"))

  }, "yy-nn"={
    k.yy<-Jcross(ppp.marks,yes,yes,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    k.nn<-Jcross(ppp.marks,no,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    j_j<-eval.fv(k.yy-k.nn)
    attr(j_j,"fname")<-c("K", paste("list(",yes,"-", yes,",", no,"-", no,")"))

  },"ny-yn"={
    k.ny<-Jcross(ppp.marks,no,yes,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    k.yn<-Jcross(ppp.marks,yes,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    j_j<-eval.fv(k.ny-k.yn)
    attr(j_j,"fname")<-c("J", paste("list(",no,"-", yes,",", yes,"-", no,")"))

  },stop("Unsupported test type!")
  )

  labl.name<-names(j_j)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(j_j,"labl")<-labls

 return(j_j)
}
