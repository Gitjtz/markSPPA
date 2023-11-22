#'
#' @title \verb{NDcross} difference estimation, tidy data
#'
#' @description
#' Estimates difference between two cross-type \emph{ND} functions (\verb{NDcross}).
#'
#' @usage
#' ND_NDcross(ppp.marks,yes="Yes",no="No",t.type="yn-nn",
#'    r = NULL, breaks = NULL, ...,
#'    k=1, correction = c("rs", "km", "han"))
#'
#' @param ppp.marks A "\verb{ppp}" object with a dichotomous mark (for example, "Yes" and "No").
#' @param yes The value of successful category in the dichotomous mark,
#' for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category in the dichotomous mark,
#' for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param t.type Character, indicating what type of test to be estimated. The value can be "yn-nn", "ny-nn", "yy-nn", and "ny-yn",see \verb{Details}
#' @param ... Same as in \verb{NDcross()}
#' @param r Same as in \verb{NDcross()}
#' @param breaks Same as in \verb{NDcross()}
#' @param k Integer, indicating the \emph{k}th nearest neighbour.
#' @param correction Same as in \verb{NDKcross()}
#'
#' @details
#' Bivariate difference \emph{ND} functions is estimated as (Raventós, et al., 2010).
#'
#' if t.type=="yn-nn", \deqn{ND_{m,l}(r)-ND_{l,l}(r)}
#' if t.type=="ny-nn", \deqn{ND_{l,m}(r)-ND_{l,l}(r)}
#' if t.type=="yy-nn", \deqn{ND_{m,m}(r)-ND_{l,l}(r)}
#' if t.type=="ny-yn", \deqn{ND_{l,l}(r)-ND_{m,l}(r)}
#' where \eqn{ND()} is the \emph{k}th nearest neighbor distance function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' This function tests, for example, the host density-dependent infection of a
#' forest disease.
#'
#' Performs the same job as \code{\link{ND_NDcross2}}, but used for tidy data,
#' i.e, with only one mark.
#'
#' Use the function \verb{rlabel()} from the package \verb{spatstat} to permute
#' the binary mark for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked as.ppp
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#'
#' @seealso \code{\link{ND_NDcross2}} which does the same job,
#' \code{\link{g_gcross}}, \code{\link{K_Kcross}}, \code{\link{J_Jcross}},\code{\link{L_Lcross}},
#' \code{\link{triNDcross}}, \code{\link{triND_ND}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(ND_NDcross(amacrine,yes="off",no="on",t.type="ny-nn"))
#' plot(ND_NDcross(ppp.ash.c3,t.type="ny-yn"))
#' nsim=2499
#' ev.nd_ndcross.ash.c3<-envelope(ppp.ash.c3,ND_NDcross,nsim=nsim,
#'     funargs=list(yes="Yes",no="No"),simulate=rlabel)
#' plot(ev.nd_ndcross.ash.c3)
#' }

ND_NDcross<-function(ppp.marks,yes="Yes",no="No",t.type="yn-nn",
                     r = NULL, breaks = NULL, ...,
                     k=1, correction = c("rs", "km", "han")){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]
  if(missing(no)) no<-levels(marks(ppp.marks))[1]
  stopifnot(all(c(yes, no ) %in% levels(marks(ppp.marks))))
  stopifnot(yes!=no)

  switch(t.type,"yn-nn"={
    nd.yn<-NDcross(ppp.marks,yes,no,...,
                 r=r, breaks=breaks, correction=correction,
                 k=k)
    nd.nn<-NDcross(ppp.marks,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 k=k)
    nd_nd<-eval.fv(nd.yn-nd.nn)
    attr(nd_nd,"fname")<-c("ND", paste("list(",yes,"-", no,",", no,")"))

  },"ny-nn"={
    nd.ny<-NDcross(ppp.marks,no,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 k=k)
    nd.nn<-NDcross(ppp.marks,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 k=k)
    nd_nd<-eval.fv(nd.ny-nd.nn)
    attr(nd_nd,"fname")<-c("ND", paste("list(",no,",", yes,"-", no,")"))

  }, "yy-nn"={
    nd.yy<-NDcross(ppp.marks,yes,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 k=k)
    nd.nn<-NDcross(ppp.marks,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 k=k)
    nd_nd<-eval.fv(nd.yy-nd.nn)
    attr(nd_nd,"fname")<-c("ND", paste("list(",yes,"-", yes,",", no,"-", no,")"))

  },"ny-yn"={
    nd.ny<-NDcross(ppp.marks,no,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 k=k)
    nd.yn<-NDcross(ppp.marks,yes,no,...,
                 r=r, breaks=breaks, correction=correction,
                 k=k)
    nd_nd<-eval.fv(nd.ny-nd.yn)
    attr(nd_nd,"fname")<-c("ND", paste("list(",no,"-", yes,",", yes,"-", no,")"))

  },stop("Unsupported test type!")
  )

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
