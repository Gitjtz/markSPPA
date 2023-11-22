#'
#' @title \verb{NDcross} difference estimation, untidy data
#'
#' @description
#' Estimates difference between two cross-type \emph{ND} functions (\verb{NDcross}).
#'
#' @usage
#' ND_NDcross2(ppp.marks,V.species="Species", host="Ash", V.bistate="Infected",
#'     yes="Yes",no="No",t.type="yn-nn",...,
#'     r = NULL, breaks = NULL,
#'     k=1, correction = c("rs", "km", "han"))
#'
#' @param ppp.marks A "\verb{ppp}" object with at least two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark (for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param t.type Character, indicating what type of test to be estimated. The value can be "yn-nn", "ny-nn", "yy-nn", and "ny-yn",see \verb{Details}
#' @param ... Same as in \verb{NDcross()}
#' @param r Same as in \verb{NDcross()}
#' @param breaks Same as in \verb{NDcross()}
#' @param k Integer, indicating the \emph{k}th nearest neighbour.
#' @param correction Same as in \verb{NDcross()}
#'
#' @details
#'
#' if t.type=="yn-nn", \deqn{ND_{m,l}(r)-ND_{l,l}(r)}
#' if t.type=="ny-nn", \deqn{ND_{l,m}(r)-ND_{l,l}(r)}
#' if t.type=="yy-nn", \deqn{ND_{m,m}(r)-ND_{l,l}(r)}
#' if t.type=="ny-yn", \deqn{ND_{l,l}(r)-ND_{m,l}(r)}
#' where \eqn{ND()} is the\emph{k}th nearest neighbor distance function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' Performs the same job as \code{\link{ND_NDcross}}, but for untidy data,
#' i.e., the mark is a data frame. The mark variable (for example, \verb{Species})
#' contains a focal tree species and other tree species.
#' Another mark variable \verb{V.bistate} is only connected to the focal species.
#'
#' This function tests, for example, the host density-dependent infection of a
#' forest disease.
#'
#' Use the function \code{\link{perm.r1fn.mkdf}} to permute \verb{V.bistate} mark of the
#' focal species for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
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
#' @seealso \code{\link{ND_NDcross}} which does the same job,
#' \code{\link{g_gcross}2}, \code{\link{K_Kcross2}}, \code{\link{J_Jcross2}},\code{\link{L_Lcross2}},
#' \code{\link{triNDcross}}, \code{\link{triND_ND}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(ND_NDcross2(ppp.marks.c3,t.type="ny-nn"))
#' nsim=2499
#' ev.ND_NDcross2.ash.c3<-envelope(ppp.marks.c3,ND_NDcross2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(simul.shuffle(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.ND_NDcross2.ash.c3)
#' }
#'
ND_NDcross2<-function(ppp.marks,V.species="Species", host="Ash",
                    V.bistate="Infected",yes="Yes",no="No",t.type="yn-nn",...,
                    r = NULL, breaks = NULL,
                    k=1, correction = c("rs", "km", "han")){

  stopifnot(is.ppp(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))
  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]
  stopifnot(host %in% marks(ppp.marks)[,V.species] && length(levels(marks(ppp.marks)[1,V.bistate]))==2)

  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  switch(t.type,"yn-nn"={
    nd.yn<-NDcross(ppp.lm,yes,no,...,
                   r=r, breaks=breaks, correction=correction,
                   k=k)
    nd.nn<-NDcross(ppp.lm,no,no,...,
                   r=r, breaks=breaks, correction=correction,
                   k=k)
    nd_nd<-eval.fv(nd.yn-nd.nn)
    attr(nd_nd,"fname")<-c("ND", paste("list(",yes,"-", no,",", no,")"))

  },"ny-nn"={
    nd.ny<-NDcross(ppp.lm,no,yes,...,
                   r=r, breaks=breaks, correction=correction,
                   k=k)
    nd.nn<-NDcross(ppp.lm,no,no,...,
                   r=r, breaks=breaks, correction=correction,
                   k=k)
    nd_nd<-eval.fv(nd.ny-nd.nn)
    attr(nd_nd,"fname")<-c("ND", paste("list(",no,",", yes,"-", no,")"))

  }, "yy-nn"={
    nd.yy<-NDcross(ppp.lm,yes,yes,...,
                   r=r, breaks=breaks, correction=correction,
                   k=k)
    nd.nn<-NDcross(ppp.lm,no,no,...,
                   r=r, breaks=breaks, correction=correction,
                   k=k)
    nd_nd<-eval.fv(nd.yy-nd.nn)
    attr(nd_nd,"fname")<-c("ND", paste("list(",yes,"-", yes,",", no,"-", no,")"))

  },"ny-yn"={
    nd.ny<-NDcross(ppp.lm,no,yes,...,
                   r=r, breaks=breaks, correction=correction,
                   k=k)
    nd.yn<-NDcross(ppp.lm,yes,no,...,
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
