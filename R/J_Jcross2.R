#'
#' @title \verb{Jcross} difference estimation, untidy data
#'
#' @description
#' Estimates difference between two cross-type \verb{J} functions (\verb{Jcross}).
#'
#' @usage
#' J_Jcross2(ppp.marks,V.species="Species", host="Ash", V.bistate="Infected",
#'     yes="Yes",no="No",t.type="yn-nn",...,
#'     eps=NULL, r=NULL, breaks=NULL, correction=NULL)
#'
#' @param ppp.marks A "\verb{ppp}" object with at least two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark (for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param t.type Character, indicating what type of test to be estimated. The value can be "yn-nn", "ny-nn", "yy-nn", and "ny-yn",see \verb{Details}
#' @param ... Same as in \verb{Jcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{Jcross()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Jcross()} from the package \verb{spatstat}
#' @param eps Same as in \verb{Jcross()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Jcross()} from the package \verb{spatstat}
#'
#' @details
#' Bivariate difference \emph{J} function is estimated as (Raventós, et al., 2010).
#'
#' if t.type=="yn-nn", \deqn{J_{m,l}(r)-J_{l,l}(r)}
#' if t.type=="ny-nn", \deqn{J_{l,m}(r)-J_{l,l}(r)}
#' if t.type=="yy-nn", \deqn{J_{m,m}(r)-J_{l,l}(r)}
#' if t.type=="ny-yn", \deqn{J_{l,l}(r)-J_{m,l}(r)}
#' where \eqn{J()} is the \verb{J} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' Performs the same job as \code{\link{J_Jcross}}, but for untidy data,
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
#' @seealso \code{\link{J_Jcross}} which does the same job,
#' \code{\link{K_Kcross2}}, \code{\link{L_Lcross2}}, \code{\link{g_gcross2}},\code{\link{ND_NDcross2}},
#' \code{\link{triJcross}}, \code{\link{triJ_J}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(J_Jcross2(ppp.marks.c3,t.type="ny-yn"))
#' nsim=2499
#' ev.j_jcross2.ash.c3<-envelope(ppp.marks.c3,J_Jcross2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.j_jcross2.ash.c3)
#' }
#'
J_Jcross2<-function(ppp.marks,V.species="Species", host="Ash",
                    V.bistate="Infected",yes="Yes",no="No",t.type="yn-nn",...,
                    eps=NULL, r=NULL, breaks=NULL, correction=NULL){

  stopifnot(is.ppp(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]) && is.factor(marks(ppp.marks)[,V.bistate]))
  stopifnot(host %in% marks(ppp.marks)[,V.species] && length(levels(marks(ppp.marks)[1,V.bistate]))==2)

  if(missing(host)) host<-levels(marks(ppp.marks)[1,V.species])[1]

  ppp.lm<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  marks(ppp.lm)<-marks(ppp.lm)[V.bistate]

  if(missing(yes)) yes<-levels(marks(ppp.marks)[1,V.bistate])[2]
  if(missing(no)) no<-levels(marks(ppp.marks)[1,V.bistate])[1]

  switch(t.type,"yn-nn"={
    k.yn<-Jcross(ppp.lm,yes,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    k.nn<-Jcross(ppp.lm,no,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    j_j<-eval.fv(k.yn-k.nn)
    attr(j_j,"fname")<-c("K", paste("list(",yes,"-", no,",", no,")"))

  },"ny-nn"={
    k.ny<-Jcross(ppp.lm,no,yes,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    k.nn<-Jcross(ppp.lm,no,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    j_j<-eval.fv(k.ny-k.nn)
    attr(j_j,"fname")<-c("K", paste("list(",no,",", yes,"-", no,")"))

  }, "yy-nn"={
    k.yy<-Jcross(ppp.lm,yes,yes,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    k.nn<-Jcross(ppp.lm,no,no,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    j_j<-eval.fv(k.yy-k.nn)
    attr(j_j,"fname")<-c("K", paste("list(",yes,"-", yes,",", no,"-", no,")"))

  },"ny-yn"={
    k.ny<-Jcross(ppp.lm,no,yes,...,
                 eps=eps, r=r, breaks=breaks, correction=correction)
    k.yn<-Jcross(ppp.lm,yes,no,...,
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
