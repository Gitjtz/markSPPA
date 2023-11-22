#'
#' @title \verb{Lcross} difference estimation, untidy data
#'
#' @description
#' Estimates difference between two cross-type  \emph{L} functions (\verb{Lcross}).
#'
#' @usage
#' L_Lcross2(ppp.marks,V.species="Species", host="Ash", V.bistate="Infected",
#'     yes="Yes",no="No",t.type="yn-nn",...,
#'     r=NULL, breaks=NULL, correction,
#'     ratio=FALSE)
#'
#' @param ppp.marks A "\verb{ppp}" object with at least two qualitative marks. One is a multiple category mark (for example, species in a forest) and the other is a dichotomous mark (for example, "Yes" and "No").
#' @param V.species A column in data represents a multiple category mark, such as tree species in a forest.
#' @param host A value of V.species, indicating the focal tree species. That is the host tree of a disease which is concerned.
#' @param V.bistate A column in data represents the bi-state of the host tree, for example, infected or uninfected.
#' @param yes The value of successful category in the mark \verb{V.bistate}, for example, "Yes". If \verb{missing}, the second level of the mark will be assigned.
#' @param no The value of failure category n the mark \verb{V.bistate}, for example, "No". If \verb{missing}, the first level of the mark will be assigned.
#' @param t.type Character, indicating what type of test to be estimated. The value can be "yn-nn", "ny-nn", "yy-nn", and "ny-yn",see \verb{Details}
#' @param ... Same as in \verb{Lcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Lcross()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{Kcross()} from the package \verb{spatstat}
#'
#'
#' @details
#' Bivariate difference \emph{L}  function is estimated as (Raventós, et al., 2010).
#'
#' if t.type=="yn-nn", \deqn{L_{m,l}(r)-L_{l,l}(r)}
#' if t.type=="ny-nn", \deqn{L_{l,m}(r)-L_{l,l}(r)}
#' if t.type=="yy-nn", \deqn{L_{m,m}(r)-L_{l,l}(r)}
#' if t.type=="ny-yn", \deqn{L_{l,l}(r)-L_{m,l}(r)}
#' where \eqn{L()} is the \emph{L} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' Performs the same job as \code{\link{L_Lcross}}, but for untidy data,
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
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' @seealso \code{\link{L_Lcross}} which does the same job,
#' \code{\link{g_gcross}}, \code{\link{K_Kcross}}, \code{\link{J_Jcross}},\code{\link{ND_NDcross}},
#' \code{\link{triKcross}}, \code{\link{triL_L}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(L_Lcross2(ppp.marks.c3,t.type="yy-nn"))
#' nsim=2499
#' ev.k_kcross2.ash.c3<-envelope(ppp.marks.c3,K_Kcross2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.k_kcross2.ash.c3)
#' }
#'
L_Lcross2<-function(ppp.marks,V.species="Species", host="Ash",
                    V.bistate="Infected",yes="Yes",no="No",t.type="yn-nn",...,
                    r=NULL, breaks=NULL, correction,
                    ratio=FALSE){

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
    l.yn<-Kcross(ppp.lm,yes,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l.nn<-Kcross(ppp.lm,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l_l<-eval.fv(l.yn-l.nn)
    attr(l_l,"fname")<-c("L", paste("list(",yes,"-", no,",", no,")"))

  },"ny-nn"={
    l.ny<-Kcross(ppp.lm,no,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l.nn<-Kcross(ppp.lm,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l_l<-eval.fv(l.ny-l.nn)
    attr(l_l,"fname")<-c("L", paste("list(",no,",", yes,"-", no,")"))

  }, "yy-nn"={
    l.yy<-Kcross(ppp.lm,yes,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l.nn<-Kcross(ppp.lm,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l_l<-eval.fv(l.yy-l.nn)
    attr(l_l,"fname")<-c("L", paste("list(",yes,"-", yes,",", no,"-", no,")"))

  },"ny-yn"={
    l.ny<-Kcross(ppp.lm,no,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    l.yn<-Kcross(ppp.lm,yes,no,...,
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
