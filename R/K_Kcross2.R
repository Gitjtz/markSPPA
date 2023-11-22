#'
#' @title \verb{Kcross} difference estimation, untidy data
#'
#' @description
#' Estimates difference between two cross-type K functions (Kcross).
#'
#' @usage
#' K_Kcross2(ppp.marks,V.species="Species", host="Ash", V.bistate="Infected",
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
#' @param ... Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param r Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param correction Same as in \verb{Kcross()} from the package \verb{spatstat}
#' @param ratio Same as in \verb{Kcross()} from the package \verb{spatstat}
#'
#'
#' @details
#' Bivariate difference \emph{K}  function is estimated as (Raventós, et al., 2010).
#'
#' if t.type=="yn-nn", \deqn{K_{m,l}(r)-K_{l,l}(r)}
#' if t.type=="ny-nn", \deqn{K_{l,m}(r)-K_{l,l}(r)}
#' if t.type=="yy-nn", \deqn{K_{m,m}(r)-K_{l,l}(r)}
#' if t.type=="ny-yn", \deqn{K_{l,l}(r)-K_{m,l}(r)}
#' where \eqn{K()} is the Ripley \emph{K} function, \emph{l} and \emph{m} are
#' two levels of a binary mark.
#'
#' Note that there is a symmetry in i and j, so that \eqn{K_{ij}(r) = K_{ji}(r)},
#' However, estimates of these functions may not be exactly symmetric in i and j,
#' because of the effect of edge corrections.
#'
#' Performs the same job as \code{\link{K_Kcross}}, but for untidy data,
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
#' Baddeley A.,Rubak E.,Turner R. (2016). Spatial Point Patterns: Methodology and Applications with R.
#' Boca Raton, FL: CRC Press.
#'
#' @seealso \code{\link{K_Kcross}} which does the same job,
#' \code{\link{g_gcross2}}, \code{\link{L_Lcross2}}, \code{\link{J_Jcross2}},\code{\link{ND_NDcross2}},
#' \code{\link{triKcross}}, \code{\link{triK_K}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' plot(K_Kcross2(ppp.marks.c3,t.type="ny-nn"))
#' nsim=2499
#' ev.k_kcross2.ash.c3<-envelope(ppp.marks.c3,K_Kcross2,nsim=nsim,
#'     funargs=list(V.species="Species", host="Ash", V.bistate="Infected"),
#'     simulate=expression(perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",
#'     host="Ash",V.bistate="Infected")))
#' plot(ev.k_kcross2.ash.c3)
#' }
#'
K_Kcross2<-function(ppp.marks,V.species="Species", host="Ash",
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
    k.yn<-Kcross(ppp.lm,yes,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    k.nn<-Kcross(ppp.lm,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    k_k<-eval.fv(k.yn-k.nn)
    attr(k_k,"fname")<-c("K", paste("list(",yes,"-", no,",", no,")"))

  },"ny-nn"={
    k.ny<-Kcross(ppp.lm,no,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    k.nn<-Kcross(ppp.lm,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    k_k<-eval.fv(k.ny-k.nn)
    attr(k_k,"fname")<-c("K", paste("list(",no,",", yes,"-", no,")"))

  }, "yy-nn"={
    k.yy<-Kcross(ppp.lm,yes,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    k.nn<-Kcross(ppp.lm,no,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    k_k<-eval.fv(k.yy-k.nn)
    attr(k_k,"fname")<-c("K", paste("list(",yes,"-", yes,",", no,"-", no,")"))

  },"ny-yn"={
    k.ny<-Kcross(ppp.lm,no,yes,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    k.yn<-Kcross(ppp.lm,yes,no,...,
                 r=r, breaks=breaks, correction=correction,
                 ratio=ratio)
    k_k<-eval.fv(k.ny-k.yn)
    attr(k_k,"fname")<-c("K", paste("list(",no,"-", yes,",", yes,"-", no,")"))

  },stop("Unsupported test type!")
  )

  labl.name<-names(k_k)
  labls<-vector(length=length(labl.name))
  labls[1]<-"r"
  labls[2]<-"{%s[%s]^{Pois}}(r)"
  for(i in 3:length(labl.name)){
    labls[i]<-paste0("{hat(%s)[%s]^{",labl.name[i],"}}(r)")
  }
  attr(k_k,"labl")<-labls

  return(k_k)
}
