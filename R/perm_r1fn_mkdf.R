#'
#' @title Permuting a part of data of a variable
#'
#' @description
#' Randomizes partial data of a mark variable corresponding to one of levels of
#' another qualitative variable
#'
#' @usage
#' perm.r1fn.mkdf(ppp.marks,V.species="Species",host="Ash",
#'     V.bistate="Infected", ...)
#'
#' @param ppp.marks A "\verb{ppp}" object with at least two mark variable and one of them is qualitative.
#' @param V.species A qualitative variable which one level is given by host
#' @param host One of category of V.species
#' @param V.bistate A mark variable for permutation. Often qualitative (binary), but quantitative is also OK.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Randomizes the binary label of a host species.
#' Often used in envelope simulation of trivariate random labeling.
#'
#'
#' @return A "\verb{ppp}" object
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#'
#' @export
#'
#' @seealso \code{\link{perm.r1r2.mkdf}}, \code{\link{perm.onemkvar}},
#' \code{\link{rpoints.r1fn.mkdf}}, \code{\link{rpoints.r1f2.1mk}}, \code{\link{simul.shift.r1f2.1mk}},
#'  \code{\link{simul.shift.r1fn.mkdf}}, \code{\link{simul.jitter.r1f2.1mk}}, \code{\link{simul.jitter.r1fn.mkdf}},
#'   \code{\link{g_gcross2}}, \code{\link{K_Kcross2}},
#' \code{\link{g_gdot2}}, \code{\link{K_Kdot2}}, \code{\link{govergdot2}},
#' \code{\link{KoverKdot2}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' perm.r1fn.mkdf(ppp.marks.c3,V.species="Species",host="Ash",V.bistate="Infected")
#'}
#'
perm.r1fn.mkdf<-function(ppp.marks,V.species="Species",host="Ash",
                        V.bistate="Infected",...){

  stopifnot(is.ppp(ppp.marks) && is.marked(ppp.marks))
  stopifnot(is.data.frame(marks(ppp.marks)))
  stopifnot(all(c(V.species, V.bistate) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]))
  stopifnot(host %in% marks(ppp.marks)[,V.species])

  ppp.marks.species<-ppp.marks[marks(ppp.marks)[V.species]==host,]
  ppp.marks.nospecies<-ppp.marks[marks(ppp.marks)[V.species]!=host,]
  marks(ppp.marks.species)[V.bistate]<-
    sample(marks(ppp.marks.species)[,V.bistate],npoints(ppp.marks.species))
  simulated.shuffle<-superimpose(ppp.marks.species,ppp.marks.nospecies)
  return(simulated.shuffle)
}
