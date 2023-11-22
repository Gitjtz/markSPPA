#'
#' @title Permuting two parts of a mark separately
#'
#' @description
#' For a mark dataframe, permuting two parts of a mark separately corresponding to two levels of another mark variable
#'
#' @usage perm.r1r2.mkdf(ppp.marks,m.qli="Species",sp1="Ash", sp2="Birch",
#'     m.qti="Disease",...)
#'
#' @param ppp.marks A "\verb{ppp}" object with at least two mark variable and one of them is binary.
#' @param m.qli A qualitative mark variable. If there is more than two categories, only the first two levels are considered.
#' @param ... further arguments passed to or from other methods.
#' @param sp1 A mark value of \verb{m.qli}, if \verb{missing}, the first level of \verb{m.qli}
#' @param sp2 Another mark value of \verb{m.qli}, if \verb{missing}, the second level of \verb{m.qli}
#' @param m.qti A mark variable of \verb{ppp.marks}, usually quantitative.
#'
#' @details
#' This function permutates \verb{m.qti} corresponding to \verb{sp1} and \verb{sp2} seperately.
#'
#' Given a category variable "Species" has two levels (Ash, Birch) and another mark variable "DBH",
#' i.e., the sizes of tree,  this function randomizes the size of each species separately.
#' For other levels of variable "Species",for example, Cork, their sizes keep unchanged.
#'
#' This function may be used for envelope simulation for data 9.
#'
#' @return A "\verb{ppp}" object.
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#'
#' @export
#'
#' @seealso \code{\link{perm.onemkvar}}, \code{\link{perm.r1fn.mkdf}},
#' \code{\link{rpoints.r1fn.mkdf}}, \code{\link{rpoints.r1f2.1mk}}, \code{\link{simul.shift.r1f2.1mk}},
#' \code{\link{simul.shift.r1fn.mkdf}}, \code{\link{simul.jitter.r1f2.1mk}}, \code{\link{simul.jitter.r1fn.mkdf}},
#' \code{\link{mkcorr.d9.t6}}, \code{\link{mkcorr.d9}}
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' #perm.r1r2.mkdf(ppp.ash.infect.dbh.c3,m.qli="Infected",m.qti="DBH")
#' }
#'
perm.r1r2.mkdf<-function(ppp.marks,m.qli="Species",sp1="Ash", sp2="Birch",
                         m.qti="Disease",...){

  stopifnot(is.ppp(ppp.marks))
  stopifnot(all(c(m.qli, m.qti ) %in% names(marks(ppp.marks))))
  stopifnot(is.factor(ppp.marks$marks[[m.qli]]) && nlevels(ppp.marks$marks[[m.qli]])>=2)

  if(missing(sp1))   sp1<-levels(marks(ppp.marks)[[m.qli]])[1]
  if(missing(sp2))   sp2<-levels(marks(ppp.marks)[[m.qli]])[2]
  ppp.marks.sp1<-subset(ppp.marks,marks(ppp.marks)[[m.qli]]==sp1)
  ppp.marks.sp2<-subset(ppp.marks,marks(ppp.marks)[[m.qli]]==sp2)
  marks(ppp.marks.sp1)[m.qti]<-sample(marks(ppp.marks.sp1)[,m.qti],npoints(ppp.marks.sp1))
  marks(ppp.marks.sp2)[m.qti]<-sample(marks(ppp.marks.sp2)[,m.qti],npoints(ppp.marks.sp2))
  simulate.sp<-superimpose(ppp.marks.sp1,ppp.marks.sp2)
  if(nlevels(marks(ppp.marks)[[m.qli]])>2){
    ppp.marks.sp3<-subset(ppp.marks,marks(ppp.marks)[[m.qli]]!=sp1 & marks(ppp.marks)[[m.qli]]!=sp2)
    simulate.sp<-superimpose(simulate.sp,ppp.marks.sp3)
  }

  return(simulate.sp)
}
