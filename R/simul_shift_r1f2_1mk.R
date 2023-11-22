#'
#' @title Randomly shifts the points corresponding to one level of a qualitative mark
#'
#' @description
#' Given a "ppp" object with a multitype qualitative mark variable, this function
#' does a shift of the points corresponding to one level of the mark.
#'
#' @usage simul.shift.r1f2.1mk(ppp.marks,yes="Yes",...,
#'     which=NULL, group, nsim=1, drop=TRUE)
#'
#' @param ppp.marks A "\verb{ppp}" object with a multitype qualitative mark variable
#' @param yes The level of which points to be randomized. If \verb{missing}, the second level of the mark will be assigned.
#' @param ... Same as in \verb{rshift()} from the package \verb{spatstat}.
#' @param which Same as in \verb{rshift()} from the package \verb{spatstat}.
#' @param group Same as in \verb{rshift()} from the package \verb{spatstat}.
#' @param nsim Same as in \verb{rshift()} from the package \verb{spatstat}.
#' @param drop Same as in \verb{rshift()} from the package \verb{spatstat}.
#'
#' @details
#' This function activates the function \verb{rshift} from the package \verb{spatstat}.
#' For a multitype ppp, this function only shifts the points corresponding to pattern \verb{yes}.
#'
#'
#' @return A "\verb{ppp}" object.
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv envelope
#' @importFrom spatstat.geom Window
#' @importFrom spatstat.random rpoint rlabel rpoispp rmpoispp runifpoint rshift rshift.ppp
#'
#' @export
#'
#' @seealso \code{\link{perm.r1r2.mkdf}}, \code{\link{perm.onemkvar}}, \code{\link{perm.r1fn.mkdf}},
#' \code{\link{rpoints.r1fn.mkdf}}, \code{\link{rpoints.r1f2.1mk}},
#'  \code{\link{simul.shift.r1fn.mkdf}}, \code{\link{simul.jitter.r1f2.1mk}}, \code{\link{simul.jitter.r1fn.mkdf}},
#' \code{\link{trig_g}},\code{\link{triK_K}}, \code{\link{heterotrl}}, \code{\link{contrl2}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' simul.shift.r1f2.1mk(amacrine,radius=0.1,edge="erode")
#' }
#'
simul.shift.r1f2.1mk<-function(ppp.marks,yes="Yes",...,
                               which=NULL, group, nsim=1, drop=TRUE){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]

  stopifnot(yes %in% marks(ppp.marks))

  ppp.marks.sp<-ppp.marks[marks(ppp.marks)==yes,]
  ppp.marks.nosp<-ppp.marks[marks(ppp.marks)!=yes,]

  simulated.sp<-rshift(ppp.marks.sp,...,
                       which=which, group=group, nsim=nsim, drop=drop)

  marks(simulated.sp)<-factor(rep(yes,npoints(simulated.sp)))

  simulated.species<-superimpose(simulated.sp,ppp.marks.nosp)

  return(simulated.species)
}

