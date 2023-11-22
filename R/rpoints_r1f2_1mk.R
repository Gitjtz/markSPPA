#'
#' @title Random simulation of the points corresponding to one level of a qualitative mark
#'
#' @description
#' Given a "ppp" object with a multitype qualitative mark variable, this function
#' does a random simulation of the points corresponding to one level of the mark.
#'
#' @usage rpoints.r1f2.1mk(ppp.marks,p.type="CSR",yes="Yes",...)
#'
#' @param ppp.marks A "\verb{ppp}" object with a multitype qualitative mark variable
#' @param p.type Indicating what type of points to be generated, see \verb{Details}.
#' @param yes The level of which points to be randomized. If \verb{missing}, the second level of the mark will be assigned.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Two type of points are supported, including "CSR" and
#' "in.unif", which generateed complete spatial random (CSR) points and
#' independent uniform random points, respectively. For "CSR", the rows of simulated marks
#' were sampled randomly from the mark data frame corresponding to "Ash"; for "in.unif",
#' the order of the rows of simulated marks did not change.
#'
#'
#' @return A "\verb{ppp}" object.
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv envelope
#' @importFrom spatstat.geom Window
#' @importFrom spatstat.random rpoint rlabel rpoispp rmpoispp runifpoint
#'
#' @export
#'
#' @seealso \code{\link{perm.r1r2.mkdf}}, \code{\link{perm.onemkvar}}, \code{\link{perm.r1fn.mkdf}},
#' \code{\link{rpoints.r1fn.mkdf}}, \code{\link{simul.shift.r1f2.1mk}},
#'  \code{\link{simul.shift.r1fn.mkdf}}, \code{\link{simul.jitter.r1f2.1mk}}, \code{\link{simul.jitter.r1fn.mkdf}},
#'   \code{\link{heterotrl}}, \code{\link{contrl2}},
#' \code{\link{trig_g}},\code{\link{triK_K}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' rpoints.r1f2.1mk(ppp.ash.c3,yes="Yes")
#' }
#'
rpoints.r1f2.1mk<-function(ppp.marks,p.type="CSR",yes="Yes",...){

  stopifnot(is.multitype(ppp.marks))
  if(missing(yes)) yes<-levels(marks(ppp.marks))[2]

  stopifnot(yes %in% marks(ppp.marks))

  ppp.marks.sp<-ppp.marks[marks(ppp.marks)==yes,]
  ppp.marks.nosp<-ppp.marks[marks(ppp.marks)!=yes,]

  switch(p.type, "CSR"={
    simulated.sp<-rpoispp(lambda=sum(intensity(ppp.marks.sp)),
                          win=Window(ppp.marks.sp),nsim=1)
    marks(simulated.sp)<-factor(rep(yes,npoints(simulated.sp)))
  },"in.unif"={
    simulated.sp<-runifpoint(n=npoints(ppp.marks.sp),win=Window(ppp.marks.sp),
                             nsim=1)
    marks(simulated.sp)<-marks(ppp.marks.sp)
  },stop("Unsupported point type!")
  )

  simulated.species<-superimpose(simulated.sp,ppp.marks.nosp)

  return(simulated.species)
}

