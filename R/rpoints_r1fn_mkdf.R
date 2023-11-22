#'
#' @title Random simulation of the points corresponding to one level of a qualitative mark
#'
#' @description
#' Given a "ppp" object with a  mark dataframe, this function
#' does a random simulation of the points corresponding to one level of a qualitative mark,
#' randomly assigns corresponding mark rows,
#' and keeps the points corresponding to other levels unchanged.
#'
#' @usage rpoints.r1fn.mkdf(ppp.marks,p.type="CSR",V.species="Species",sp="Cork",...)
#'
#' @param ppp.marks A "\verb{ppp}" object with a mark dataframe including at least one qualitative variable
#' @param p.type Indicating what type of points to be generated, see \verb{Details}.
#' @param V.species The name of a qualitative mark variable
#' @param sp The level of which points to be simulated. If missing, the first level of \verb{V.species}
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
#' @seealso \code{\link{perm.r1r2.mkdf}}, \code{\link{perm.onemkvar}},
#' \code{\link{rpoints.r1fn.mkdf}}, \code{\link{rpoints.r1f2.1mk}}, \code{\link{simul.shift.r1f2.1mk}},
#'  \code{\link{simul.shift.r1fn.mkdf}}, \code{\link{simul.jitter.r1f2.1mk}}, \code{\link{simul.jitter.r1fn.mkdf}},
#'   \code{\link{heterotrl}}, \code{\link{contrl2}},
#' \code{\link{trig_g}},\code{\link{triK_K}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' rpoints.r1fn.mkdf(ppp.marks.c3,V.species="Species",sp="Cork")
#' }
#'
rpoints.r1fn.mkdf<-function(ppp.marks,p.type="CSR",V.species="Species",sp="Cork",...){

  stopifnot(is.ppp(ppp.marks) & is.data.frame(marks(ppp.marks)))
  stopifnot(V.species %in% names(marks(ppp.marks)))
  stopifnot(is.factor(marks(ppp.marks)[,V.species]))
  if(missing(sp)) sp<-levels(marks(ppp.marks)[V.species])[1]
  stopifnot(sp %in% marks(ppp.marks)[,V.species])

  ppp.marks.sp<-ppp.marks[marks(ppp.marks)[V.species]==sp,]
  ppp.marks.nosp<-ppp.marks[marks(ppp.marks)[V.species]!=sp,]

  switch(p.type, "CSR"={
    simulated.sp<-rpoispp(lambda=sum(intensity(ppp.marks.sp)),
                          win=Window(ppp.marks.sp),nsim=1)
    n.row<-sample.int(npoints(ppp.marks.sp),npoints(simulated.sp),replace=TRUE)
    marks(simulated.sp)<-marks(ppp.marks.sp)[n.row,]
  },"in.unif"={
    simulated.sp<-runifpoint(n=npoints(ppp.marks.sp),win=Window(ppp.marks.sp),
                             nsim=1)
    marks(simulated.sp)<-marks(ppp.marks.sp)
  },stop("Unsupported point type!")
  )

  simulate.species<-superimpose(simulated.sp,ppp.marks.nosp)
  return(simulate.species)
}

