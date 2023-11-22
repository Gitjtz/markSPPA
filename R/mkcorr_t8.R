#' @title Density Correlation Function
#'
#' @description
#' Estimates density correlation function for univariate pattern with one quantitative marks.
#'
#' @usage
#' mkcorr.t8(X, r=NULL, delta=NULL, stoyan=0.15, method="density",
#'     weights=NULL,normalise=TRUE, internal=NULL,...)
#'
#' @param X An object of class "\verb{ppp}"，with one quantitative mark.
#' @param r Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param delta Same as in \verb{localpcf()} from the package \verb{spatstat}.
#' @param stoyan Same as in \verb{localpcf()} from the package \verb{spatstat}.
#' @param method Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param weights Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param normalise Same as in \verb{markcorr()} from the package \verb{spatstat}
#' @param internal Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param ... Same as in \verb{markcorr()} from the package \verb{spatstat}.
#'
#' @details
#' The test function (test 8), is
#' \deqn{C_{m,g}(r)=(m_i-\mu)(g_i(r)-g(r))}
#' where the function \eqn{g_i()} is the local pair-correlation function of each
#' point \emph{i} and \eqn{g()} is the pair-correlation function which equals
#' to the mean of \eqn{g_i()}, \eqn{m_i} is the mark value of point \emph{i}
#' and \eqn{\mu} is the mean mark of all points.
#'
#' For envelope simulation, user may do a permutation of the mark or a simulation
#' of a CSR of all points.
#'
#' This function is time consuming.
#'
#' This function is similar to the function \eqn{C_{m,K}(r)} which is a cumulative one.
#'
#' @return A object of "\verb{fv}".
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#'
#' @importFrom spatstat.geom pointweights area handle.r.b.args pickoption
#' @importFrom spatstat.explore check.testfun rmax.rule implemented.for.K fv
#'
#' @importFrom spatstat.geom closepairs ppp unitname
#' @importFrom spatstat.explore sewsmod bind.fv edge.Trans edge.Ripley formula.fv
#'
#' @importFrom spatstat.explore localK localpcf
#'
#' @importFrom stats complete.cases sd weighted.mean var
#'
#' @export
#' @references
#' Fedriani, J, T. Wiegand, G. Calvo, A. Suarez-Esteban, M. Jacome, M. Zywiec,
#' and M.Delibes. 2015. Unravelling conflicting density- and distance-dependent
#' effects on plant reproduction using a spatially explicit approach.
#' Journal of Ecology: 103: 1344-1353
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' @seealso \code{\link{mkcorr.t7}}, \code{\link{mkcorr.d9.t6}}, \code{\link{mkcorr.d7.t6}},
#' \code{\link{mkcorr.d8.t6}}, \code{\link{mkcorr.d7.t8}}, \code{\link{mkcorr.d8.t8}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' the data spruces from package spatstat
#' nsim=2499
#' mkc.t8.spruces<-mkcorr.t8(spruces)
#' plot(mkc.t8.spruces)
#' ev.mkc.t8.spruces<-envelope(spruces,mkcorr.t8, nsim=nsim,simulate=rlabel)
#' plot(ev.mkc.t8.spruces)
#' ev.mkc.t8.spruces<-envelope(spruces,mkcorr.t8, nsim=nsim)
#' }
#'
mkcorr.t8<-function(X, r=NULL, delta=NULL, stoyan=0.15,
                    method="density", weights=NULL,
                    normalise=TRUE, internal=NULL,...) {

  stopifnot(is.ppp(X) && is.marked(X))
  stopifnot(is.numeric(X$marks))
  stopifnot(all(complete.cases(X$marks)))

  f = function(m1, m2) {m1*m2}
  f1<-NULL
  h <- check.testfun(f, f1, X)
  f <- h$f
  f1 <- h$f1
  ftype <- h$ftype

  if(ftype!="general") stop("The fucntion input is not supported in current version!")

  marx <- marks(X)
  npts <- npoints(X)

  if (unweighted <- is.null(weights)) {
    weights <- rep(1, npts)
  }  else {
    weights <- pointweights(X, weights = weights, parent = parent.frame())
    stopifnot(all(weights > 0))
  }

  W <- X$window
  rmaxdefault <- rmax.rule("K", W, npts/area(W))
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  m.I<-marx-mean(marx)
  sd.I<-sd(marx)

  close<-closepairs(X, rmax=rmax)
  dIJ<-close$d
  dIJ.fac<-factor(dIJ)
  dIJ.pos<-as.numeric(dIJ.fac)
  dIJ.lev<-levels(dIJ.fac)
  dIJ.uni<-as.numeric(dIJ.lev)
  m.J.localpcf<-matrix(NA,nrow=npts,ncol=length(dIJ.uni))
  for (i in 1: length(dIJ.uni)){
    m.J.localpcf[,i]<-localpcf(X, rmax=rmax, delta=delta,stoyan=stoyan,nr=length(r),
                               rvalue=dIJ.uni[i])
  }
  m.J.localpcf<-m.J.localpcf[,dIJ.pos]
  m.J.Kest<-apply(m.J.localpcf,2,mean)
  m.J.m<-m.J.localpcf-rep(m.J.Kest,each=nrow(m.J.localpcf))

  sd.J<-sd(m.J.m)
  mI<-m.I[close$i]
  mJ<-m.J.m[close$j,]
  mJ<-diag(mJ)
  ff<-mI*mJ

  if (!unweighted) {
    ff <- ff * weights
    sd.I.w<-sd(marx*weights)
    sd.J.w<-sd.J
  }


  if (is.logical(ff))
    ff <- as.numeric(ff)
  if (!is.numeric(ff))
    stop("function f did not return numeric values")

  if (anyNA(ff))
    switch(ftype, mul = , equ = stop("some marks were NA"),
           product = , general = stop("function f returned some NA values"))
  # if (any(ff < 0))
  #   switch(ftype, mul = , equ = stop("negative marks are not permitted"),
  #          product = , general = stop("negative values of function f are not permitted"))

  Ef <- internal$Ef
  if (is.null(Ef)) {
    Ef <- if (unweighted){
      sd.I*sd.J
    } else {
      sd.I.w*sd.J.w
    }
  }


  if (normalise) {
    theory <- 1
    Efdenom <- Ef
  }    else {
    theory <- Ef
    Efdenom <- 1
  }

  if (normalise) {
    if (Efdenom == 0)
      stop("Cannot normalise the mark correlation; the denominator is zero")
    else if (Efdenom < 0)
      warning(paste("Problem when normalising the mark correlation:",
                    "the denominator is negative"))
  }

  result <- data.frame(r=r, theo= rep.int(theory, length(r)))
  desc <- c("distance argument r", "theoretical value (independent marks) for %s")
  alim <- c(0, min(rmax, rmaxdefault))

  ylab <- quote(C[mg](r))
  yexp <- quote(C[mg](r))
  fnam <-c("C","mg")

  result <- fv(result, "r", ylab, "theo", , alim,
               c("r", "{%s[%s]^{iid}}(r)"), desc, fname = fnam)


  edgewt <- rep.int(1, length(dIJ))
  kmm<- sewsmod(dIJ, ff, edgewt, Efdenom, r, method,...)
  result <- bind.fv(result, data.frame(un = kmm), "{hat(%s)[%s]^{un}}(r)",
                    "uncorrected estimate of %s", "un")

  nama2 <- names(result)
  corrxns <- rev(nama2[nama2 != "r"])
  formula(result) <- (. ~ r)
  fvnames(result, ".") <- corrxns
  unitname(result) <- unitname(X)
  attr(result,"yexp")<-yexp
  return(result)

}
