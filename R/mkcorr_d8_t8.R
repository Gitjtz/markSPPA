#' @title Density Correlation Function, data 8
#'
#' @description
#' Estimates density correlation function for univariate pattern with one quantitative marks.
#'
#' @usage
#' mkcorr.d8.t8(X, m.qli="Infected", m.qti="DBH", ind=TRUE,
#'     r=NULL, delta=NULL, stoyan=0.15,
#'     method="density",
#'     weights=NULL,normalise=TRUE, internal=NULL, ...,
#'     mu.global=TRUE, dens.global=TRUE)
#'
#' @param X An object of class "\verb{ppp}"，with one quantitative mark.
#' @param m.qli A mark of data "\verb{X}" indicating the qualitative mark.
#' @param m.qti A mark of data "\verb{X}" indicating the quantitative mark.
#' @param ind logical, indicating the mark \verb{m.qti} of which level of \verb{m.qli} is used. If \verb{TRUE}, the mark corresponding to the second level of \verb{m.qli} is used;  else, the first level.
#' @param r Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param delta Same as in \verb{localpcf()} from the package \verb{spatstat}.
#' @param stoyan Same as in \verb{localpcf()} from the package \verb{spatstat}.
#' @param method Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param weights Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param normalise Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param internal Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param ... Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param mu.global Logical, indicating whether mean is global or local.  If TRUE, \eqn{\mu=\mu_{l+m}} ; else \eqn{\mu=\mu_l} or\eqn{\mu=\mu_m}
#' @param dens.global Logical, indicating whether desnity is global or local.  If TRUE, \eqn{K(r)=g_{l+m}(r)} ; else \eqn{g(r)=g_l(r)} or\eqn{g(r)=g_m(r)}
#'
#'
#' @details
#' According to the variables \verb{mu.global} and \verb{dens.global}, there are
#' four types of the test function (test 7) for data 8:
#' \deqn{C_{m,g}(r)=(m_{il}-\mu_l)(g_{jm}(r)-g_m(r))}
#' \deqn{C_{m,g}(r)=(m_{il}-\mu_l)(g_{jm}(r)-g_{l+m}(r))}
#' \deqn{C_{m,g}(r)=(m_{il}-\mu_{l+m})(g_{jm}(r)-g_m(r))}
#' \deqn{C_{m,g}(r)=(m_{il}-\mu_{l+m})(g_{jm}(r)-g_{l+m}(r))}
#' where the function \eqn{g_i()} is the local pair correlation function of
#' point \emph{i} and \eqn{g()} is the pair correlation function which equals
#' to the mean of \eqn{g_i()}, \eqn{m_{il}} is the mark value of point \emph{i}
#' of mark \emph{l}, and \eqn{\mu_l} is mean of mark \emph{l} of all points.
#'
#' This function is also suitable for data 9. In the case of data 9, both
#' \verb{mu.global} and \verb{dens.global} should be set to \verb{FALSE}.
#'
#' This function is time consuming.
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
#' @importFrom spatstat.explore localpcf localpcf
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
#' @seealso \code{\link{mkcorr.d8.t7}}, \code{\link{mkcorr.d7.t8}}, \code{\link{mkcorr.d8.t6}},
#' \code{\link{mkcorr.t8}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' nsim=29
#' ev.mkc.d8t8.c3<-envelope(ppp.ash.infect.dbh.c3, m.qli="Infected", m.qti="DBH",
#'     mkcorr.d8.t8,
#'     funargs=list(stoyan=0.15), nsim=nsim,
#'     simulate=NULL)
#' plot(ev.mkc.d8t8.c3)
#' }
#'
mkcorr.d8.t8<-function(X, m.qli="Infected", m.qti="DBH", ind=TRUE,
                   r=NULL, delta=NULL, stoyan=0.15,
                   method="density", weights=NULL,
                   normalise=TRUE, internal=NULL, ...,
                   mu.global=TRUE, dens.global=TRUE) {

  stopifnot(is.ppp(X) && is.marked(X))
  stopifnot(is.data.frame(marks(X)))
  stopifnot(all(c(m.qli, m.qti) %in% names(marks(X))))
  stopifnot(is.numeric(X$marks[[m.qti]]) && is.factor(marks(X)[[m.qli]]))
  stopifnot(all(complete.cases(X$marks)))
  if(nlevels(X$marks[[m.qli]])>2)
    warning("Multiple factor, only the first two are calculated!")

  f = function(m1, m2) {m1*m2}
  f1<-NULL
  h <- check.testfun(f, f1, X)
  f <- h$f
  f1 <- h$f1
  ftype <- h$ftype

  if(ftype!="general") stop("The fucntion input is not supported in current version!")

  if (ind){
    yes<-levels(marks(X)[1,m.qli])[2]
    no<-levels(marks(X)[1,m.qli])[1]
  } else {
    yes<-levels(marks(X)[1,m.qli])[1]
    no<-levels(marks(X)[1,m.qli])[2]
  }

  X.no<-subset(X,marks(X)[[m.qli]]==no)
  X.yes<-subset(X,marks(X)[[m.qli]]==yes)

  marx <- marks(X)[[m.qti]]
  marx.yes <- marks(X)[[m.qti]][ marks(X)[[m.qli]]==yes]
  npts <- npoints(X)
  npts.no<-npoints(X.no)


  if (unweighted <- is.null(weights)) {
    weights <- rep(1, npts)
  }  else {
    weights <- pointweights(X, weights = weights, parent = parent.frame())
    stopifnot(all(weights > 0))
  }

  weights.yes<-weights[marks(X)[[m.qli]]==yes]

  W <- X$window
  rmaxdefault <- rmax.rule("K", W, npts/area(W))
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  # correction.given <- !missing(correction)
  # correction <- pickoption("correction", correction, c(none = "none",
  #                                                      isotropic = "isotropic", Ripley = "isotropic", trans = "translate",
  #                                                      translate = "translate", translation = "translate", best = "best"),
  #                          multi = FALSE)
  # correction <- implemented.for.K(correction, W$type, correction.given)

  if(mu.global) {
    m.I<-marx.yes-mean(marx)
    sd.I<-sd(marx.yes)

  } else {
    m.I<-marx.yes-mean(marx.yes)
    sd.I<-sd(marx.yes)

  }

  if(dens.global){
    cross<-crosspairs(X.yes, X.no, rmax=rmax)
    dIJ<-cross$d
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
    m.J.m<-m.J.m[marks(X)[[m.qli]]==no,]
    sd.J<-sd(m.J.m)
    mI<-m.I[cross$i]
    mJ<-m.J.m[cross$j,]
    mJ<-diag(mJ)
    ff<-mI*mJ

  } else {
    cross<-crosspairs(X.yes, X.no, rmax=rmax)
    dIJ<-cross$d
    dIJ.fac<-factor(dIJ)
    dIJ.pos<-as.numeric(dIJ.fac)
    dIJ.lev<-levels(dIJ.fac)
    dIJ.uni<-as.numeric(dIJ.lev)
    m.J.localpcf<-matrix(NA,nrow=npts.no,ncol=length(dIJ.uni))
    for (i in 1: length(dIJ.uni)){
      m.J.localpcf[,i]<-localpcf(X.no, rmax=rmax, delta=delta,stoyan=stoyan,nr=length(r),
                             rvalue=dIJ.uni[i])
    }
    m.J.localpcf<-m.J.localpcf[,dIJ.pos]
    m.J.Kest<-apply(m.J.localpcf,2,mean)
    m.J.m<-m.J.localpcf-rep(m.J.Kest,each=nrow(m.J.localpcf))

    sd.J<-sd(m.J.m)
    mI<-m.I[cross$i]
    mJ<-m.J.m[cross$j,]
    mJ<-diag(mJ)
    ff<-mI*mJ

  }


  if (!unweighted) {
    sd.I.w<-sd(marx.yes*weights.yes)
    sd.J.w<-sd.J
    ff<-ff*weights.yes
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

  ylab <- quote(C[m1K2](r))
  yexp <- quote(C[m1K2](r))
  fnam <-c("C","m1K2")

  result <- fv(result, "r", ylab, "theo", , alim,
               c("r", "{%s[%s]^{iid}}(r)"), desc, fname = fnam)

  edgewt <- rep.int(1, length(dIJ))
  kmm<- sewsmod(dIJ, ff, edgewt, Efdenom, r, method,...)

  result <- bind.fv(result, data.frame(un = kmm), "{hat(%s)[%s]^{un}}(r)",
                    "uncorrected estimate of %s", "un")

  # switch(correction,none={
  #   result <- bind.fv(result, data.frame(un = kmm), "{hat(%s)[%s]^{un}}(r)",
  #                     "uncorrected estimate of %s", "un")
  # },translate={
  #   result <- bind.fv(result, data.frame(trans = kmm),
  #                     "{hat(%s)[%s]^{trans}}(r)", "translation-corrected estimate of %s",
  #                     "trans")
  # },isotropic={
  #   result <- bind.fv(result, data.frame(iso = kmm), "{hat(%s)[%s]^{iso}}(r)",
  #                     "Ripley isotropic correction estimate of %s",
  #                     "iso")
  # })

  nama2 <- names(result)
  corrxns <- rev(nama2[nama2 != "r"])
  formula(result) <- (. ~ r)
  fvnames(result, ".") <- corrxns
  unitname(result) <- unitname(X)
  attr(result,"yexp")<-yexp
  return(result)

}
