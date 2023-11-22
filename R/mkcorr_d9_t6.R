#'
#' @title \verb{markcorr} for data 9 test 6
#'
#' @description
#' Estimating Schlather I mark correlation function for a bivariate pattern
#' with one quantitative mark.
#'
#' @usage
#' mkcorr.d9.t6(X, m.qli="Species", m.qti="Disease",r=NULL,
#'     correction=c("isotropic", "Ripley", "translate"),
#'     method="density", ...,
#'     weights=NULL, normalise=TRUE, internal=NULL)
#'
#' @param X An object of class "\verb{ppp}"，with one qualitative mark and one quantitative mark.
#' @param m.qli A mark of data "\verb{X}" indicating the qualitative mark.
#' @param m.qti A mark of data "\verb{X}" indicating the quantitative mark.
#' @param r Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param correction Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param method Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param ... Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param weights Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param normalise Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param internal Same as in \verb{markcorr()} from the package \verb{spatstat}.
#'
#' @details
#' The data type and test function are referred to (Wiegand and Moloney, 2013).
#' Data 9 is a bivariate pattern with one quantitative mark.
#'
#' This function tests (test 6)
#' \deqn{(m_{il}-\mu_l(r))(m_{jm}-\mu_m(r))}
#' The result will be normalised by sd1*sd2, if \verb{normalised==TRUE}.
#'
#' In the case of data 9, a bivariate pattern with one quantitative mark,
#' the marks must only be randomized within each pattern (i.e., a point of
#' pattern 1 should not receive a mark of a point of pattern 2 and vice versa)
#' (Wiegand and Moloney,2013). So employ the function \code{\link{perm.r1r2.mkdf}} to do
#' envelope simulation.
#'
#' @return A object of "\verb{fv}".
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#'
#' @importFrom spatstat.geom pointweights area handle.r.b.args pickoption
#' @importFrom spatstat.explore check.testfun rmax.rule implemented.for.K fv
#'
#' @importFrom spatstat.geom closepairs ppp unitname crosspairs
#' @importFrom spatstat.explore sewsmod bind.fv edge.Trans edge.Ripley formula.fv
#'
#' @importFrom stats complete.cases sd weighted.mean var
#'
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' @seealso \code{\link{mkcorr.d9}}, \code{\link{mkcorr.d7.t6}},
#' \code{\link{mkcorr.d8.t6}}, \code{\link{mkcorr.t7}} ,\code{\link{mkcorr.t8}}
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' nsim=2499
#' ev.mkc.t6.ash.birch.disease.c4<-envelope(ppp.marks.ash.birch.c4,mkcorr.d9.t6,
#'     funargs=list(m.qli="Species",m.qti="Disease"), nsim=nsim,
#'     simulate=expression(perm.r1r2.mkdf(ppp.marks.ash.birch.c4,
#'     m.qli="Species",m.qti="Disease")))
#' plot(ev.mkc.t6.ash.birch.disease.c4)
#'
#' ev.mkc.t6.c4<-anylapply(ppp.marks.sp4c.c4,envelope,mkcorr.d9.t6,
#'     funargs=list(m.qli="Species",m.qti="Disease"), nsim=nsim,
#'     simulate=simul.bivar.sep)
#' plot(ev.mkc.t6.c4, main="", main.panel=letters[1:6], legend=FALSE)
#' }
#'
mkcorr.d9.t6<-function(X, m.qli="Species", m.qti="Disease", r=NULL,
                       correction=c("isotropic", "Ripley", "translate"),
                       method="density", ..., weights=NULL,
                       normalise=TRUE, internal=NULL){

  stopifnot(is.ppp(X) && is.marked(X))
  stopifnot(all(c(m.qli, m.qti ) %in% names(marks(X))))

  X<-subset(X,select=c(m.qli,m.qti))
  stopifnot(is.factor(X$marks[[1]]) && is.numeric(X$marks[[2]]))
  stopifnot(all(complete.cases(X$marks[[1]])) && all(complete.cases(X$marks[[2]])))
  if(nlevels(X$marks[[1]])>2)
    warning("Multiple factor, only the first two are calculated!")

  f = function(m1, m2) {m1 * m2}
  f1<-NULL
  h <- check.testfun(f, f1, X)
  f <- h$f
  f1 <- h$f1
  ftype <- h$ftype

  marx <- marks(X, dfok = TRUE)
  mari <- marx[[m.qli]]
  levi <- levels(mari)
  dumi <- 1 * outer(mari, levi, "==")
  colnames(dumi) <- paste(m.qli, levi, sep=".")
  marx <- as.data.frame(append(marx[, m.qti, drop = FALSE],
                               list(dumi), after = 1 - 1))
  npts <- npoints(X)

  if (unweighted <- is.null(weights)) {
    weights <- rep(1, npts)
  }  else {
    weights <- pointweights(X, weights = weights, parent = parent.frame())
    stopifnot(all(weights > 0))
  }

  w.I<-weights[as.logical(marx[[1]])]
  w.J<-weights[as.logical(marx[[2]])]
  m.I<-marx[[ncol(marx)]][as.logical(marx[[1]])]
  m.J<-marx[[ncol(marx)]][as.logical(marx[[2]])]

  W <- X$window
  rmaxdefault <- rmax.rule("K", W, npts/area(W))
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  correction.given <- !missing(correction) && !is.null(correction)
  if (is.null(correction))
    correction <- c("isotropic", "Ripley", "translate")

  correction <- pickoption("correction", correction,
                           c(none = "none", border = "border", bord.modif = "bord.modif",
                             isotropic = "isotropic", Ripley = "isotropic",
                             translate = "translate", translation = "translate",
                             best = "best"), multi = TRUE)
  correction <- implemented.for.K(correction, W$type, correction.given)

  Ef <- internal$Ef
  if (is.null(Ef)) {
    Ef <- if (unweighted){
      sd(m.I)*sd(m.J)
    } else {
      sd(m.I*w.I)*sd(m.J*w.J)
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
    # else if (Efdenom < 0)
    #   warning(paste("Problem when normalising the mark correlation:",
    #                 "the denominator is negative"))
  }


  result <- data.frame(r=r, theo= rep.int(theory, length(r)))
  desc <- c("distance argument r", "theoretical value (independent marks) for %s")
  alim <- c(0, min(rmax, rmaxdefault))

  ylab <- substitute(I[mm](r), list(mm = paste(levels(mari),collapse=",")))
  yexp <- substitute(I[mm](r), list(mm = paste(levels(mari),collapse=",")))
  fnam <-c("I",paste0("list(",levels(mari)[1],",",levels(mari)[2],")"))


  result <- fv(result, "r", ylab, "theo", , alim,
               c("r", "{%s[%s]^{iid}}(r)"), desc, fname = fnam)

  X.I<-subset(X,marks(X)[[m.qli]]==levels(marks(X)[[m.qli]])[1])
  X.J<-subset(X,marks(X)[[m.qli]]==levels(marks(X)[[m.qli]])[2])
  close<-crosspairs(X.I,X.J,rmax)
  dIJ <- close$d
  I <- close$i
  J <- close$j

  mI<-m.I[I]
  mJ<-m.J[J]

  mI.t6<-mI
  mJ.t6<-mJ
  for(i in 2:length(r)){
    mI[which(dIJ>r[i-1]) &  which(dIJ<=r[i])]<-mI.t6[which(dIJ>r[i-1]) &  which(dIJ<=r[i])]-
      mean(mI.t6[which(dIJ>r[i-1]) &  which(dIJ<=r[i])])
    mJ[which(dIJ>r[i-1]) &  which(dIJ<=r[i])]<-mJ.t6[which(dIJ>r[i-1]) &  which(dIJ<=r[i])]-
      mean(mJ.t6[which(dIJ>r[i-1]) &  which(dIJ<=r[i])])
  }

  fargs<-NULL
  XI <- ppp(close$xi, close$yi, window = W, check = FALSE)
  ff <- switch(ftype, mul = mI * mJ, equ = (mI == mJ), product = {
    if (is.null(fargs)) {
      fI <- f1(mI)
      fJ <- f1(mJ)
    } else {
      fI <- do.call(f1, append(list(mI), fargs))
      fJ <- do.call(f1, append(list(mJ), fargs))
    }
    fI * fJ
  }, general = {
    if (is.null(fargs)) f(mI, mJ) else do.call(f,
                                               append(list(mI, mJ), fargs))
  })


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
  if (!unweighted)
    ff <- ff * weights[I] * weights[J]
  if (any(correction == "none")) {
    edgewt <- rep.int(1, length(dIJ))
    Mnone <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method,
                     ...)
    result <- bind.fv(result, data.frame(un = Mnone), "{hat(%s)[%s]^{un}}(r)",
                      "uncorrected estimate of %s", "un")
  }
  if (any(correction == "translate")) {
    XJ <- ppp(close$xj, close$yj, window = W, check = FALSE)
    edgewt <- edge.Trans(XI, XJ, paired = TRUE)
    Mtrans <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method,...)
    result <- bind.fv(result, data.frame(trans = Mtrans),
                      "{hat(%s)[%s]^{trans}}(r)", "translation-corrected estimate of %s",
                      "trans")
  }
  if (any(correction == "isotropic")) {
    edgewt <- edge.Ripley(XI, matrix(dIJ, ncol = 1))
    Miso <- sewsmod(dIJ, ff, edgewt, Efdenom, r, method,...)
    result <- bind.fv(result, data.frame(iso = Miso), "{hat(%s)[%s]^{iso}}(r)",
                      "Ripley isotropic correction estimate of %s",
                      "iso")
  }
  nama2 <- names(result)
  corrxns <- rev(nama2[nama2 != "r"])
  formula(result) <- (. ~ r)
  fvnames(result, ".") <- corrxns
  unitname(result) <- unitname(X)
  attr(result,"yexp")<-yexp
  return(result)

}


