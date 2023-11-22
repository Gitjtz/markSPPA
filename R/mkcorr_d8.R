
#' Mark correlation functions for an univariate pattern with one quantitative mark and one qualitative mark
#'
#' @title \verb{markcorr} for data 8 tests 1-5
#'
#' @description
#' Estimates mark correlation functions to an univariate pattern with
#' one quantitative mark and one qualitative mark. The test functions include
#' mark correlation function, r-mark correlation function,mark variogram,
#' and Moran's I.
#'
#' @usage
#' mkcorr.d8(X, m.qli="Infected", m.qti="DBH",
#'     f = function(m1, m2) {m1*m2}, r=NULL,
#'     correction=c("isotropic", "Ripley", "translate"),
#'     method="density", ...,
#'     weights=NULL, f1=NULL, normalise=TRUE, fargs=NULL, internal=NULL)
#'
#' @param X An object of class "\verb{ppp}"，with one qualitative mark and one quantitative mark.
#' @param m.qli A mark of data "\verb{X}" indicates the qualitative mark.
#' @param m.qti A mark of data "\verb{X}" indicates the quantitative mark.
#' @param f Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param r Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param correction Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param method Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param ... Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param weights Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param f1 Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param normalise Same as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param fargs SSame as in \verb{markcorr()} from the package \verb{spatstat}.
#' @param internal Same as in \verb{markcorr()} from the package \verb{spatstat}.
#'
#' @details
#' An example for this data structure is a pattern of surviving and dead trees
#' of a given species (qualitative mark), where the size of the trees is also
#' known (the quantitative mark). The basic interest in analyzing these types
#' of data is to find out whether the quantitative mark (e.g., size)
#' shows a spatial correlation with the qualitatively marked points
#' (i.e., surviving and dead)(Wiegand and Moloney,2013).
#'
#' The data type and test function are referred to (Wiegand and Moloney, 2013).
#' Data 8 is an univariate pattern with one quantitative mark and one qualitative mark.
#'
#' The following test functions (tests 1-5) have been tested.
#' \itemize{
#'     \item f = function(m1,m2) {m1*m2}
#'     \item f = function(m1,m2) {m1}
#'     \item f = function(m1,m2) {m2}
#'     \item f = function(m1,m2) {0.5*(m1-m2)^2}
#'     \item f = function(m1,m2) {(m1-mean(c(m1,m2)))*(m2-mean(c(m1,m2)))}
#'     \item f = function(m1,m2) {(m1-mean(m1))*(m2-mean(m2))}
#' }
#' The last test can be a second version of test 5.
#' Other f(m1,m2) functions may also work.
#'
#' The function \code{\link{perm.onemkvar}} may be used to
#' permute one the two marks for envelope simulation.
#'
#' @return A object of "\verb{fv}".
#'
#' @importFrom spatstat.geom is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#'
#' @importFrom spatstat.geom pointweights area handle.r.b.args pickoption
#' @importFrom spatstat.explore check.testfun rmax.rule implemented.for.K fv
#'
#' @importFrom spatstat.geom closepairs ppp unitname crosspairs unitname<-
#' @importFrom spatstat.explore sewsmod bind.fv edge.Trans edge.Ripley formula.fv
#'
#' @importFrom spatstat.explore fvnames formula<- fvnames<-
#'
#' @importFrom stats complete.cases sd weighted.mean var
#'
#' @export
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Raventós, J., Wiegand, T., and de Luis, M. 2010. Evidence for the spatial
#' segregation hypothesis: a test with nine-year survivorship data in a
#' Mediterranean shrubland. Ecology. 91 7:2110-2120.
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' @seealso \code{\link{mkcorr.d8.t6}}, \code{\link{mkcorr.d7}},
#' \code{\link{mkcorr.d9}}, \code{\link{mkcorr.t7}} ,\code{\link{mkcorr.t8}}
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' nsim=2499
#' ev.mkc.m1m2.ash.infect.dbh.c3<-envelope(ppp.ash.infect.dbh.c3,mkcorr.d8,
#'     funargs=list(m.qli="Infected",m.qti="DBH",
#'     f = function(m1, m2) {m1 * m2}), nsim=nsim,
#'     simulate=expression(perm.onemkvar(ppp.ash.infect.dbh.c3,V.simul="DBH")))
#' plot(ev.mkc.m1m2.ash.infect.dbh.c3)
#'
#' ev.mkc.MI.ash.infect.dbh.c3<-envelope(ppp.ash.infect.dbh.c3,mkcorr.d8,
#'     funargs=list(m.qli="Infected",m.qti="DBH",
#'     f = function(m1, m2) {(m1-mean(c(m1,m2)))*(m2-mean(c(m1,m2)))}), nsim=nsim,
#'     simulate=expression(perm.onemkvar(ppp.ash.infect.dbh.c3,V.simul="DBH")))
#' plot(ev.mkc.MI.ash.infect.dbh.c3)
#'
#' ev.mkc.vario.ash.infect.dbh.c3<-envelope(ppp.ash.infect.dbh.c3,mkcorr.d8,
#'     funargs=list(m.qli="Infected",m.qti="DBH",
#'     f=function(m1,m2){0.5*(m1-m2)^2}), nsim=nsim,
#'     imulate=expression(perm.onemkvar(ppp.ash.infect.dbh.c3,V.simul="DBH")))
#' plot(ev.mkc.vario.ash.infect.dbh.c3)
#' }
#'
mkcorr.d8<-function(X, m.qli="Infected", m.qti="DBH",
                   f = function(m1, m2) {m1*m2}, r=NULL,
                   correction=c("isotropic", "Ripley", "translate"),
                   method="density", ..., weights=NULL,
                   f1=NULL, normalise=TRUE, fargs=NULL, internal=NULL){

  stopifnot(is.ppp(X) && is.marked(X))
  stopifnot(all(c(m.qli, m.qti ) %in% names(marks(X))))
  stopifnot(!I(is.null(f) && is.null(f1)))
  X<-subset(X,select=c(m.qli,m.qti))
  stopifnot(is.factor(X$marks[[m.qli]]))
  stopifnot(all(complete.cases(X$marks[[m.qli]])) && all(complete.cases(X$marks[[m.qti]])))
  if(nlevels(X$marks[[m.qli]])>2)
    warning("Multiple factor, only the first two are calculated!")

  h <- check.testfun(f, f1, X)
  f <- h$f
  f1 <- h$f1
  ftype <- h$ftype

  if(ftype!="general") stop("The fucntion input is not supported in current version!")

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
  m.I<-marx[[m.qti]][as.logical(marx[[1]])]
  m.J<-marx[[m.qti]][as.logical(marx[[2]])]

  W <- X$window
  rmaxdefault <- rmax.rule("K", W, npts/area(W))
  breaks <- handle.r.b.args(r, NULL, W, rmaxdefault = rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  # correction.given<-TRUE
  correction.given <- !missing(correction) && !is.null(correction)
  if (is.null(correction))
    correction <- c("isotropic", "Ripley", "translate")

  correction <- pickoption("correction", correction,
                           c(none = "none", border = "border", bord.modif = "bord.modif",
                             isotropic = "isotropic", Ripley = "isotropic",
                             translate = "translate", translation = "translate",
                             best = "best"), multi = TRUE)
  correction <- implemented.for.K(correction, W$type, correction.given)

  m12<-gsub(" ","",deparse(f)[3])

  Ef <- internal$Ef
  if (is.null(Ef)) {
    Ef <- switch(ftype, mul = {
      stop("Error: invalid function type!")
    }, equ = {
      stop("Error: invalid function type!")
    },product = {
      stop("Error: invalid function type!")
    }, general = {
      if (is.null(fargs)) {

        if (unweighted) {
          switch(m12, "m1*m2"={
            mean(marx[[m.qti]])^2
          }, "m1"={
            mean(marx[[m.qti]])
          }, "m2"={
            mean(marx[[m.qti]])
          }, "0.5*(m1-m2)^2"={
            var(marx[[m.qti]])
          }, "(m1-mean(m1))*(m2-mean(m2))"={
            var(marx[[m.qti]])
          }, "(m1-mean(c(m1,m2)))*(m2-mean(c(m1,m2)))"={
            var(marx[[m.qti]])
          }, mean(outer(marx[[m.qti]],marx[[m.qti]],f))
          )#switch end
        } else {
          switch(m12, "m1*m2"={
            weighted.mean(marx[[m.qti]], weights)^2
          }, "m1"={
            weighted.mean(marx[[m.qti]], weights)
          }, "m2"={
            weighted.mean(marx[[m.qti]], weights)
          }, "0.5*(m1-m2)^2"={
            var(marx[[m.qti]]*weights)
          }, "(m1-mean(m1))*(m2-mean(m2))"={
            var(marx[[m.qti]]*weights)
          }, "(m1-mean(c(m1,m2)))*(m2-mean(c(m1,m2)))"={
            var(marx[[m.qti]]*weights)
          }, mean(outer(marx[[m.qti]]*weights,marx[[m.qti]]*weights,f))
          )#switch end

        }

      } else {
        do.call(outer, append(list(marx[[m.qti]], marx[[m.qti]], f), fargs))
      }

    }, stop("Internal error: invalid ftype"))#switch end
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

  if (ftype %in% c("mul", "equ")) {
    ylab <- substitute(k[mm](r), list(mm = paste(levels(mari),collapse=","))) #quote(k[mm](r))
    yexp <- substitute(k[mm](r), list(mm = paste(levels(mari),collapse=",")))
    fnam <-c("k",paste0("list(",levels(mari)[1],",",levels(mari)[2],")")) #c("k", "mm")

  }    else {

    dot<-'~symbol(\"\\267\")'

    ylab<-switch(m12,"m1*m2"={
      substitute(k[mm](r), list(mm = paste(levels(mari),collapse=",")))
    },"m1"={
      substitute(k[m,d](r),list(m=levels(mari)[1],d=dot))
    },"m2"={
      substitute(k[m,d](r),list(m=levels(mari)[2],d=dot))
    },"0.5*(m1-m2)^2"={
      substitute(gamma[mm](r), list(mm = paste(levels(mari),collapse=",")))
    },"(m1-mean(m1))*(m2-mean(m2))"={
      substitute(I[mm](r), list(mm = paste(levels(mari),collapse=",")))
    },"(m1-mean(c(m1,m2)))*(m2-mean(c(m1,m2)))"={
      substitute(I[mm](r), list(mm = paste(levels(mari),collapse=",")))
    },quote(k[f](r))
    )

    yexp<-switch(m12,"m1*m2"={
      substitute(k[mm](r), list(mm = paste(levels(mari),collapse=",")))
    },"m1"={
      substitute(k[m~symbol("\267")](r),list(m=levels(mari)[1]))
    },"m2"={
      substitute(k[m~symbol("\267")](r),list(m=levels(mari)[2]))
    },"0.5*(m1-m2)^2"={
      substitute(gamma[mm](r), list(mm = paste(levels(mari),collapse=",")))
    },"(m1-mean(m1))*(m2-mean(m2))"={
      substitute(I[mm](r), list(mm = paste(levels(mari),collapse=",")))
    },"(m1-mean(c(m1,m2)))*(m2-mean(c(m1,m2)))"={
      substitute(I[mm](r), list(mm = paste(levels(mari),collapse=",")))
    },quote(k[f](r))
    )

    fnam.2<-switch(m12,"m1*m2"={
      substitute(mm, list(mm = paste(levels(mari),collapse=",")))
    },"m1"={
      substitute(m,list(m=levels(mari)[1]))
    },"m2"={
      substitute(m,list(m=levels(mari)[2]))
    },"0.5*(m1-m2)^2"={
      substitute(mm, list(mm = paste(levels(mari),collapse=",")))
    },"(m1-mean(m1))*(m2-mean(m2))"={
      substitute(mm, list(mm = paste(levels(mari),collapse=",")))
    },"(m1-mean(c(m1,m2)))*(m2-mean(c(m1,m2)))"={
      substitute(mm, list(mm = paste(levels(mari),collapse=",")))
    }
    )

    fnam<-switch(m12,"m1*m2"={
      c("k",paste0("list(",fnam.2,")"))
    },"m1"={
      c("k", paste0("list(",fnam.2,",",dot,")"))
    },"m2"={
      c("k", paste0("list(",fnam.2,",",dot,")"))
    },"0.5*(m1-m2)^2"={
      c("gamma",paste0("list(",fnam.2,")"))
    },"(m1-mean(m1))*(m2-mean(m2))"={
      c("I",paste0("list(",fnam.2,")"))
    },"(m1-mean(c(m1,m2)))*(m2-mean(c(m1,m2)))"={
      c("I",paste0("list(",fnam.2,")"))
    },c("k", "f")
    )

  }


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
