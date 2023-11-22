#' @title \emph{k}th Neighbour Distance Function \eqn{ND^k(r)}
#'
#' @description
#' Estimates the \emph{k}th nearest neighbour distance distribution function
#' \eqn{ND^k(r)} from a point pattern in a window of arbitrary shape.
#'
#' @usage
#' NDest(X, r = NULL, breaks = NULL, ...,
#'    k=1, correction = c("rs", "km", "han"), domain = NULL)
#'
#' @param X A object of class "\verb{ppp}", same as that in the funciton \verb{Gest} from the package \verb{spatstat.}
#' @param r same as that in the funciton \verb{Gest} from the package \verb{spatstat.}
#' @param breaks same as that in the funciton \verb{Gest} from the package \verb{spatstat.}
#' @param ... same as that in the funciton \verb{Gest} from the package \verb{spatstat.}
#' @param k Integer, indicating the \emph{k}th nearest neighbour.
#' @param correction same as that in the funciton \verb{Gest} from the package \verb{spatstat.}
#' @param domain same as that in the funciton \verb{Gest} from the package \verb{spatstat.}
#'
#' @details
#' This function is the same as \verb{Gest} in \verb{spatatat} when \verb{k=1}.
#'
#' @return A object of "\verb{fv}".
#'
#' @importFrom spatstat.geom verifyclass is.multitype marks npoints superimpose intensity is.ppp is.marked
#' @importFrom spatstat.explore pcfcross pcfdot Kcross Kdot eval.fv
#'
#' @importFrom spatstat.geom pointweights area handle.r.b.args pickoption
#' @importFrom spatstat.explore check.testfun rmax.rule implemented.for.K fv
#'
#' @importFrom spatstat.geom closepairs ppp unitname crosspairs ppsubset
#' @importFrom spatstat.explore sewsmod bind.fv edge.Trans edge.Ripley formula.fv
#'
#' @importFrom spatstat.geom nndist nncross bdist.points eroded.areas inside.owin is.subset.owin
#' @importFrom spatstat.explore rmax.rule makefvlabel km.rs km.rs.opt
#'
#' @importFrom graphics hist
#' @export
#'
#' @references
#' R package \verb{spatstat}
#'
#' @seealso \code{\link{NDcross}}, \code{\link{NDmulti}},
#' \code{\link{NDdot}},
#' \code{\link{triNDcross}}.
#'
#' @author Tianzhong Jing <Jingtianzhong@163.com>
#'
#' @references
#' Wiegand, T. and Moloney, K. A. 2013. Handbook of spatial point-pattern
#' analysis in ecology. Chapman and Hall/CRC, Boca Raton.
#'
#' Wiegand, T. 2018. User Manual for the Programita software.
#'
#' @seealso \code{\link{NDcross}}, \code{\link{NDmulti}},
#' \code{\link{NDdot}},
#' \code{\link{triNDcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' # data amacrine is from package spatstat
#' plot(NDest(amacrine,k=2))
#' }
#'
NDest<-function (X, r = NULL, breaks = NULL, ...,
                 k=1, correction = c("rs", "km", "han"), domain = NULL)
{
  verifyclass(X, "ppp")
  if (!is.null(domain))
    stopifnot(is.subset.owin(domain, Window(X)))
  W <- X$window
  npts <- npoints(X)
  lambda <- npts/area(W)
  rmaxdefault <- rmax.rule("G", W, lambda)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault = rmaxdefault)
  rvals <- breaks$r
  rmax <- breaks$max
  zeroes <- numeric(length(rvals))
  if (is.null(correction)) {
    correction <- c("rs", "km", "han")
  }
  else correction <- pickoption("correction", correction, c(none = "none",
                                                            border = "rs", rs = "rs", KM = "km", km = "km", Kaplan = "km",
                                                            han = "han", Hanisch = "han", cs = "han", ChiuStoyan = "han",
                                                            best = "km"), multi = TRUE)
  nnd <- nndist(X$x, X$y,k=k)
  bdry <- bdist.points(X)
  if (!is.null(domain)) {
    ok <- inside.owin(X, w = domain)
    nnd <- nnd[ok]
    bdry <- bdry[ok]
  }
  o <- pmin.int(nnd, bdry)
  d <- (nnd <= bdry)
  df <- data.frame(r = rvals, theo = 1 - exp(-lambda * pi *
                                               rvals^2))
  Z <- fv(df, "r", substitute({ND^k}(r), list(k=k)), "theo", . ~ r, c(0,
                                                            rmax), c("r", "%s[pois](r)"), c("distance argument r",
                                                                                            "theoretical Poisson %s"), fname = "ND")
  if ("none" %in% correction) {
    if (npts <= 1)
      edf <- zeroes
    else {
      hh <- hist(nnd[nnd <= rmax], breaks = breaks$val,
                 plot = FALSE)$counts
      edf <- cumsum(hh)/length(nnd)
    }
    Z <- bind.fv(Z, data.frame(raw = edf), "hat(%s)[raw](r)",
                 "uncorrected estimate of %s", "raw")
  }
  if ("han" %in% correction) {
    if (npts <= 1)
      G <- zeroes
    else {
      x <- nnd[d]
      a <- eroded.areas(W, rvals, subset = domain)
      h <- hist(x[x <= rmax], breaks = breaks$val, plot = FALSE)$counts
      G <- cumsum(h/a)
      G <- G/max(G[is.finite(G)])
    }
    Z <- bind.fv(Z, data.frame(han = G), "hat(%s)[han](r)",
                 "Hanisch estimate of %s", "han")
    attr(Z, "alim") <- range(rvals[G <= 0.9])
  }
  if (any(correction %in% c("rs", "km"))) {
    want.rs <- "rs" %in% correction
    want.km <- "km" %in% correction
    if (npts == 0) {
      result <- list(rs = zeroes, km = zeroes, hazard = zeroes,
                     theohaz = zeroes)
    }
    else {
      result <- km.rs.opt(o, bdry, d, breaks, KM = want.km,
                          RS = want.rs)
      if (want.km)
        result$theohaz <- 2 * pi * lambda * rvals
    }
    wanted <- c(want.rs, rep(want.km, 3L))
    wantednames <- c("rs", "km", "hazard", "theohaz")[wanted]
    result <- as.data.frame(result[wantednames])
    Z <- bind.fv(Z, result, c("hat(%s)[bord](r)", "hat(%s)[km](r)",
                              "hat(h)[km](r)", "h[pois](r)")[wanted], c("border corrected estimate of %s",
                                                                        "Kaplan-Meier estimate of %s", "Kaplan-Meier estimate of hazard function h(r)",
                                                                        "theoretical Poisson hazard function h(r)")[wanted],
                 if (want.km)
                   "km"
                 else "rs")
    attr(Z, "alim") <- with(Z, range(.x[.y <= 0.9]))
  }
  nama <- names(Z)
  fvnames(Z, ".") <- rev(setdiff(nama, c("r", "hazard", "theohaz")))
  unitname(Z) <- unitname(X)
  return(Z)
}
