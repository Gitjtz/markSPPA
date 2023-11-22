
#' @title Marked \emph{k}th Nearest Neighbour Distance Function
#'
#' @description
#' For a marked point pattern, estimate the distribution of the distance from a typical point in subset \emph{I}
#' to the \emph{k}th nearest point of subset \emph{J}
#'
#' @usage
#' NDmulti(X, I, J, r = NULL, breaks = NULL, ...,
#'     k=1, disjoint = NULL,
#'     correction = c("rs", "km", "han"))
#'
#' @param X A "\verb{ppp}" object. Same as in \verb{Gmulti} from the package \verb{spatstat}
#' @param I Subset of points of X from which distances are measured. Same as in \verb{Gmulti} from the package \verb{spatstat}
#' @param J Subset of points in X to which distances are measured. Same as in \verb{Gmulti} from the package \verb{spatstat}
#' @param r Same as in \verb{Gmulti} from the package \verb{spatstat}
#' @param breaks Same as in \verb{Gmulti} from the package \verb{spatstat}
#' @param ... Same as in \verb{Gmulti} from the package \verb{spatstat}
#' @param k Integer, indicating the \emph{k}th nearest neighbour.
#' @param disjoint Same as in \verb{Gmulti} from the package \verb{spatstat}
#' @param correction Same as in \verb{Gmulti} from the package \verb{spatstat}
#'
#' @details
#' This function is the same as \verb{Gmulti} in \verb{spatatat} when \verb{k=1}.
#'
#' @return An object of class "\verb{fv}"
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
#' @importFrom spatstat.geom nncross bdist.points eroded.areas
#' @importFrom spatstat.explore rmax.rule makefvlabel km.rs
#'
#' @export
#' @references
#' R package \verb{spatstat}
#'
#'
#' @seealso \code{\link{NDcross}}, \code{\link{ND_NDcross}},
#' \code{\link{ND_NDcross2}},\code{\link{ND_NDdot}},\code{\link{NDoverNDdot}},
#' \code{\link{triND_ND}},\code{\link{triNDcross}}.
#'
#' @author Tianzhong Jing \email{Jingtianzhong@163.com}
#'
#' @examples
#' \dontrun{
#' # data amacrine is from package spatstat
#' plot(NDmulti(amacrine, marks(amacrine)=="off", marks(amacrine)=="on"))
#' plot(NDmulti(amacrine, marks(amacrine)=="off", marks(amacrine)=="on",k=2))
#' }
#'
NDmulti<-function (X, I, J, r = NULL, breaks = NULL, ...,
                   k=1, disjoint = NULL,
                   correction = c("rs", "km", "han")){
  verifyclass(X, "ppp")
  W <- X$window
  npts <- npoints(X)
  areaW <- area(W)
  I <- ppsubset(X, I, "I")
  J <- ppsubset(X, J, "J")
  if (is.null(I) || is.null(J))
    stop("I and J must be valid subset indices")
  nI <- sum(I)
  nJ <- sum(J)
  if (nI == 0)
    stop("No points satisfy condition I")
  if (nJ == 0)
    stop("No points satisfy condition J")
  if (is.null(disjoint))
    disjoint <- !any(I & J)
  if (is.null(correction))
    correction <- c("rs", "km", "han")
  correction <- pickoption("correction", correction, c(none = "none",
                                                       border = "rs", rs = "rs", KM = "km", km = "km", Kaplan = "km",
                                                       han = "han", Hanisch = "han", best = "km"), multi = TRUE)
  lamJ <- nJ/areaW
  rmaxdefault <- rmax.rule("G", W, lamJ)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault = rmaxdefault)
  rmax <- breaks$max
  rvals <- breaks$r
  zeroes <- numeric(length(rvals))
  df <- data.frame(r = rvals, theo = 1 - exp(-lamJ * pi * rvals^2))
  fname <- c("ND", "list(I,J)")
  Z <- fv(df, "r", substitute({ND[I, J]^k}(r), list(k=k)), "theo", . ~ r, c(0, rmax),
          c("r", makefvlabel(NULL, NULL, fname, "pois")), c("distance argument r",
                                                            "theoretical Poisson %s"), fname = fname, yexp = substitute({ND[list(I, J)]^{k}}(r),list(k=k)))
  XI <- X[I]
  XJ <- X[J]
  if (disjoint)
    nnd <- nncross(XI, XJ, k=k, what = "dist")
  else {
    seqnp <- seq_len(npts)
    iX <- seqnp[I]
    iY <- seqnp[J]
    nnd <- nncross(XI, XJ, iX, iY, k=k, what = "dist")
  }
  bdry <- bdist.points(XI)
  o <- pmin.int(nnd, bdry)
  d <- (nnd <= bdry)
  if ("none" %in% correction) {
    if (npts == 0)
      edf <- zeroes
    else {
      hh <- hist(nnd[nnd <= rmax], breaks = breaks$val,
                 plot = FALSE)$counts
      edf <- cumsum(hh)/length(nnd)
    }
    Z <- bind.fv(Z, data.frame(raw = edf), makefvlabel(NULL,
                                                       "hat", fname, "raw"), "uncorrected estimate of %s",
                 "raw")
  }
  if ("han" %in% correction) {
    if (npts == 0)
      G <- zeroes
    else {
      x <- nnd[d]
      a <- eroded.areas(W, rvals)
      h <- hist(x[x <= rmax], breaks = breaks$val, plot = FALSE)$counts
      G <- cumsum(h/a)
      G <- G/max(G[is.finite(G)])
    }
    Z <- bind.fv(Z, data.frame(han = G), makefvlabel(NULL,
                                                     "hat", fname, "han"), "Hanisch estimate of %s", "han")
    attr(Z, "alim") <- range(rvals[G <= 0.9])
  }
  if (any(correction %in% c("rs", "km"))) {
    if (npts == 0)
      result <- data.frame(rs = zeroes, km = zeroes, hazard = zeroes)
    else {
      result <- km.rs(o, bdry, d, breaks)
      result <- as.data.frame(result[c("rs", "km", "hazard")])
    }
    Z <- bind.fv(Z, result, c(makefvlabel(NULL, "hat", fname,
                                          "bord"), makefvlabel(NULL, "hat", fname, "km"), "hazard(r)"),
                 c("border corrected estimate of %s", "Kaplan-Meier estimate of %s",
                   "Kaplan-Meier estimate of hazard function lambda(r)"),
                 "km")
    attr(Z, "alim") <- range(rvals[result$km <= 0.9])
  }
  nama <- names(Z)
  fvnames(Z, ".") <- rev(nama[!(nama %in% c("r", "hazard"))])
  unitname(Z) <- unitname(X)
  return(Z)
}
