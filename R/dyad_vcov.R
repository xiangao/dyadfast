#' Fast dyadic-robust variance estimation
#'
#' Computes the Aronow, Samii & Assenova (2015) cluster-robust variance
#' estimator for dyadic data via a single O(nK) scatter-add pass, replacing
#' the O(m * n * K) agent loop in dyadRobust.
#'
#' The efficient formula is:
#'   Vhat = brd \%*\% (H'H - D'D) / n^2 \%*\% brd
#' where S is the n x K score matrix (estfun), D is the d x K matrix of
#' dyad-level score sums (one row per unique dyad), H is the m x K matrix of
#' agent-level score sums (H[i,] = sum of scores for all dyads involving i),
#' and brd = n * (X'X)^(-1) is the bread.
#'
#' For cross-sectional data where each dyad appears once, D'D = S'S and the
#' \code{dyad} argument can be omitted. For panel data where the same dyad
#' appears across multiple time periods, \code{dyad} must be provided so that
#' scores are correctly summed within each dyad before subtraction.
#'
#' @param fit Model object. Supports any class compatible with
#'   \code{sandwich::estfun()} and \code{sandwich::bread()} (e.g., \code{lm}).
#' @param ego Length-n vector of "i" agent identifiers for each observation.
#' @param alter Length-n vector of "j" agent identifiers for each observation.
#' @param dyad Optional length-n vector of dyad (pair) identifiers. Required
#'   for panel data where the same pair appears in multiple observations.
#'   When NULL (default), each observation is treated as its own dyad
#'   (equivalent to cross-sectional data).
#'
#' @importFrom sandwich estfun bread sandwich
#' @importFrom stats coef
#'
#' @return A list with components:
#'   \item{bhat}{Named vector of coefficient estimates.}
#'   \item{sehat}{Named vector of dyadic-robust standard errors.}
#'   \item{Vhat}{K x K dyadic-robust covariance matrix.}
#'   \item{N}{Number of agents m.}
#'   \item{n}{Number of dyadic observations n.}
#'
#' @examples
#' # Simulate a small undirected dyadic dataset (N=10 agents, 45 dyads)
#' set.seed(1)
#' N <- 10
#' pairs <- subset(expand.grid(ego = 1:N, alter = 1:N), ego < alter)
#' pairs$x <- rnorm(nrow(pairs))
#' pairs$y <- 0.5 * pairs$x + rnorm(nrow(pairs))
#' fit <- lm(y ~ x, data = pairs)
#' res <- dyad_vcov(fit, ego = pairs$ego, alter = pairs$alter)
#' res$sehat
#'
#' @export
dyad_vcov <- function(fit, ego, alter, dyad = NULL) {
  S   <- sandwich::estfun(fit)   # n x K score matrix (X_i * e_i)
  brd <- sandwich::bread(fit)    # K x K: n * (X'X)^{-1}
  n   <- nrow(S)
  K   <- ncol(S)

  # Convert identifiers to character for rowsum() grouping
  ego_ch   <- as.character(ego)
  alter_ch <- as.character(alter)
  agents   <- sort(unique(c(ego_ch, alter_ch)))
  m        <- length(agents)

  # Scatter-add: H[i,] = sum of scores for all dyads involving agent i
  # rowsum(M, group) returns one row per unique group value
  H_ego   <- rowsum(S, ego_ch,   reorder = TRUE)   # rows labelled by ego id
  H_alter <- rowsum(S, alter_ch, reorder = TRUE)   # rows labelled by alter id

  H <- matrix(0.0, nrow = m, ncol = K, dimnames = list(agents, colnames(S)))
  H[rownames(H_ego),   ] <- H[rownames(H_ego),   ] + H_ego
  H[rownames(H_alter), ] <- H[rownames(H_alter), ] + H_alter

  # Core K x K products
  HtH <- crossprod(H)   # H'H

  # Dyad-level subtraction: D'D where D = scores summed within each dyad
  # For cross-sectional data (no dyad arg), D = S so D'D = S'S
  if (is.null(dyad)) {
    DtD <- crossprod(S)
  } else {
    D <- rowsum(S, as.character(dyad), reorder = FALSE)
    DtD <- crossprod(D)
  }

  # meat = (H'H - D'D) / n  so that sandwich() gives brd %*% meat %*% brd / n
  # = brd %*% (H'H - D'D) / n^2 %*% brd  (the desired formula)
  meat <- (HtH - DtD) / n
  Vhat <- sandwich::sandwich(fit, bread. = brd, meat. = meat)

  # Enforce positive semi-definiteness
  ev <- eigen(Vhat, symmetric = TRUE)
  n_neg <- sum(ev$values < 0)
  if (n_neg > 0) {
    warning(sprintf("%d negative eigenvalue(s) set to zero.", n_neg))
    ev$values[ev$values < 0] <- 0
    Vhat <- ev$vectors %*% diag(ev$values, K) %*% t(ev$vectors)
    colnames(Vhat) <- rownames(Vhat) <- colnames(brd)
  }

  sehat <- sqrt(diag(Vhat))
  bhat  <- coef(fit)

  list(bhat = bhat, sehat = sehat, Vhat = Vhat, N = m, n = n)
}
