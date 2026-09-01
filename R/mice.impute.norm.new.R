#' Imputation by Bayesian linear regression
#'
#' Calculates imputations for univariate missing data by Bayesian linear
#' regression, also known as the normal model.
#'
#' @aliases mice.impute.norm norm
#' @inheritParams mice.impute.pmm
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @author Stef van Buuren, Karin Groothuis-Oudshoorn
#' @details
#' Imputation of \code{y} by the normal model by the method defined by
#' Rubin (1987, p. 167). The procedure is as follows:
#'
#' \enumerate{
#' \item{Calculate the cross-product matrix \eqn{S=X_{obs}'X_{obs}}.}
#' \item{Calculate \eqn{V = (S+{diag}(S)\kappa)^{-1}}, with some small ridge
#' parameter \eqn{\kappa}.}
#' \item{Calculate regression weights \eqn{\hat\beta = VX_{obs}'y_{obs}.}}
#' \item{Draw a random variable \eqn{\dot g \sim \chi^2_\nu} with \eqn{\nu=n_1 - q}.}
#' \item{Calculate \eqn{\dot\sigma^2 = (y_{obs} - X_{obs}\hat\beta)'(y_{obs} - X_{obs}\hat\beta)/\dot g.}}
#' \item{Draw \eqn{q} independent \eqn{N(0,1)} variates in vector \eqn{\dot z_1}.}
#' \item{Calculate \eqn{V^{1/2}} by Cholesky decomposition.}
#' \item{Calculate \eqn{\dot\beta = \hat\beta + \dot\sigma\dot z_1 V^{1/2}}.}
#' \item{Draw \eqn{n_0} independent \eqn{N(0,1)} variates in vector \eqn{\dot z_2}.}
#' \item{Calculate the \eqn{n_0} values \eqn{y_{imp} = X_{mis}\dot\beta + \dot z_2\dot\sigma}.}
#' }
#'
#' Using \code{mice.impute.norm} for all columns emulates Schafer's NORM method (Schafer, 1997).
#' @references
#' Rubin, D.B (1987). Multiple Imputation for Nonresponse in Surveys. New York: John Wiley & Sons.
#'
#' Schafer, J.L. (1997). Analysis of incomplete multivariate data. London: Chapman & Hall.
#' @family univariate imputation functions
#' @keywords datagen
#' @export
mice.impute.norm.new <- function(y, ry, x, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .norm.draw.new(y, ry, x, ...)
  x[wy, ] %*% parm$beta + rnorm(sum(wy)) * parm$sigma
}


#' Draws values of beta and sigma by Bayesian linear regression
#'
#' This function draws random values of beta and sigma under the Bayesian
#' linear regression model as described in Rubin (1987, p. 167). This function
#' can be called by user-specified imputation functions.
#'
#' @aliases norm.draw .norm.draw
#' @param y Incomplete data vector of length \code{n}
#' @param ry Vector of missing data pattern (\code{FALSE}=missing,
#' \code{TRUE}=observed)
#' @param x Matrix (\code{n} x \code{p}) of complete covariates.
#' @param rank.adjust Argument that specifies whether \code{NA}'s in the
#' coefficients need to be set to zero. Only relevant when \code{ls.meth = "qr"}
#' AND the predictor matrix is rank-deficient.
#' @param ... Other named arguments.
#' @return A \code{list} containing components \code{coef} (least squares estimate),
#' \code{beta} (drawn regression weights) and \code{sigma} (drawn value of the
#' residual standard deviation).
#' @references
#' Rubin, D.B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
#' @author Gerko Vink, 2018, for this version, based on earlier versions written
#' by Stef van Buuren, Karin Groothuis-Oudshoorn, 2017
#' @export
norm.draw.new <- function(y, ry, x, ...) {
  return(.norm.draw.new(y, ry, x, ...))
}

###' @rdname norm.draw
###' @export
.norm.draw.new <- function(y, ry, x, ...) {
  return(estimice2(x[ry, , drop = FALSE], y[ry], ...))
}


#' Computes least squares parameters
#'
#' This function computes least squares estimates, variance/covariance matrices,
#' residuals and degrees of freedom according to ridge regression, QR decomposition
#' or Singular Value Decomposition. This function is internally called by .norm.draw(),
#' but can be called by any user-specified imputation function.
#'
#' When calculating the inverse of the crossproduct of the predictor matrix,
#' problems may arise. For example, taking the inverse is not possible when the
#' predictor matrix is rank deficient, or when the estimation problem is
#' computationally singular. This function detects such error cases and
#' automatically falls back to adding a ridge penalty to the diagonal of the
#' crossproduct to allow for proper calculation of the inverse.
#'
#' @note
#' This functions adds a star to variable names in the mice iteration
#' history to signal that a ridge penalty was added. In that case, it
#' also adds an entry to \code{loggedEvents}.
#'
#' @aliases estimice
#' @param x Matrix (\code{n} x \code{p}) of complete covariates.
#' @param y Incomplete data vector of length \code{n}
#' @param ls.meth the method to use for obtaining the least squares estimates. By
#' default parameters are drawn by means of QR decomposition.
#' @param ridge A small numerical value specifying the size of the ridge used.
#' The default value \code{ridge = 1e-05} represents a compromise between stability
#' and unbiasedness. Decrease \code{ridge} if the data contain many junk variables.
#' Increase \code{ridge} for highly collinear data.
#' @param ... Other named arguments.
#' @return A \code{list} containing components \code{c} (least squares estimate),
#' \code{r} (residuals), \code{v} (variance/covariance matrix) and \code{df}
#' (degrees of freedom).
#' @author Thom Volker, 2025
#' @export
estimice2 <- function(x, y, ls.meth = "qr", ridge = 1e-05, ...) {
  dim <- dim(x)

  if (ls.meth == "qr") {
    fit <- .lm.fit(x = x, y = y)
    df <- dim[1] - fit$rank
    rankdef <- is_rank_deficient(fit$rank, dim[1], dim[2])
    resids <- fit$residuals

    if (rankdef) {
      if (fit$rank < dim[1]) {
        # collinear system but p < n; kick out collinear predictors
        piv <- match(seq_len(dim[2]), fit$pivot)
        coefs <- fit$coefficients[piv] # re-order coefficients such that correct ones are set to zero
        S <- matrix(0, dim[2], dim[2])
        S[seq_len(fit$rank), seq_len(fit$rank)] <- chol(QR_to_inv(
          fit$qr,
          fit$rank
        ))
        Sigma <- S[piv, piv, drop = FALSE]
      } else {
        ls.meth <- "ridge"
      }
      message <- paste0(
        "Your data are super-collinear. See `?is_rank_deficient` for more information."
      )
      #updateLog(out = mess, frame = 6)
      #if (get.printFlag()) {
      #  cat("*")
      #} # indicator of added ridge penalty in the printed iteration history
    } else {
      coefs <- fit$coefficients
      Sigma <- chol(QR_to_inv(fit$qr, fit$rank))
    }
  }
  if (ls.meth == "ridge") {
    pen <- ridge * apply(x, 2, var) * sqrt(dim[1] * (dim[1]-1) / var(y)) # same scaling as glmnet
    fit <- .lm.fit(
      rbind(x, diag(sqrt(pen))),
      c(y, rep(0, dim[2]))
    )
    coefs <- fit$coefficients
    inv <- QR_to_inv(fit$qr, fit$rank)
    Sigma <- chol(inv)
    df <- pmax(1, dim[1] - sum(diag(x %*% inv %*% t(x)))) # penalty matrix is diagonal, so t(pen) == pen
    resids <- fit$residuals
  }
  if (ls.meth == "svd") {
    fit <- NULL
    eps <- sqrt(.Machine$double.eps)
    s <- svd(x)
    rank <- sum(s$d > eps)
    if (rank == dim[1]) {
      dinv <- ifelse(s$d > eps, s$d / (s$d^2 + ridge), 0)
    } else {
      dinv <- ifelse(s$d > eps, 1/s$d, 0)
    }
    coefs <- s$v %*% (crossprod(s$u, y) * dinv)
    resids <- y - x %*% coefs
    Sigma <- t(s$v[,seq_len(rank),drop = FALSE] %*% diag(dinv[seq_len(rank)], nrow = rank))
    df <- pmax(dim[1] - rank, 1)
    dim[2] <- nrow(Sigma) # draw coefficients in subspace
  }
  sigma.star <- sqrt(c(crossprod(resids)) / rchisq(1, df))
  Sigma <- sigma.star * Sigma
  out <- list(
    coef = as.matrix(coefs),
    beta = as.matrix(coefs + c(crossprod(Sigma, rnorm(dim[2])))),
    sigma = sigma.star,
    estimation = ls.meth
  )
  out
}


get.printFlag <- function(start = 4) {
  while (inherits(
    try(get("printFlag", parent.frame(start)), silent = TRUE),
    "try-error"
  )) {
    start <- start + 1
  }
  get("printFlag", parent.frame(start))
}

is_rank_deficient <- function(rank, n, p) {
  rank < p | rank == n
}


QR_to_inv <- function(qr, rank) {
  R <- qr[1:rank, 1:rank, drop = FALSE]
  R[lower.tri(R)] <- 0
  chol2inv(R)
}


