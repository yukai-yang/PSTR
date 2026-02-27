#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
#################################################################################


# compute the first order derivative of dg/dgamma, dg/dc
# input:
#	vg, vector of the transition functions
#	vs, vector of the transition variables
#	vp, est. of nonlinear parameters, vp[1] gamma, otherwise c
# output: matrix of derivatives, row=length
DerGFunc <- function(vg,vs,vp)
{
  gamma = vp[1]; cc = vp[2:length(vp)]
  tmp1 = vg * (1-vg)
  tmp2 = matrix(vs, length(vs), length(cc))
  tmp2 = t(tmp2) - cc
  
  ret = tmp1 * apply(tmp2,2,prod)
  
  ftmp <- function(iter){
    tmp3 = tmp2; tmp3[iter,] = 1
    return(- tmp1 * apply(tmp3,2,prod) * gamma)
  }
  
  ret = cbind(ret, sapply(1:length(cc),ftmp))
  return(ret)
}



# compute the first and second order derivative of dg/dgamma, dg/dc
# input:
#	vg, vector of the transition functions
#	vs, vector of the transition variables
#	vp, est. of nonlinear parameters, vp[1] gamma, otherwise c
# output: matrix of derivatives, row=length
Der2GFunc <- function(vg,vs,vp)
{
  gamma = vp[1]; cc = vp[2:length(vp)]
  tmp1 = vg * (1-vg) # g^2 * zeta
  tmp2 = t(matrix(vs, length(vs), length(cc))) - cc # s - c
  tmp3 = apply(tmp2,2,prod) # prod all
  
  de1 = tmp1 * tmp3
  
  ftmp <- function(iter){
    tmp = tmp2; tmp[iter,] = 1
    return(apply(tmp,2,prod))
  }
  tmp4 = sapply(1:length(cc),ftmp) # prod without k
  tmp4c = c(tmp4) # vector version of tmp4
  
  de1 = cbind(de1, - tmp1 * tmp4c * gamma) # columns are the parameters
  
  de2 = de1[,1] * (1-2*vg) * tmp3 # d^2 g / d gamma^2
  
  # d^2 g / d gamma d c
  de2 = cbind(de2, 2*(1-vg) * de1[,2:ncol(de1)] * tmp3 + tmp1 * tmp3 * gamma * tmp4c - tmp1 * tmp4c)
  
  # d^2 g / dc dc' vec half
  for(iter in 2:ncol(de1)){
    de2 = cbind(de2, (2*vg-1) * de1[,iter] * gamma * tmp4[,iter-1]) # d^2 g / d c^2
    if(iter<ncol(de1)) for(jter in (iter+1):ncol(de1)){
      de2 = cbind(de2,-2*(1-vg)*de1[,jter]*gamma*tmp4[,iter-1]+de1[,iter]*gamma*tmp4[,jter-1]+(1-vg)*gamma*tmp4[,iter-1]/tmp2[jter-1,])
    }
  }
  
  ret = list(de1=de1,de2=de2)
  return(ret)
}


#' Estimate a PSTR model by nonlinear least squares
#'
#' \code{EstPSTR} estimates either a nonlinear PSTR model (when \code{iq} is provided) or a
#' linear fixed-effects panel regression (when \code{iq = NULL}).
#'
#' Two equivalent interfaces are available:
#' \enumerate{
#'   \item Wrapper function: \code{EstPSTR(use = obj, ...)}.
#'   \item R6 method: \code{obj$EstPSTR(...)}.
#' }
#' The wrapper calls the corresponding R6 method and returns \code{use} invisibly.
#'
#' The transition function is logistic and depends on a transition variable \eqn{q_{it}} and
#' nonlinear parameters \eqn{\gamma > 0} and switching locations \eqn{c_1 < \cdots < c_m}:
#' \deqn{g(q_{it}; \gamma, c_1,\ldots,c_m) = \left(1 + \exp\left[-\gamma \prod_{j=1}^{m}(q_{it}-c_j)\right]\right)^{-1}.}
#' The smoothness parameter is internally reparametrised as \eqn{\gamma = \exp(\delta)}, where
#' \eqn{\delta \in \mathbb{R}}. The optimisation is always carried out in \eqn{\delta} and \eqn{c}.
#'
#' If \code{par = NULL}, the function constructs default initial values from quantiles of the
#' selected transition variable and treats the first element as \eqn{\delta}.
#'
#' @param use An object of class \code{"PSTR"} created by \code{\link{NewPSTR}}.
#' @param im Integer. Number of switches \eqn{m} in the transition function. Default is \code{1}.
#' @param iq Either an integer index (column number in the transition-variable matrix) or a
#'   character string (transition-variable name) specifying which transition variable to use.
#'   If \code{NULL}, a linear fixed-effects panel regression is estimated.
#' @param par Numeric vector of length \code{im + 1} giving initial values for the nonlinear
#'   parameters. The expected order is \code{c(delta, c_1, ..., c_m)} if \code{useDelta = TRUE},
#'   or \code{c(gamma, c_1, ..., c_m)} if \code{useDelta = FALSE}. If \code{NULL}, defaults are
#'   constructed automatically and \code{useDelta} is ignored.
#' @param useDelta Logical. If \code{TRUE}, the first element of \code{par} is interpreted as
#'   \eqn{\delta}. If \code{FALSE}, it is interpreted as \eqn{\gamma} and internally converted
#'   to \eqn{\delta = \log(\gamma)} before optimisation.
#' @param vLower Numeric scalar or vector. Lower offsets defining the lower bounds in the optimiser.
#'   Bounds are applied to the internal parameter vector used in optimisation (with the first
#'   element being \eqn{\delta}).
#' @param vUpper Numeric scalar or vector. Upper offsets defining the upper bounds in the optimiser.
#'   Bounds are applied to the internal parameter vector used in optimisation (with the first
#'   element being \eqn{\delta}).
#' @param method Character. Optimisation method passed to \code{\link[stats:optim]{stats::optim}}.
#'   Default is \code{"L-BFGS-B"} (bounded optimisation).
#'
#' @return Invisibly returns \code{use} with estimation results added. In particular, for a
#'   nonlinear PSTR model (\code{iq} not \code{NULL}), the object contains (among others):
#' \describe{
#'   \item{\code{delta}}{Estimate of \eqn{\delta}.}
#'   \item{\code{gamma}}{Estimate of \eqn{\gamma = \exp(\delta)}.}
#'   \item{\code{c}}{Estimates of \eqn{c_1,\ldots,c_m}.}
#'   \item{\code{vg}}{Estimated transition-function values \eqn{g_{it}}.}
#'   \item{\code{beta}}{Estimated coefficients (named as \code{var_0} for linear-part coefficients and \code{var_1} for nonlinear-part coefficients).}
#'   \item{\code{vU}}{Residuals.}
#'   \item{\code{vM}}{Estimated individual effects.}
#'   \item{\code{s2}}{Estimated residual variance.}
#'   \item{\code{cov}}{Cluster-robust and heteroskedasticity-consistent covariance matrix of all estimates.}
#'   \item{\code{se}}{Standard errors corresponding to \code{est}.}
#'   \item{\code{est}}{Vector of all estimates (coefficients followed by nonlinear parameters).}
#'   \item{\code{mbeta}}{Estimates of coefficients in the second extreme regime (when available).}
#'   \item{\code{mse}}{Standard errors for \code{mbeta} (when available).}
#' }
#' For a linear fixed-effects model (\code{iq = NULL}), the object contains \code{beta}, \code{vU},
#' \code{vM}, \code{s2}, \code{cov}, \code{se}, and \code{est}.
#'
#' @seealso \code{\link{NewPSTR}}, \code{\link{LinTest}}, \code{\link{WCB_LinTest}},
#'   \code{\link{EvalTest}}, \code{\link[stats:optim]{stats::optim}}.
#'
#' @examples
#' \donttest{
#' pstr <- NewPSTR(Hansen99, dep = "inva", indep = 4:20,
#'                indep_k = c("vala","debta","cfa","sales"),
#'                tvars = c("vala"), iT = 14)
#'
#' # 1) Linear fixed-effects model
#' pstr <- EstPSTR(use = pstr)
#' print(pstr, mode = "estimates", digits = 6)
#'
#' # 2) Nonlinear PSTR model
#' pstr <- EstPSTR(use = pstr, im = 1, iq = 1, useDelta = TRUE,
#'                par = c(.63, 0), vLower = 4, vUpper = 4)
#' print(pstr, mode = "estimates", digits = 6)
#'
#' # R6 method interface (equivalent)
#' pstr$EstPSTR(im = 1, iq = 1, useDelta = TRUE, par = c(.63, 0), method = "CG")
#' }
#'
#' @name EstPSTR
NULL


PSTR$set("public", "EstPSTR", function(im=1, iq=NULL, par=NULL, useDelta=FALSE, vLower=2, vUpper=2, method='L-BFGS-B'){
  iT <- private$iT
  iN <- private$iN
  
  # data
  vY  <- private$vY
  vYb <- private$vYb
  mX  <- private$mX
  mXb <- private$mXb
  mK  <- private$mK
  ik  <- ncol(mK)
  
  ftmp <- function(vx) vx - mean(vx)
  
  # store estimation settings/results in private
  private$imm <- im
  private$iq  <- iq
  
  if(!is.null(iq)){
    
    if(im < 1) stop(simpleError("The number of switches is invalid."))
    
    # resolve iq if it is a name
    if(!is.numeric(iq)) private$iq <- which(private$mQ_name == iq)
    if(length(private$iq) > 1) stop(simpleError("Sorry! We only support the one transition variable case."))
    
    # IMPORTANT: use resolved iq
    vQ <- private$mQ[, private$iq]
    mQ <- t(matrix(vQ, iT*iN, im))
    
    ResiduleSumSquare <- function(vp){
      # vp[1] = log(gamma) or delta
      vg <- fTF(vx=mQ, gamma=exp(vp[1]), vc=vp[2:length(vp)])
      mXX <- mK * vg
      aXX <- array(c(mXX), dim=c(iT,iN,ik))
      mXXb <- cbind(mXb, matrix(c(apply(aXX, c(2,3), ftmp)), iT*iN, ik))
      tmp <- qr.solve(mXXb, vYb)
      vE <- c(vYb - mXXb %*% tmp)
      sum(vE*vE) / iT / iN
    }
    
    if(is.null(par)){
      useDelta <- TRUE
      tmp <- unname(quantile(vQ, (1:im) / (im+1)))
      par <- c(log(8/min(diff(c(0,tmp)))), tmp)
    }
    
    if(!useDelta) par[1] <- log(par[1])
    private$par <- par
    
    if(method == 'L-BFGS-B'){
      opt <- optim(par=par, fn=ResiduleSumSquare, method="L-BFGS-B",
                   lower=par-vLower, upper=par+vUpper)
    } else {
      opt <- optim(par=par, fn=ResiduleSumSquare, method=method)
    }
    
    # return value (store to private)
    self$delta <- opt$par[1]
    self$gamma <- exp(self$delta)
    self$c <- opt$par[2:length(opt$par)]
    private$convergence <- opt$convergence
    
    vg <- fTF(vx=mQ, gamma=self$gamma, vc=self$c) # g_it
    self$vg <- vg
    
    mXX <- mK * vg
    aXX <- array(c(mXX), dim=c(iT,iN,ik))
    mXXb <- cbind(mXb, matrix(c(apply(aXX, c(2,3), ftmp)), iT*iN, ik))
    
    tmp <- qr.solve(mXXb, vYb)
    self$beta <- c(tmp)
    names(self$beta) <- c(paste0(private$mX_name,'_0'),
                             paste0(private$mK_name,'_1'))
    
    mXX <- cbind(mX, mXX)
    private$mXX <- mXX
    
    mtmp <- matrix(c(vY - mXX %*% tmp), iT, iN)
    self$vM <- c(apply(mtmp, 2, mean))
    self$vU <- c(apply(mtmp, 2, ftmp))
    self$s2 <- c(self$vU %*% self$vU) / (iT*iN)
    
    # computing standard errors
    tg <- Der2GFunc(vg=vg, vs=vQ, vp=c(self$gamma, self$c))
    de1 <- tg$de1
    de2 <- tg$de2
    
    beta1 <- self$beta[(ncol(mX)+1):length(self$beta)]
    
    dedp <- -mXXb
    d2edp2 <- array(0, dim=c(iT*iN, length(self$beta)+1+im, length(self$beta)+1+im))
    
    tcnt <- 1
    for(iter in 1:ncol(de1)){
      mKK <- mK * de1[,iter]
      aKK <- array(c(mKK), dim=c(iT,iN,ik))
      mKK <- matrix(c(apply(aKK, c(2,3), ftmp)), iT*iN, ik)
      
      dedp <- cbind(dedp, -mKK %*% beta1)
      
      d2edp2[,(ncol(mX)+1):length(self$beta), length(self$beta)+iter] <- -mKK
      d2edp2[, length(self$beta)+iter, (ncol(mX)+1):length(self$beta)] <- -mKK
      
      for(jter in iter:ncol(de1)){
        mKK <- mK * de2[,tcnt]
        aKK <- array(c(mKK), dim=c(iT,iN,ik))
        mKK <- matrix(c(apply(aKK, c(2,3), ftmp)), iT*iN, ik)
        
        d2edp2[, length(self$beta)+iter, length(self$beta)+jter] <- -mKK %*% beta1
        d2edp2[, length(self$beta)+jter, length(self$beta)+iter] <- -mKK %*% beta1
        tcnt <- tcnt + 1
      }
    }
    
    mh <- 2 * self$vU * dedp
    ah <- array(c(mh), dim=c(iT,iN,ncol(dedp)))
    hi <- matrix(c(apply(ah, c(2,3), sum)), iN, ncol(dedp))
    
    mB <- 0
    for(iter in 1:iN) mB <- mB + hi[iter,] %*% t(hi[iter,])
    
    invA <- 0
    for(iter in 1:(iT*iN))
      invA <- invA + (dedp[iter,] %*% t(dedp[iter,]) + d2edp2[iter,,] * self$vU[iter]) * 2
    
    ttmp <- try(solve(invA), silent=TRUE)
    if(inherits(ttmp, 'try-error')){
      s <- svd(invA)
      invA <- s$u %*% diag(1/s$d) %*% t(s$u)
    } else {
      invA <- ttmp
    }
    
    self$cov <- invA %*% mB %*% t(invA)
    self$se <- sqrt(diag(self$cov))
    names(self$se) <- c(names(self$beta), 'gamma', paste0('c_', 1:im))
    
    self$est <- c(self$beta, self$gamma, self$c)
    names(self$est) <- names(self$se)
    
    mM <- NULL
    mname <- NULL
    mTmp <- diag(length(private$mX_name))
    
    for(iter in 1:length(private$mX_name)){
      idx <- private$mX_name[iter] == private$mK_name
      if(any(idx)){
        mM <- rbind(mM, c(mTmp[iter,], idx))
        mname <- c(mname, private$mX_name[iter])
      }
    }
    
    if(!is.null(mM)){
      mM <- cbind(mM, matrix(0, nrow(mM), 1+im))
      self$mbeta <- c(mM %*% self$est)
      names(self$mbeta) <- mname
      self$mse <- sqrt(diag(mM %*% self$cov %*% t(mM)))
      names(self$mse) <- mname
    }
    
  } else {
    
    tmp <- qr.solve(mXb, vYb)
    
    self$beta <- c(tmp)
    names(self$beta) <- private$mX_name
    
    mtmp <- matrix(c(vY - mX %*% tmp), iT, iN)
    self$vM <- c(apply(mtmp, 2, mean))
    self$vU <- c(apply(mtmp, 2, ftmp))
    self$s2 <- c(self$vU %*% self$vU) / (iT*iN)
    
    # computing standard errors
    dedp <- -mXb
    mh <- 2 * self$vU * dedp
    ah <- array(c(mh), dim=c(iT,iN,ncol(dedp)))
    hi <- matrix(c(apply(ah, c(2,3), sum)), iN, ncol(dedp))
    
    mB <- 0
    for(iter in 1:iN) mB <- mB + hi[iter,] %*% t(hi[iter,])
    
    invA <- 0
    for(iter in 1:(iT*iN))
      invA <- invA + dedp[iter,] %*% t(dedp[iter,]) * 2
    
    ttmp <- try(solve(invA), silent=TRUE)
    if(inherits(ttmp, 'try-error')){
      s <- svd(invA)
      invA <- s$u %*% diag(1/s$d) %*% t(s$u)
    } else {
      invA <- ttmp
    }
    
    self$cov <- invA %*% mB %*% t(invA)
    self$se <- sqrt(diag(self$cov))
    names(self$se) <- names(self$beta)
    
    self$est <- self$beta
    names(self$est) <- names(self$se)
  }
  
  invisible(self)
})

#' @rdname EstPSTR
#' @export
EstPSTR <- function(use, im=1, iq=NULL, par=NULL, useDelta=FALSE, vLower=2, vUpper=2, method='L-BFGS-B')
{
  if(!inherits(use, 'PSTR'))
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  
  use$EstPSTR(im=im, iq=iq, par=par, useDelta=useDelta, vLower=vLower, vUpper=vUpper, method=method)
  invisible(use)
}
           