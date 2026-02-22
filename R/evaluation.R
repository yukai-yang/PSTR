#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
#################################################################################


#' Evaluate an estimated PSTR model
#'
#' \code{EvalTest} provides post-estimation evaluation tests for an estimated PSTR model.
#' It supports two null hypotheses:
#' \describe{
#'   \item{Parameter constancy}{No time variation in parameters (labelled \code{"time-varying"}).}
#'   \item{No remaining nonlinearity}{No remaining nonlinearity/heterogeneity given a candidate transition variable (labelled \code{"heterogeneity"}).}
#' }
#'
#' Wild bootstrap (WB) and wild cluster bootstrap (WCB) versions are available via
#' \code{WCB_TVTest} (parameter constancy) and \code{WCB_HETest} (no remaining nonlinearity).
#'
#' Two equivalent interfaces are available for each test:
#' \enumerate{
#'   \item Wrapper function, for example \code{EvalTest(use = obj, ...)}.
#'   \item R6 method, for example \code{obj$EvalTest(...)}.
#' }
#' Each wrapper calls the corresponding R6 method and returns \code{use} invisibly.
#'
#' The bootstrap variants are computationally intensive. WB is robust to heteroskedasticity,
#' while WCB is additionally robust to within-individual dependence (cluster dependence).
#' Parallel execution can be enabled via \code{parallel} and \code{cpus}.
#'
#' @param use An object of class \code{"PSTR"} returned by \code{\link{EstPSTR}}.
#'   The model must be estimated (nonlinear PSTR) before evaluation tests can be run.
#' @param type Character vector. Which evaluation tests to run in \code{EvalTest}.
#'   Must be a subset of \code{c("time-varying","heterogeneity")}. Default is both.
#' @param vq Numeric vector. Candidate transition variable used by the no remaining nonlinearity
#'   test. Required if \code{"heterogeneity"} is included in \code{type}, and required for
#'   \code{WCB_HETest}.
#' @param iB Integer. Number of bootstrap replications. Default is \code{100}.
#' @param parallel Logical. Whether to use parallel computation (via the \pkg{snowfall} backend).
#' @param cpus Integer. Number of CPU cores used if \code{parallel = TRUE}.
#'
#' @return Invisibly returns \code{use} with evaluation results added.
#' \describe{
#'   \item{\code{tv}}{A list of parameter-constancy (time-varying) test results, one element per \eqn{m}.}
#'   \item{\code{ht}}{A list of no remaining nonlinearity (heterogeneity) test results, one element per \eqn{m}.}
#'   \item{\code{wcb_tv}}{A numeric matrix of WB/WCB p-values for parameter-constancy tests (one row per \eqn{m}).}
#'   \item{\code{wcb_ht}}{A numeric matrix of WB/WCB p-values for no remaining nonlinearity tests (one row per \eqn{m}).}
#' }
#' The individual list elements in \code{tv} and \code{ht} contain LM-type test statistics and
#' p-values (including HAC variants), consistent with the output from \code{\link{LinTest}}.
#'
#' @seealso \code{\link{NewPSTR}}, \code{\link{LinTest}}, \code{\link{WCB_LinTest}},
#'   \code{\link{EstPSTR}}.
#'
#' @examples
#' \donttest{
#' pstr <- NewPSTR(Hansen99, dep = "inva", indep = 4:20,
#'                indep_k = c("vala","debta","cfa","sales"),
#'                tvars = c("vala"), iT = 14)
#'
#' # estimate first
#' pstr <- EstPSTR(use = pstr, im = 1, iq = 1, useDelta = TRUE, par = c(.63, 0), method = "CG")
#'
#' # evaluation tests
#' pstr <- EvalTest(
#'   use = pstr,
#'   type = c("time-varying","heterogeneity"),
#'   vq = as.matrix(Hansen99[,'vala'])[,1]
#' )
#' print(pstr, mode = "evaluation")
#'
#' # bootstrap variants (requires snowfall)
#' library(snowfall)
#' pstr <- WCB_TVTest(
#'     use = pstr, iB = 4,
#'     parallel = TRUE, cpus = 2)
#' pstr <- WCB_HETest(
#'     use = pstr,
#'     vq = as.matrix(Hansen99[,'vala'])[,1],
#'     iB = 4, parallel = TRUE, cpus = 2)
#' print(pstr, mode = "evaluation")
#' }
#'
#' @name EvalTest
NULL


PSTR$set("public", "EvalTest", function(type = c("time-varying", "heterogeneity"), vq = NULL) {
  
  type <- match.arg(type, several.ok = TRUE)
  
  if (is.null(private$iq)) {
    stop(simpleError("Estimate the PSTR model first!"))
  }
  
  im = private$im
  
  mD = diag(1,private$iN) %x% rep(1,private$iT)
  mM = diag(1, private$iN*private$iT) - tcrossprod(mD)/private$iT
  
  tmp = c(private$mK %*% self$beta[(ncol(private$mX)+1):length(self$beta)])
  tmp = mD * tmp ## pp.14
  mV = cbind(private$mXX, tmp)
  mV2 = mM %*% mV
  invVV = svd_pinv(crossprod(mV2))
  
  if(length(grep("time-varying",type))>0){
    self$tv = list(); length(self$tv) = im
    
    vt = 1:private$iT/private$iT
    
    mW = NULL
    for(mter in 1:im){
      mW = cbind(mW, private$mXX*(vt**mter))
      self$tv[[mter]] = LMTEST(iT=private$iT,iN=private$iN,vU=self$vU,mX=mV,
                              mW=mW,mM=mM,s2=self$s2,mX2=mV2,invXX=invVV)
    }
  }
  
  if(length(grep("heterogeneity",type))>0){
    self$ht = list(); length(self$ht) = im
    
    mW = NULL
    for(mter in 1:im){
      mW = cbind(mW, private$mXX*(vq**mter))
      self$ht[[mter]] = LMTEST(iT=private$iT,iN=private$iN,vU=self$vU,mX=mV,
                              mW=mW,mM=mM,s2=self$s2,mX2=mV2,invXX=invVV)
    }
  }
  
  cli::cli_alert_success("Done!")
  invisible(self)
})


#' @rdname EvalTest
#' @export
EvalTest <- function(use, type = c("time-varying", "heterogeneity"), vq = NULL) {
  if (!inherits(use, "PSTR")) {
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  }
  use$EvalTest(type = type, vq = vq)
  invisible(use)
}

PSTR$set("public", ".set_vY", function(vY_new) {
  private$vY <- as.numeric(vY_new)
  invisible(self)
})


PSTR$set("public", ".get_mK", function() { private$mK })
PSTR$set("public", ".get_mX", function() { private$mX })
PSTR$set("public", ".get_mXX", function() { private$mXX })


PSTR$set("public", "WCB_TVTest", function(iB = 100, parallel = FALSE, cpus = 4) {
  
  if (is.null(private$iq)) {
    stop(simpleError("Estimate the PSTR model first!"))
  }
  
  # snowfall is only needed when parallel = TRUE
  if (parallel) {
    if (!requireNamespace("snowfall", quietly = TRUE)) {
      stop(
        "Parallel bootstrap requires the 'snowfall' package. ",
        "Please install it via install.packages('snowfall').",
        call. = FALSE
      )
    }
  }
  
  # helper: run bootstrap either in parallel (snowfall) or serial (base sapply)
  .run_boot <- function(FUN) {
    if (parallel) {
      snowfall::sfInit(parallel = TRUE, cpus = cpus)
      on.exit(snowfall::sfStop(), add = TRUE)
      snowfall::sfSapply(1:iB, FUN)
    } else {
      sapply(1:iB, FUN)
    }
  }
  
  # IMPORTANT: use a deep clone so we don't overwrite the original object
  ruse <- self$clone(deep = TRUE)
  
  im <- private$im
  iT <- private$iT
  iN <- private$iN
  
  vU <- self$vU
  eY <- private$vY - vU
  
  mD <- diag(1, iN) %x% rep(1, iT)
  mM <- diag(1, iN * iT) - tcrossprod(mD) / iT
  
  tmp <- c(private$mK %*% self$beta[(ncol(private$mX) + 1):length(self$beta)])
  tmp <- mD * tmp
  mV <- cbind(private$mXX, tmp)
  mV2 <- mM %*% mV
  invVV <- svd_pinv(crossprod(mV2))
  
  ftmp_wb <- function(bter) { # WB
    ve1 <- sample(c(1, -1), iT * iN, replace = TRUE) * vU
    ruse$.set_vY(eY + ve1)
    
    EST <- EstPSTR(
      use = ruse, im = 1, iq = ruse$iq,
      par = c(self$delta, self$c),
      useDelta = TRUE, vLower = 1, vUpper = 1
    )
    vu1 <- EST$vU
    ss1 <- EST$s2
    
    mK <- EST$.get_mK()
    beta <- EST$beta
    beta_k <- tail(beta, ncol(mK))
    tmp1 <- c(mK %*% beta_k)
    
    tmp1 <- mD * tmp1
    mV11 <- cbind(EST$.get_mXX(), tmp1)
    mV12 <- mM %*% mV11
    invVV1 <- svd_pinv(crossprod(mV12))
    
    sLMTEST(iT = iT, iN = iN, vU = vu1, mX = mV11, mW = mW, mM = mM,
            s2 = ss1, mX2 = mV12, invXX = invVV1)
  }
  
  ftmp_wcb <- function(bter) { # WCB
    ve2 <- c(t(matrix(sample(c(1, -1), iN, replace = TRUE), iN, iT))) * vU
    ruse$.set_vY(eY + ve2)
    
    EST <- EstPSTR(
      use = ruse, im = 1, iq = ruse$iq,
      par = c(self$delta, self$c),
      useDelta = TRUE, vLower = 1, vUpper = 1
    )
    vu2 <- EST$vU
    ss2 <- EST$s2
    
    mK <- EST$.get_mK()
    beta <- EST$beta
    beta_k <- tail(beta, ncol(mK))
    tmp2 <- c(mK %*% beta_k)
    
    tmp2 <- mD * tmp2
    mV21 <- cbind(EST$.get_mXX(), tmp2)
    mV22 <- mM %*% mV21
    invVV2 <- svd_pinv(crossprod(mV22))
    
    sLMTEST(iT = iT, iN = iN, vU = vu2, mX = mV21, mW = mW, mM = mM,
            s2 = ss2, mX2 = mV22, invXX = invVV2)
  }
  
  self$wcb_tv <- NULL
  vt <- 1:iT / iT
  mW <- NULL
  
  for (mter in 1:im) {
    
    mW <- cbind(mW, private$mXX * (vt^mter))
    LM <- sLMTEST(iT = iT, iN = iN, vU = vU, mX = mV, mW = mW, mM = mM,
                  s2 = self$s2, mX2 = mV2, invXX = invVV)
    
    qLM1 <- .run_boot(ftmp_wb)
    qLM2 <- .run_boot(ftmp_wcb)
    
    self$wcb_tv <- rbind(self$wcb_tv, c(LM, mean(LM <= qLM1), mean(LM <= qLM2)))
  }
  
  self$wcb_tv <- matrix(self$wcb_tv, nrow = im)
  
  cli::cli_alert_success("Done!")
  invisible(self)
})



#' @rdname EvalTest
#' @export
WCB_TVTest <- function(use, iB = 100, parallel = FALSE, cpus = 4) {
  if (!inherits(use, "PSTR")) {
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  }
  use$WCB_TVTest(iB = iB, parallel = parallel, cpus = cpus)
  invisible(use)
}


PSTR$set("public", "WCB_HETest", function(vq, iB = 100, parallel = FALSE, cpus = 4) {
  
  if (is.null(private$iq)) {
    stop(simpleError("Estimate the PSTR model first!"))
  }
  
  if (is.null(vq)) {
    stop(simpleError("Please provide 'vq' for the heterogeneity test."))
  }
  
  # snowfall is only needed when parallel = TRUE
  if (parallel) {
    if (!requireNamespace("snowfall", quietly = TRUE)) {
      stop(
        "Parallel bootstrap requires the 'snowfall' package. ",
        "Please install it via install.packages('snowfall').",
        call. = FALSE
      )
    }
  }
  
  # helper: run bootstrap either in parallel (snowfall) or serial (base sapply)
  .run_boot <- function(FUN) {
    if (parallel) {
      snowfall::sfInit(parallel = TRUE, cpus = cpus)
      on.exit(snowfall::sfStop(), add = TRUE)
      snowfall::sfSapply(1:iB, FUN)
    } else {
      sapply(1:iB, FUN)
    }
  }
  
  # IMPORTANT: use a deep clone so we don't overwrite the original object
  ruse <- self$clone(deep = TRUE)
  
  im <- private$im
  iT <- private$iT
  iN <- private$iN
  
  vU <- self$vU
  eY <- private$vY - vU
  
  mD <- diag(1, iN) %x% rep(1, iT)
  mM <- diag(1, iN * iT) - tcrossprod(mD) / iT
  
  # build V for the auxiliary regression (pp.14)
  mK0 <- private$mK
  beta0 <- self$beta
  beta_k0 <- tail(beta0, ncol(mK0))
  tmp <- c(mK0 %*% beta_k0)
  
  tmp <- mD * tmp
  mV <- cbind(private$mXX, tmp)
  mV2 <- mM %*% mV
  invVV <- svd_pinv(crossprod(mV2))
  
  ftmp_wb <- function(bter) { # WB
    ve1 <- sample(c(1, -1), iT * iN, replace = TRUE) * vU
    ruse$.set_vY(eY + ve1)
    
    EST <- EstPSTR(
      use = ruse, im = 1, iq = ruse$iq,
      par = c(self$delta, self$c),
      useDelta = TRUE, vLower = 1, vUpper = 1
    )
    
    vu1 <- EST$vU
    ss1 <- EST$s2
    
    mK <- EST$.get_mK()
    beta <- EST$beta
    beta_k <- tail(beta, ncol(mK))
    tmp1 <- c(mK %*% beta_k)
    
    tmp1 <- mD * tmp1
    mV11 <- cbind(EST$.get_mXX(), tmp1)
    mV12 <- mM %*% mV11
    invVV1 <- svd_pinv(crossprod(mV12))
    
    sLMTEST(
      iT = iT, iN = iN, vU = vu1,
      mX = mV11, mW = mW, mM = mM,
      s2 = ss1, mX2 = mV12, invXX = invVV1
    )
  }
  
  ftmp_wcb <- function(bter) { # WCB
    ve2 <- c(t(matrix(sample(c(1, -1), iN, replace = TRUE), iN, iT))) * vU
    ruse$.set_vY(eY + ve2)
    
    EST <- EstPSTR(
      use = ruse, im = 1, iq = ruse$iq,
      par = c(self$delta, self$c),
      useDelta = TRUE, vLower = 1, vUpper = 1
    )
    
    vu2 <- EST$vU
    ss2 <- EST$s2
    
    mK <- EST$.get_mK()
    beta <- EST$beta
    beta_k <- tail(beta, ncol(mK))
    tmp2 <- c(mK %*% beta_k)
    
    tmp2 <- mD * tmp2
    mV21 <- cbind(EST$.get_mXX(), tmp2)
    mV22 <- mM %*% mV21
    invVV2 <- svd_pinv(crossprod(mV22))
    
    sLMTEST(
      iT = iT, iN = iN, vU = vu2,
      mX = mV21, mW = mW, mM = mM,
      s2 = ss2, mX2 = mV22, invXX = invVV2
    )
  }
  
  self$wcb_ht <- NULL
  mW <- NULL
  
  for (mter in 1:im) {
    
    mW <- cbind(mW, private$mXX * (vq ^ mter))
    
    LM <- sLMTEST(
      iT = iT, iN = iN, vU = vU,
      mX = mV, mW = mW, mM = mM,
      s2 = self$s2, mX2 = mV2, invXX = invVV
    )
    
    qLM1 <- .run_boot(ftmp_wb)
    qLM2 <- .run_boot(ftmp_wcb)
    
    self$wcb_ht <- rbind(self$wcb_ht, c(LM, mean(LM <= qLM1), mean(LM <= qLM2)))
  }
  
  self$wcb_ht <- matrix(self$wcb_ht, nrow = im)
  
  cli::cli_alert_success("Done!")
  invisible(self)
})


#' @rdname EvalTest
#' @export
WCB_HETest <- function(use, vq, iB = 100, parallel = FALSE, cpus = 4) {
  if (!inherits(use, "PSTR")) {
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  }
  use$WCB_HETest(vq = vq, iB = iB, parallel = parallel, cpus = cpus)
  invisible(use)
}