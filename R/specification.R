#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## Aug 2023
#################################################################################


# Evaluate the transition function.
#
# This function evaluate the transition function values. It is used by other functions in the package.
#
# If \code{vx} is a matrix, its row number must be equal to the length of \code{vc}.
#
# vx a vector or matrix of the transition variables.
# gamma the smoothness parameter.
# vc a vector of the location parameters, whose length is the number of switches in the transition function.
# return If vx is a vector, then a scalor retured, otherwise a vector.
fTF <- function(vx, gamma, vc)
{
  tmp = matrix(vx-vc, nrow=length(vc))
  tmp = -apply(tmp,2,prod)*gamma
  return(1 / ( exp(tmp) + 1 ))
}
 

# Compute the LM tests and the p-values.
#
# This function computes the LM tests and the p-values. It is used by other functions in the package.
LMTEST <- function(iT, iN, vU, mX, mW, mM, s2, mX2, invXX)
{
  df1 = ncol(mW)
  df2 = iT*iN - df1 - iN - ncol(mX)
  
  mW2 = mM %*% mW
  mXW2 = crossprod(mX2, mW2)
  S1 = ( crossprod(mW2) - t(mXW2) %*% invXX %*% mXW2 ) * s2
  invS1 = try(chol2inv(chol(S1)),silent=T)
  if(inherits(invS1,'try-error')){
    ttmp = svd(S1);	invS1 = ttmp$u %*% diag(1/ttmp$d, length(ttmp$d)) %*% t(ttmp$u)
  }
  
  vW2U = crossprod(mW2, vU)
  LM1_X = c(t(vW2U) %*% invS1 %*% vW2U)
  PV1_X = 1-pchisq(LM1_X,df=df1)
  LM1_F = LM1_X * df2 / (iT*iN*df1)
  PV1_F = 1-pf(LM1_F,df1=df1,df2=df2)
  
  mZ = cbind(mX2, mW2)
  Delta = 0
  for(nter in  0:(iN-1)*iT){
    itmp = (nter+1):(nter+iT)
    tmp = t(mZ[itmp,]) %*% vU[itmp] %*% t(vU[itmp]) %*% mZ[itmp,]
    Delta = Delta + tmp
  }
  tmp = cbind( -crossprod(mXW2, invXX), diag(1, df1) )
  S2 = tmp %*% Delta %*% t(tmp)
  invS2 = try(chol2inv(chol(S2)),silent=T)
  if(inherits(invS2, 'try-error')){
    ttmp = svd(S2); invS2 = ttmp$u %*% diag(1/ttmp$d, length(ttmp$d)) %*% t(ttmp$u)
  }
  
  LM2_X = c(t(vW2U) %*% invS2 %*% vW2U)
  PV2_X = 1-pchisq(LM2_X,df=df1)
  LM2_F = LM2_X * df2 / (iT*iN*df1)
  PV2_F = 1-pf(LM2_F,df1=df1,df2=df2)
  
  return(list(LM1_X=LM1_X, PV1_X=PV1_X, LM1_F=LM1_F, PV1_F=PV1_F,
              LM2_X=LM2_X, PV2_X=PV2_X, LM2_F=LM2_F, PV2_F=PV2_F))
}  


# Compute the LM tests and the p-values.
#
# This function computes the LM tests and the p-values. It is a simplified version of \code{LMTEST} and is used by other functions in the package.
sLMTEST <- function(iT, iN, vU, mX, mW, mM, s2, mX2, invXX)
{
  df1 = ncol(mW)
  df2 = iT*iN - df1 - iN - ncol(mX)
  
  mW2 = mM %*% mW
  mXW2 = crossprod(mX2, mW2)
  
  S1 = ( crossprod(mW2) - t(mXW2) %*% invXX %*% mXW2 ) * s2
  invS1 = try(chol2inv(chol(S1)),silent=T)
  if(inherits(invS1, 'try-error')){
    ttmp = svd(S1); invS1 = ttmp$u %*% diag(1/ttmp$d, length(ttmp$d)) %*% t(ttmp$u)
  }
  
  vW2U = crossprod(mW2, vU)
  LM1_X = c(t(vW2U) %*% invS1 %*% vW2U)
  
  return(LM1_X)
}


#' Linearity (homogeneity) tests for PSTR models
#'
#' These functions conduct linearity (homogeneity) tests against the alternative of
#' a logistic smooth transition component in a Panel Smooth Transition Regression (PSTR) model.
#'
#' \strong{Two equivalent interfaces are available:}
#' \enumerate{
#'   \item \strong{Wrapper functions:} \code{LinTest(use = obj)} and \code{WCB_LinTest(use = obj, ...)}.
#'   \item \strong{R6 methods:} \code{obj$LinTest()} and \code{obj$WCB_LinTest(...)}.
#' }
#' The wrapper functions call the corresponding R6 methods and return the (mutated) object invisibly.
#'
#' The tests are carried out for each potential transition variable specified in \code{tvars}
#' when creating the model via \code{\link{NewPSTR}}. For each transition variable, tests are computed
#' for the number of switches \eqn{m = 1, \ldots, im}, where \eqn{im} is the maximal number of switches.
#'
#' The procedures produce two families of tests:
#' \describe{
#'   \item{(i) Linearity tests for each \eqn{m}}{
#'     For a fixed \eqn{m}, the null hypothesis is
#'     \deqn{H_0^i: \beta_{i} = \beta_{i-1} = \cdots = \beta_{1} = 0, \qquad i = 1, \ldots, m.}
#'   }
#'   \item{(ii) Sequence tests for selecting \eqn{m}}{
#'     These are conditional tests with null
#'     \deqn{H_0^i: \beta_{i} = 0 \mid \beta_{i+1} = \cdots = \beta_{m} = 0, \qquad i = 1, \ldots, m.}
#'   }
#' }
#'
#' For each hypothesis, four asymptotic LM-type tests are reported:
#' \itemize{
#'   \item \eqn{\chi^2}-version LM test.
#'   \item F-version LM test.
#'   \item \eqn{\chi^2}-version HAC LM test (heteroskedasticity and autocorrelation consistent).
#'   \item F-version HAC LM test.
#' }
#'
#' \code{WCB_LinTest} additionally reports wild bootstrap (WB) and wild cluster bootstrap (WCB) p-values.
#' WB is robust to heteroskedasticity, while WCB is robust to both heteroskedasticity and within-individual
#' dependence (cluster dependence). The bootstrap routines can be computationally expensive; parallel execution
#' can be enabled via \code{parallel = TRUE}.
#'
#' Results are stored in the returned object (see \strong{Value}).
#'
#' @param use An object of class \code{"PSTR"} created by \code{\link{NewPSTR}}.
#' @param iB Integer. Number of bootstrap repetitions. Default is \code{100}.
#' @param parallel Logical. Whether to use parallel computation in bootstrap routines.
#' @param cpus Integer. Number of CPU cores to use when \code{parallel = TRUE}. Ignored otherwise.
#'
#' @return
#' Both functions return \code{use} invisibly, after adding the following components:
#' \describe{
#'   \item{\code{test}}{List. Asymptotic linearity test results for each transition variable and \eqn{m}.}
#'   \item{\code{sqtest}}{List. Asymptotic sequence test results for each transition variable and \eqn{m}.}
#'   \item{\code{wcb_test}}{List (only for \code{WCB_LinTest}). WB and WCB p-values for the linearity tests.}
#'   \item{\code{wcb_sqtest}}{List (only for \code{WCB_LinTest}). WB and WCB p-values for the sequence tests.}
#' }
#'
#' @seealso \code{\link{NewPSTR}}, \code{\link{EstPSTR}}, \code{\link{EvalTest}}.
#'
#' @examples
#' pstr <- NewPSTR(Hansen99, dep = "inva", indep = 4:20,
#'                indep_k = c("vala","debta","cfa","sales"),
#'                tvars = c("vala"), iT = 14)
#'
#' # R6 method interface
#' pstr$LinTest()
#'
#' # Wrapper interface (equivalent)
#' pstr <- LinTest(pstr)
#'
#' # Show results
#' print(pstr, mode = "tests")
#'
#' \donttest{
#' # Bootstrap tests (can be slow)
#' pstr$WCB_LinTest(iB = 200, parallel = TRUE, cpus = 2)
#' # or
#' pstr <- WCB_LinTest(use = pstr, iB = 200, parallel = TRUE, cpus = 2)
#'
#' print(pstr, mode = "tests")
#' }
#'
#' @name LinTest
NULL


PSTR$set("public", "LinTest", function(){
  iT = private$iT; iN = private$iN
  im = private$im
  
  # get the data here
  vY = private$vY; vYb = private$vYb
  mX = private$mX; mXb = private$mXb
  mK = private$mK
  
  self$test = list(); length(self$test) = ncol(private$mQ)
  self$sqtest = list(); length(self$sqtest) = ncol(private$mQ)
  
  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - tcrossprod(mD)/iT
  
  mX2 = mM %*% mX; invXX = chol2inv(chol(crossprod(mX2)))
  tmp = chol2inv(chol(crossprod(mXb))) %*% crossprod(mXb,vYb)
  vU = matrix(c(vY-mX%*%tmp), iT, iN)
  vU = c(t(t(vU)-apply(t(vU), 1, mean)))
  s2 = sum((vU-mean(vU))**2)/(iT*iN) # sigma^2
  
  coln = c(t(matrix(1:iN, iN, iT)))
  
  for(qter in 1:ncol(private$mQ)){
    vQ = private$mQ[,qter]
    
    self$test[[qter]] = list(); length(self$test[[qter]]) = im
    self$sqtest[[qter]] = list(); length(self$sqtest[[qter]]) = im
    
    mW = mK*vQ
    self$test[[qter]][[1]] = LMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)
    self$sqtest[[qter]][[1]] = self$test[[qter]][[1]]
    
    if(im>1) for(mter in 2:im){
      mXK = cbind(mX,mW); mX2K = mM %*% mXK; invXK = chol2inv(chol(crossprod(mX2K)))
      mXKb = NULL; for(nter in 1:iN){
        tmp = coln==nter; mXKb = rbind(mXKb, t(t(mW[tmp,])-apply(t(mW[tmp,]),1,mean)))
      }
      mXKb = cbind(mXb, mXKb)
      tmp = chol2inv(chol(crossprod(mXKb))) %*% crossprod(mXKb,vYb)
      vUK = matrix(c(vY-mXK%*%tmp), iT, iN)
      vUK = c(t(t(vUK)-apply(t(vUK), 1, mean)))
      s2K = sum((vUK-mean(vUK))**2)/(iT*iN) # sigma^2
      mWK = mK*(vQ**mter)
      
      self$sqtest[[qter]][[mter]] = LMTEST(iT=iT,iN=iN,vU=vUK,mX=mXK,mW=mWK,mM=mM,s2=s2K,mX2=mX2K,invXX=invXK)
      
      mW = cbind(mW, mWK)
      self$test[[qter]][[mter]] = LMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)
    }
  }
  
  cli::cli_alert_success("Done!")
  
  invisible(self)
})


#' @rdname LinTest
#' @export
LinTest <- function(use)
{
  if(!inherits(use, 'PSTR'))
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  
  use$LinTest()
  invisible(use)
}


PSTR$set("public", "WCB_LinTest", function(iB = 100, parallel = FALSE, cpus = 2) {
  
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
  
  iT <- private$iT
  iN <- private$iN
  im <- private$im
  
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
  
  # get the data here
  vY <- private$vY
  vYb <- private$vYb
  mX <- private$mX
  mXb <- private$mXb
  mK <- private$mK
  
  self$wcb_test <- vector("list", ncol(private$mQ))
  self$wcb_sqtest <- vector("list", ncol(private$mQ))
  
  beta <- chol2inv(chol(crossprod(mXb))) %*% crossprod(mXb, vYb)
  vU <- matrix(vY - mX %*% beta, iT, iN)
  mu <- apply(vU, 2, mean)
  vU <- c(t(t(vU) - mu))
  s2 <- sum((vU - mean(vU))^2) / (iT * iN) # sigma^2
  
  mD <- diag(1, iN) %x% rep(1, iT)
  mM <- diag(1, iN * iT) - tcrossprod(mD) / iT
  mX2 <- mM %*% mX
  invXX <- chol2inv(chol(crossprod(mX2)))
  eY <- mD %*% mu + mX %*% beta
  
  coln <- c(t(matrix(1:iN, iN, iT)))
  
  ftmp_wb <- function(bter) { # WB
    ve <- sample(c(1, -1), iT * iN, replace = TRUE) * vU
    my <- matrix(eY + ve, iT, iN)
    vyb <- c(t(t(my) - apply(my, 2, mean)))
    tmp <- chol2inv(chol(crossprod(mXb))) %*% crossprod(mXb, vyb)
    vu <- matrix(c(c(my) - mX %*% tmp), iT, iN)
    vu <- c(t(t(vu) - apply(vu, 2, mean)))
    ss <- sum((vu - mean(vu))^2) / (iT * iN)
    sLMTEST(iT = iT, iN = iN, vU = vu, mX = mX, mW = mW, mM = mM, s2 = ss, mX2 = mX2, invXX = invXX)
  }
  
  ftmp_wcb <- function(bter) { # WCB
    ve <- c(t(matrix(sample(c(1, -1), iN, replace = TRUE), iN, iT))) * vU
    my <- matrix(eY + ve, iT, iN)
    vyb <- c(t(t(my) - apply(my, 2, mean)))
    tmp <- chol2inv(chol(crossprod(mXb))) %*% crossprod(mXb, vyb)
    vu <- matrix(c(c(my) - mX %*% tmp), iT, iN)
    vu <- c(t(t(vu) - apply(vu, 2, mean)))
    ss <- sum((vu - mean(vu))^2) / (iT * iN)
    sLMTEST(iT = iT, iN = iN, vU = vu, mX = mX, mW = mW, mM = mM, s2 = ss, mX2 = mX2, invXX = invXX)
  }
  
  sqftmp_wb <- function(bter) { # WB
    ve <- sample(c(1, -1), iT * iN, replace = TRUE) * vUK
    my <- matrix(eYK + ve, iT, iN)
    vyb <- c(t(t(my) - apply(my, 2, mean)))
    tmp <- chol2inv(chol(crossprod(mXKb))) %*% crossprod(mXKb, vyb)
    vu <- matrix(c(c(my) - mXK %*% tmp), iT, iN)
    vu <- c(t(t(vu) - apply(vu, 2, mean)))
    ss <- sum((vu - mean(vu))^2) / (iT * iN)
    sLMTEST(iT = iT, iN = iN, vU = vu, mX = mXK, mW = mWK, mM = mM, s2 = ss, mX2 = mX2K, invXX = invXK)
  }
  
  sqftmp_wcb <- function(bter) { # WCB
    ve <- c(t(matrix(sample(c(1, -1), iN, replace = TRUE), iN, iT))) * vUK
    my <- matrix(eYK + ve, iT, iN)
    vyb <- c(t(t(my) - apply(my, 2, mean)))
    tmp <- chol2inv(chol(crossprod(mXKb))) %*% crossprod(mXKb, vyb)
    vu <- matrix(c(c(my) - mXK %*% tmp), iT, iN)
    vu <- c(t(t(vu) - apply(vu, 2, mean)))
    ss <- sum((vu - mean(vu))^2) / (iT * iN)
    sLMTEST(iT = iT, iN = iN, vU = vu, mX = mXK, mW = mWK, mM = mM, s2 = ss, mX2 = mX2K, invXX = invXK)
  }
  
  for (qter in 1:ncol(private$mQ)) {
    
    vQ <- private$mQ[, qter]
    
    mW <- mK * vQ
    LM <- sLMTEST(iT = iT, iN = iN, vU = vU, mX = mX, mW = mW, mM = mM, s2 = s2, mX2 = mX2, invXX = invXX)
    
    qLM1 <- .run_boot(ftmp_wb)
    qLM2 <- .run_boot(ftmp_wcb)
    
    rtmp <- c(LM, mean(LM <= qLM1), mean(LM <= qLM2))
    rrtmp <- rtmp
    
    if (im > 1) for (mter in 2:im) {
      
      mXK <- cbind(mX, mW)
      mX2K <- mM %*% mXK
      invXK <- chol2inv(chol(crossprod(mX2K)))
      
      mXKb <- NULL
      for (nter in 1:iN) {
        tmp <- coln == nter
        mXKb <- rbind(mXKb, t(t(mW[tmp, ]) - apply(t(mW[tmp, ]), 1, mean)))
      }
      mXKb <- cbind(mXb, mXKb)
      
      tmp <- chol2inv(chol(crossprod(mXKb))) %*% crossprod(mXKb, vYb)
      vUK <- matrix(c(vY - mXK %*% tmp), iT, iN)
      muK <- apply(vUK, 2, mean)
      vUK <- c(t(t(vUK) - muK))
      s2K <- sum((vUK - mean(vUK))^2) / (iT * iN)
      eYK <- mD %*% muK + mXK %*% tmp
      
      mWK <- mK * (vQ^mter)
      
      sqLM <- sLMTEST(iT = iT, iN = iN, vU = vUK, mX = mXK, mW = mWK, mM = mM, s2 = s2K, mX2 = mX2K, invXX = invXK)
      
      sqLM1 <- .run_boot(sqftmp_wb)
      sqLM2 <- .run_boot(sqftmp_wcb)
      
      rrtmp <- rbind(rrtmp, c(sqLM, mean(sqLM <= sqLM1), mean(sqLM <= sqLM2)))
      
      # update joint test with expanded W
      mW <- cbind(mW, mWK)
      LM <- sLMTEST(iT = iT, iN = iN, vU = vU, mX = mX, mW = mW, mM = mM, s2 = s2, mX2 = mX2, invXX = invXX)
      
      qLM1 <- .run_boot(ftmp_wb)
      qLM2 <- .run_boot(ftmp_wcb)
      
      rtmp <- rbind(rtmp, c(LM, mean(LM <= qLM1), mean(LM <= qLM2)))
    }
    
    self$wcb_test[[qter]] <- matrix(rtmp, nrow = im)
    self$wcb_sqtest[[qter]] <- matrix(rrtmp, nrow = im)
  }
  
  cli::cli_alert_success("Done!")
  invisible(self)
})


#' @rdname LinTest
#' @export
WCB_LinTest <- function(use, iB=100, parallel=FALSE, cpus=2)
{
  if(!inherits(use, 'PSTR'))
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  
  use$WCB_LinTest(iB, parallel, cpus)
  invisible(use)
}
