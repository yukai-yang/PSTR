#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## Sep 2017
#################################################################################


#' Conduct the evaluation tests.
#'
#' These functions conduct the evaluation tests against two alternatives: 1. the parameters are time-varying and 2. there is remaining nonlinearity (remaining heterogeneity).
#'
#' \code{EvalTest} implements the evaluation tests.
#'
#' \code{WCB_TVTest} implements the wild bootstrap (WB) and the wild cluster bootstrap (WCB) evaluation test of no time-varying parameters.
#'
#' \code{WCB_HETest} implements the wild bootstrap (WB) and the wild cluster bootstrap (WCB) evaluation test of no remaining nonlinearity (no remaining heterogeneity).
#'
#' The functions need the return value (an object of the class PSTR) from the \code{\link{EstPSTR}}. Note that the PSTR model should be estimated before conducting the evaluation tests. They copy the object, reuse its contents, especially the estimates, to produce the evaluation test results, and then return a new object of the class PSTR. The user can choose to save the return value to a new object or simply to overwrite the object returned from \code{\link{EstPSTR}}. See the example below.
#'
#' The functions conduct two kinds of evaluation tests.
#' The first kind of tests does the time-varying evaluation tests.
#' The second kind of tests does the no remaining nonlinearity (no remaining heterogeneity) evaluation tests based on the vector of a new transition variable that the user input in the arguments.
#'
#' The results of the evaluation tests include four kinds of tests
#' \itemize{
#'   \item \eqn{\chi^2}-version LM test: the LM test with asymptotically \eqn{\chi^2} distribution under the null hypothesis.
#'   \item F-version LM test: the LM test with asymptotically \eqn{F} distribution under the null hypothesis. The finite sample actual size is supposed to be improved.
#'   \item \eqn{\chi^2}-version HAC test: the HAC LM test with asymptotically \eqn{\chi^2} distribution under the null hypothesis, which is heteroskedasticity and autocorrelation consistent.
#'   \item F-version HAC test: the HAC LM test with asymptotically \eqn{F} distribution under the null hypothesis, which is heteroskedasticity and autocorrelation consistent. The finite sample actual size is supposed to be improved.
#' }
#'
#' The wild bootstrap (WB) evaluation tests are heteroskedasticity robust, while the wild cluster bootstrap (WCB) ones are both cluster-dependency and heteroskedasticity robust. Cluster-dependency implies that there can be dependency (autocorrelation) within individual, but no correlation across individuals. The WB and WCB tests may take quite a long time to run which depends on the model specification and the number of repetitions \code{iB}. It is strongly recommended to use super-computation server with many cores to run the code instead of a personal computer. The user may first try a small number of repetitions \code{iB} and estimate the time consumed for a larger number of \code{iB}.
#'
#' The functions never change the existing values in the input PSTR object. They add more values (attributes) into the input object and return.
#'
#' @param use an object of the class PSTR, created by \code{\link{EstPSTR}} function.
#' @param type a character vector specifying the types of the evaluation tests to be conducted. The value can be taken either or both of \code{"time-varying"} \code{"heterogeneity"}. By default, do both.
#' @param vq a vector of a new transition variable for the no remaining nonlinearity test.
#' @param iB specify the number of repetitions in the bootstrap procedure. By default, it is 100.
#' @param parallel a boolean value showing if the parallel computation is applied.
#' @param cpus number of cores used in the parallel computation. The value will be ignored if \code{parallel=F}.
#'
#' @return a new object of the class PSTR containing the results from the evaluation tests.
#'
#' The return object from \code{EvalTest} contains the following new components:
#' \item{tv}{a list of the time-varying evaluation tests. The length of the list is the maximal number of switches. Each element of the list corresponds to the time-varying evaluation test results based on different number of switches.}
#' \item{ht}{a list of the no remaining nonlinearity (no remaining heterogeneity) evaluation tests. The length of the list is the maximal number of switches. Each element of the list corresponds to the time-varying evaluation test results based on different number of switches. The input vector of a new transition variable is used to compute the tests.}
#'
#' The return object from \code{WCB_TVTest} contains the following new component:
#' \item{wcb_tv}{a matrix containing the results from the WB and WCB time-varying tests. Each row of the matrix contains the p-value of the WB and WCB tests.}
#'
#' The return object from \code{WCB_HETest} contains the following new component:
#' \item{wcb_ht}{a matrix containing the results from the WB and WCB no remaining nonlinearity (heterogeneity) tests. Each row of the matrix contains the p-value of the WB and WCB tests.}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{NewPSTR}}, \code{\link{LinTest}}, \code{\link{WCB_LinTest}} and \code{\link{EstPSTR}}
#' @keywords evaluation
#'
#' @examples
#' \donttest{
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala'), iT=14) # create a new PSTR object
#'
#' # Estimate the model first
#' pstr = EstPSTR(use=pstr, im=1, iq=1, useDelta=TRUE, par=c(.63,0), method='CG')
#'
#' # Then you can evaluate the model
#' pstr = EvalTest(use=pstr, vq=pstr$mQ[,1])
#'
#' print(pstr, "eval")
#'
#' # You can do the wild bootstrap and wild cluster bootstrap
#'
#' library(snowfall)
#'
#' pstr = WCB_TVTest(use=pstr, iB=4, parallel=TRUE, cpus=2)
#'
#' # pstr$mQ[,1] is the transition variable stored in the object
#' # You can also try other transition variables.
#' pstr = WCB_HETest(use=pstr, vq=pstr$mQ[,1], iB=4, parallel=TRUE, cpus=2)
#'
#' print(pstr, "eval")
#'
#' # Don't forget to change the values of iB and cpus during experiments.
#' }
#'
#' @name EvalTest
NULL



#' @rdname EvalTest
#' @export
EvalTest <- function(use, type=c("time-varying","heterogeneity"), vq=NULL)
{
  if(!inherits(use, 'PSTR'))
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  if(is.null(use$iq))
    stop(simpleError("Estimate the PSTR model first!"))
  ret = use
  im = use$im
  
  mD = diag(1,use$iN) %x% rep(1,use$iT)
  mM = diag(1, use$iN*use$iT) - tcrossprod(mD)/use$iT
  
  tmp = c(use$mK %*% use$beta[(ncol(use$mX)+1):length(use$beta)])
  tmp = use$mD * tmp ## pp.14
  mV = cbind(use$mXX, tmp)
  mV2 = mM %*% mV
  invVV = chol2inv(chol(crossprod(mV2)))
  
  if(length(grep("time-varying",type))>0){
    ret$tv = list(); length(ret$tv) = im
    
    vt = 1:use$iT/use$iT
    
    mW = NULL
    for(mter in 1:im){
      mW = cbind(mW, use$mXX*(vt**mter))
      ret$tv[[mter]] = LMTEST(iT=use$iT,iN=use$iN,vU=use$vU,mX=mV,mW=mW,mM=mM,s2=use$s2,mX2=mV2,invXX=invVV)
    }
  }
  
  if(length(grep("heterogeneity",type))>0){
    ret$ht = list(); length(ret$ht) = im
    
    mW = NULL
    for(mter in 1:im){
      mW = cbind(mW, use$mXX*(vq**mter))
      ret$ht[[mter]] = LMTEST(iT=use$iT,iN=use$iN,vU=use$vU,mX=mV,mW=mW,mM=mM,s2=use$s2,mX2=mV2,invXX=invVV)
    }
  }
  
  return(ret)
}



#' @rdname EvalTest
#' @export
WCB_TVTest <- function(use, iB=100, parallel=F, cpus=4)
{
  if(!inherits(use, 'PSTR'))
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  if(is.null(use$iq))
    stop(simpleError("Estimate the PSTR model first!"))
  ret = use; ruse = use
  im = use$im
  
  iT = use$iT; iN = use$iN
  vU = use$vU; eY = use$vY - vU
  
  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - tcrossprod(mD)/iT
  
  tmp = c(use$mK %*% use$beta[(ncol(use$mX)+1):length(use$beta)])
  tmp = use$mD * tmp
  mV = cbind(use$mXX, tmp)
  mV2 = mM %*% mV
  invVV = chol2inv(chol(crossprod(mV2)))
  
  ftmp_wb <- function(bter){# WB
    ve1 = sample(c(1,-1),iT*iN,replace=T)*vU
    ruse$vY = eY + ve1
    EST = EstPSTR(use=ruse,im=1,iq=ruse$iq,par=c(use$delta,use$c),useDelta=T,vLower=1,vUpper=1)
    vu1 = EST$vU; ss1 = EST$s2 # sigma^2
    tmp = c(EST$mK%*%EST$beta[(ncol(EST$mX)+1):length(EST$beta)])
    tmp = EST$mD * tmp
    mV11 = cbind(EST$mXX, tmp)
    mV12 = mM %*% mV11
    invVV1 = chol2inv(chol(t(mV12)%*%mV12))
    return(sLMTEST(iT=iT,iN=iN,vU=vu1,mX=mV11,mW=mW,mM=mM,s2=ss1,mX2=mV12,invXX=invVV1))
  }
  
  ftmp_wcb <- function(bter){# WCB
    ve2 = c(t(matrix(sample(c(1,-1),iN,replace=T),iN,iT)))*vU
    ruse$vY = eY + ve2
    EST = EstPSTR(use=ruse,im=1,iq=ruse$iq,par=c(use$delta,use$c),useDelta=T,vLower=1,vUpper=1)
    vu2 = EST$vU; ss2 = EST$s2 # sigma^2
    tmp = c(EST$mK%*%EST$beta[(ncol(EST$mX)+1):length(EST$beta)])
    tmp = EST$mD * tmp
    mV21 = cbind(EST$mXX, tmp)
    mV22 = mM %*% mV21
    invVV2 = chol2inv(chol(t(mV22)%*%mV22))
    return(sLMTEST(iT=iT,iN=iN,vU=vu2,mX=mV21,mW=mW,mM=mM,s2=ss2,mX2=mV22,invXX=invVV2))
  }
  
  ret$wcb_tv = NULL
  vt = 1:iT/iT
  mW = NULL
  
  for(mter in 1:im){
    mW = cbind(mW, use$mXX*(vt**mter))
    LM = sLMTEST(iT=iT,iN=iN,vU=vU,mX=mV,mW=mW,mM=mM,s2=use$s2,mX2=mV2,invXX=invVV)
    
    sfInit(parallel=parallel,cpus=cpus)
    qLM1 = sfSapply(1:iB,ftmp_wb)
    qLM2 = sfSapply(1:iB,ftmp_wcb)
    sfStop()
    
    ret$wcb_tv = rbind(ret$wcb_tv,c(LM, mean(LM<=qLM1), mean(LM<=qLM2)))
  }
  
  ret$wcb_tv = matrix(ret$wcb_tv, nrow=im)
  
  return(ret)
}



#' @rdname EvalTest
#' @export
WCB_HETest <- function(use, vq, iB=100, parallel=F, cpus=4)
{
  if(!inherits(use, 'PSTR'))
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  if(is.null(use$iq))
    stop(simpleError("Estimate the PSTR model first!"))
  ret = use; ruse = use
  im = use$im
  
  iT = use$iT; iN = use$iN
  vU = use$vU; eY = use$vY - vU
  
  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - tcrossprod(mD)/iT
  
  tmp = c(use$mK %*% use$beta[(ncol(use$mX)+1):length(use$beta)])
  tmp = use$mD * tmp
  mV = cbind(use$mXX, tmp)
  mV2 = mM %*% mV
  invVV = chol2inv(chol(crossprod(mV2)))
  
  ftmp_wb <- function(bter){# WB
    ve1 = sample(c(1,-1),iT*iN,replace=T)*vU
    ruse$vY = eY + ve1
    EST = EstPSTR(use=ruse,im=1,iq=ruse$iq,par=c(use$delta,use$c),useDelta=T,vLower=1,vUpper=1)
    vu1 = EST$vU; ss1 = EST$s2 # sigma^2
    tmp = c(EST$mK%*%EST$beta[(ncol(EST$mX)+1):length(EST$beta)])
    tmp = EST$mD * tmp
    mV11 = cbind(EST$mXX, tmp)
    mV12 = mM %*% mV11
    invVV1 = chol2inv(chol(t(mV12)%*%mV12))
    return(sLMTEST(iT=iT,iN=iN,vU=vu1,mX=mV11,mW=mW,mM=mM,s2=ss1,mX2=mV12,invXX=invVV1))
  }
  
  ftmp_wcb <- function(bter){# WCB
    ve2 = c(t(matrix(sample(c(1,-1),iN,replace=T),iN,iT)))*vU
    ruse$vY = eY + ve2
    EST = EstPSTR(use=ruse,im=1,iq=ruse$iq,par=c(use$delta,use$c),useDelta=T,vLower=1,vUpper=1)
    vu2 = EST$vU; ss2 = EST$s2 # sigma^2
    tmp = c(EST$mK%*%EST$beta[(ncol(EST$mX)+1):length(EST$beta)])
    tmp = EST$mD * tmp
    mV21 = cbind(EST$mXX, tmp)
    mV22 = mM %*% mV21
    invVV2 = chol2inv(chol(t(mV22)%*%mV22))
    return(sLMTEST(iT=iT,iN=iN,vU=vu2,mX=mV21,mW=mW,mM=mM,s2=ss2,mX2=mV22,invXX=invVV2))
  }
  
  ret$wcb_ht = NULL
  mW = NULL
  
  for(mter in 1:im){
    mW = cbind(mW, use$mXX*(vq**mter))
    LM = sLMTEST(iT=iT,iN=iN,vU=vU,mX=mV,mW=mW,mM=mM,s2=use$s2,mX2=mV2,invXX=invVV)
    
    sfInit(parallel=parallel,cpus=cpus)
    qLM1 = sfSapply(1:iB,ftmp_wb)
    qLM2 = sfSapply(1:iB,ftmp_wcb)
    sfStop()
    
    ret$wcb_ht = rbind(ret$wcb_ht,c(LM, mean(LM<=qLM1), mean(LM<=qLM2)))
  }
  
  ret$wcb_ht = matrix(ret$wcb_ht, nrow=im)
  
  return(ret)
}