#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## Sep 2017
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
  if(class(invS1)=='try-error'){
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
  if(class(invS2)=='try-error'){
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
  if(class(invS1)=='try-error'){
    ttmp = svd(S1); invS1 = ttmp$u %*% diag(1/ttmp$d, length(ttmp$d)) %*% t(ttmp$u)
  }
  
  vW2U = crossprod(mW2, vU)
  LM1_X = c(t(vW2U) %*% invS1 %*% vW2U)
  
  return(LM1_X)
}


#' Conduct the linearity tests.
#'
#' These functions conduct the linearity tests against the alternative of a logistic smooth transition nonlinear component.
#'
#' \code{LinTest} implements the linearity tests.
#'
#' \code{WCB_LinTest} implements the wild bootstrap (WB) and the wild cluster bootstrap (WCB) linearity tests.
#'
#' The functions need the return value (an object of the class PSTR) from the \code{\link{NewPSTR}}. They copy the object, reuse its contents to produce the linearity test results, and then return a new object of the class PSTR. The user can choose to save the return value to a new object or simply to overwrite the object returned from \code{NewPSTR}. See the example below.
#'
#' The functions conduct two kinds of linearity tests.
#'
#' The first kind of tests does the linearity tests based on each potential transition variable specified in the argument \code{tvars} when the user calls the \code{\link{NewPSTR}} function. For each potential transition variable, the function conducts linearity tests for numbers of switches from 1 up to \code{im}. The linearity tests has the null hypothesis
#' \deqn{H_0^i: \beta_{i} = \beta_{i-1} = \beta_{i-2} = ... = \beta_{1} = 0}
#' for \eqn{i = 1, ..., m}, where \eqn{m} is the maximal number of switches \code{im}.
#'
#' The second kind does the linearity tests for selecting the number of switches based on each potential transition variable. The linearity tests for selecting the number of switches has the null hypothesis
#' \deqn{H_0^i: \beta_{i} = 0 | \beta_{i+1} = \beta_{i+2} = ... = \beta_{m} = 0}
#' for \eqn{i = 1, ..., m}, where \eqn{m} is the maximal number of switches \code{im}.
#'
#' The results of the linearity tests include four kinds of tests
#' \itemize{
#'   \item \eqn{\chi^2}-version Linearity test: the linearity LM test with asymptotically \eqn{\chi^2} distribution under the null hypothesis of linearity.
#'   \item F-version Linearity test: the linearity LM test with asymptotically \eqn{F} distribution under the null hypothesis of linearity. The finite sample actual size is supposed to be improved.
#'   \item \eqn{\chi^2}-version HAC Linearity test: the linearity LM test with asymptotically \eqn{\chi^2} distribution under the null hypothesis of linearity, which is heteroskedasticity and autocorrelation consistent.
#'   \item F-version HAC Linearity test: the linearity LM test with asymptotically \eqn{F} distribution under the null hypothesis of linearity, which is heteroskedasticity and autocorrelation consistent. The finite sample actual size is supposed to be improved.
#' }
#'
#' The wild bootstrap (WB) tests are heteroskedasticity robust, while the wild cluster bootstrap (WCB) ones are both cluster-dependency and heteroskedasticity robust. Cluster-dependency implies that there can be dependency (autocorrelation) within individual, but no correlation across individuals. The WB and WCB tests may take quite a long time to run which depends on the model specification and the number of repetitions \code{iB}. It is strongly recommended to use super-computation server with many cores to run the code instead of a personal computer. The user may first try a small number of repetitions \code{iB} and estimate the time consumed for a larger number of \code{iB}.
#'
#' The two functions never change the existing values in the input PSTR object. They add more values (attributes) into the input object and return.
#'
#' @param use an object of the class PSTR, created by \code{\link{NewPSTR}} function.
#' @param iB specify the number of repetitions in the bootstrap procedure. By default, it is 100.
#' @param parallel a boolean value showing if the parallel computation is applied.
#' @param cpus number of cores used in the parallel computation. The value will be ignored if \code{parallel=F}.
#'
#' @return a new object of the class PSTR containing the results from the linearity tests.
#'
#' The object is a list containing the components made in \code{\link{NewPSTR}} and the following new components:
#' \item{test}{a list of the linearity test results. The length is the number of potential transition variables specified when creating the object of the class PSTR by calling \code{NewPSTR}. See argument \code{tvars} in \code{\link{NewPSTR}}. Each element of the list corresponds to the linearity test results based on the corresponding transition variable, and the element is also a list whose elements correspond to different numbers of switches.}
#' \item{sqtest}{a list of the linearity test results for selecting number of switches. It has the same length as \code{test}. Each element of the list corresponds to the linearity test results based on the corresponding transition variable, and the element is also a list whose elements correspond to different numbers of switches.}
#' \item{wcb_test}{a list of the linearity test results. The length is the number of potential transition variables specified when creating the object of the class PSTR by calling \code{NewPSTR}. See argument \code{tvars} in \code{\link{NewPSTR}}. Each element of the list is a matrix containing the linearity test results (p-values) based on the corresponding transition variable. The rows are different numbers of switches, and two columns from WB to WCB.}
#' \item{wcb_sqtest}{a list of the linearity test results for selecting number of switches. It has the same length as \code{test}. Each element of the list is a matrix containing the linearity test results based on the corresponding transition variable. The rows are different numbers of switches, and two columns from WB to WCB.}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{NewPSTR}}
#' @keywords specification
#'
#' @examples
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala'), iT=14) # create a new PSTR object
#'
#' pstr = LinTest(pstr)
#'
#' print(pstr, "tests")
#'
#' \donttest{
#' # Don't forget to attach the package for the parallel computation.
#' library(snowfall)
#'
#' # you should not run this on personal computer!
#' # pstr = WCB_LinTest(use=pstr, iB=5000, parallel=TRUE, cpus=50)
#'
#' # a light version for checking on your personal computer.
#' pstr = WCB_LinTest(use=pstr, iB=4, parallel=TRUE, cpus=2)
#'
#' print(pstr, "tests")
#' }
#'
#' @name LinTest
NULL


#' @rdname LinTest
#' @export
LinTest <- function(use)
{
  if(class(use)!="PSTR")
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  ret = use
  iT = use$iT; iN = use$iN
  im = use$im
  
  # get the data here
  vY = use$vY; vYb = use$vYb
  mX = use$mX; mXb = use$mXb
  mK = use$mK
  
  ret$test = list(); length(ret$test) = ncol(use$mQ)
  ret$sqtest = list(); length(ret$sqtest) = ncol(use$mQ)
  
  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - tcrossprod(mD)/iT
  
  mX2 = mM %*% mX; invXX = chol2inv(chol(crossprod(mX2)))
  tmp = chol2inv(chol(crossprod(mXb))) %*% crossprod(mXb,vYb)
  vU = matrix(c(vY-mX%*%tmp), iT, iN)
  vU = c(t(t(vU)-apply(t(vU), 1, mean)))
  s2 = sum((vU-mean(vU))**2)/(iT*iN) # sigma^2
  
  coln = c(t(matrix(1:iN, iN, iT)))
  
  for(qter in 1:ncol(use$mQ)){
    vQ = use$mQ[,qter]
    
    ret$test[[qter]] = list(); length(ret$test[[qter]]) = im
    ret$sqtest[[qter]] = list(); length(ret$sqtest[[qter]]) = im
    
    mW = mK*vQ
    ret$test[[qter]][[1]] = LMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)
    ret$sqtest[[qter]][[1]] = ret$test[[qter]][[1]]
    
    
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
      
      ret$sqtest[[qter]][[mter]] = LMTEST(iT=iT,iN=iN,vU=vUK,mX=mXK,mW=mWK,mM=mM,s2=s2K,mX2=mX2K,invXX=invXK)
      
      mW = cbind(mW, mWK)
      ret$test[[qter]][[mter]] = LMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)
    }
  }
  
  return(ret)
  
}


#' @rdname LinTest
#' @export
WCB_LinTest <- function(use, iB=100, parallel=F, cpus=4)
{
  if(class(use)!="PSTR")
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  ret = use
  iT = use$iT; iN = use$iN
  im = use$im
  
  # get the data here
  vY = use$vY; vYb = use$vYb
  mX = use$mX; mXb = use$mXb
  mK = use$mK
  
  ret$wcb_test = list(); length(ret$wcb_test) = ncol(use$mQ)
  ret$wcb_sqtest = list(); length(ret$wcb_sqtest) = ncol(use$mQ)
  
  beta = chol2inv(chol(crossprod(mXb))) %*% crossprod(mXb,vYb)
  vU = matrix(vY-mX%*%beta, iT, iN)
  mu = apply(vU, 2, mean)
  vU = c(t(t(vU) - mu))
  s2 = sum((vU-mean(vU))**2)/(iT*iN) # sigma^2
  
  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - tcrossprod(mD)/iT
  mX2 = mM %*% mX
  invXX = chol2inv(chol(crossprod(mX2)))
  eY = mD%*%mu + mX%*%beta
  
  coln = c(t(matrix(1:iN, iN, iT)))
  
  ftmp_wb <- function(bter){# WB
    ve = sample(c(1,-1),iT*iN,replace=T)*vU
    my = matrix(eY + ve, iT, iN)
    vyb = c(t(t(my) - apply(my, 2, mean)))
    tmp = chol2inv(chol(crossprod(mXb))) %*% crossprod(mXb,vyb)
    vu = matrix(c(c(my)-mX%*%tmp), iT, iN)
    vu = c(t(t(vu)-apply(vu, 2, mean)))
    ss = sum((vu-mean(vu))**2)/(iT*iN) # sigma^2
    return(sLMTEST(iT=iT,iN=iN,vU=vu,mX=mX,mW=mW,mM=mM,s2=ss,mX2=mX2,invXX=invXX))
  }
  
  ftmp_wcb <- function(bter){# WCB
    ve = c(t(matrix(sample(c(1,-1),iN,replace=T), iN, iT)))*vU
    my = matrix(eY + ve, iT, iN)
    vyb = c(t(t(my) - apply(my, 2, mean)))
    tmp = chol2inv(chol(crossprod(mXb))) %*% crossprod(mXb,vyb)
    vu = matrix(c(c(my)-mX%*%tmp), iT, iN)
    vu = c(t(t(vu)-apply(vu, 2, mean)))
    ss = sum((vu-mean(vu))**2)/(iT*iN) # sigma^2
    return(sLMTEST(iT=iT,iN=iN,vU=vu,mX=mX,mW=mW,mM=mM,s2=ss,mX2=mX2,invXX=invXX))
  }
  
  sqftmp_wb <- function(bter){# WB
    ve = sample(c(1,-1),iT*iN,replace=T)*vUK
    my = matrix(eYK + ve, iT, iN)
    vyb = c(t(t(my) - apply(my, 2, mean)))
    tmp = chol2inv(chol(crossprod(mXKb))) %*% crossprod(mXKb,vyb)
    vu = matrix(c(c(my)-mXK%*%tmp), iT, iN)
    vu = c(t(t(vu)-apply(vu, 2, mean)))
    ss = sum((vu-mean(vu))**2)/(iT*iN) # sigma^2
    return(sLMTEST(iT=iT,iN=iN,vU=vu,mX=mXK,mW=mWK,mM=mM,s2=ss,mX2=mX2K,invXX=invXK))
  }
  
  sqftmp_wcb <- function(bter){# WCB
    ve = c(t(matrix(sample(c(1,-1),iN,replace=T), iN, iT)))*vUK
    my = matrix(eYK + ve, iT, iN)
    vyb = c(t(t(my) - apply(my, 2, mean)))
    tmp = chol2inv(chol(crossprod(mXKb))) %*% crossprod(mXKb,vyb)
    vu = matrix(c(c(my)-mXK%*%tmp), iT, iN)
    vu = c(t(t(vu)-apply(vu, 2, mean)))
    ss = sum((vu-mean(vu))**2)/(iT*iN) # sigma^2
    return(sLMTEST(iT=iT,iN=iN,vU=vu,mX=mXK,mW=mWK,mM=mM,s2=ss,mX2=mX2K,invXX=invXK))
  }
  
  for(qter in 1:ncol(use$mQ)){
    vQ = use$mQ[,qter]
    
    mW = mK*vQ
    LM = sLMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)
    
    sfInit(parallel=parallel,cpus=cpus)
    #sfExport(list=c('sLMTEST'))
    qLM1 = sfSapply(1:iB,ftmp_wb)
    qLM2 = sfSapply(1:iB,ftmp_wcb)
    sfStop()
    
    rtmp = c(LM, mean(LM<=qLM1), mean(LM<=qLM2))
    rrtmp = rtmp
    
    if(im>1) for(mter in 2:im){
      
      mXK = cbind(mX,mW); mX2K = mM %*% mXK; invXK = chol2inv(chol(crossprod(mX2K)))
      mXKb = NULL; for(nter in 1:iN){
        tmp = coln==nter; mXKb = rbind(mXKb, t(t(mW[tmp,])-apply(t(mW[tmp,]),1,mean)))
      }
      mXKb = cbind(mXb, mXKb)
      tmp = chol2inv(chol(crossprod(mXKb))) %*% crossprod(mXKb,vYb)
      vUK = matrix(c(vY-mXK%*%tmp), iT, iN)
      muK = apply(vUK, 2, mean)
      vUK = c(t(t(vUK) - muK))
      s2K = sum((vUK-mean(vUK))**2)/(iT*iN) # sigma^2
      eYK = mD%*%muK + mXK%*%tmp
      
      mWK = mK*(vQ**mter)
      
      sqLM = sLMTEST(iT=iT,iN=iN,vU=vUK,mX=mXK,mW=mWK,mM=mM,s2=s2K,mX2=mX2K,invXX=invXK)
      
      sfInit(parallel=parallel,cpus=cpus)
      #sfExport(list=c('sLMTEST'))
      sqLM1 = sfSapply(1:iB,sqftmp_wb)
      sqLM2 = sfSapply(1:iB,sqftmp_wcb)
      sfStop()
      
      rrtmp = rbind( rrtmp, c(sqLM, mean(sqLM<=sqLM1), mean(sqLM<=sqLM2)) )
      
      ##
      mW = cbind(mW, mWK)
      LM = sLMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)
      
      sfInit(parallel=parallel,cpus=cpus)
      #sfExport(list=c('sLMTEST'))
      qLM1 = sfSapply(1:iB,ftmp_wb)
      qLM2 = sfSapply(1:iB,ftmp_wcb)
      sfStop()
      
      rtmp = rbind( rtmp, c(LM, mean(LM<=qLM1), mean(LM<=qLM2)) )
    }
    
    ret$wcb_test[[qter]] = matrix(rtmp, nrow=im)
    ret$wcb_sqtest[[qter]] = matrix(rrtmp, nrow=im)
    
  }
  
  return(ret)
}