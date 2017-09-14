#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## Aug 2017
#################################################################################


#' PSTR: A package implementing the Panel Smooth Transition Regression (PSTR) modelling.
#'
#' The package implements the Panel Smooth Transition Regression (PSTR) modelling.
#'
#' The modelling procedure consists of three stages: Specification, Estimation and Evaluation.
#' The package offers tools helping the package users to conduct model specification tests,
#' to do PSTR model estimation, and to do model evaluation.
#'
#' The cluster-dependency and heteroskedasticity-consistent tests are implemented in the package.
#'
#' The wild bootstrap and cluster wild bootstrap tests are also implemented.
#'
#' Parallel computation (as an option) is implemented in some functions, especially the bootstrap tests.
#' Therefore, the package suits tasks running many cores on super-computation servers.
#'
#' The Panel Smooth Transition Regression (PSTR) model is defined to be
#' \deqn{y_{it} = \mu_i + \beta_0' x_{it} + \beta_1' z_{it} g_{it} + u_{it}}
#' where \eqn{g_{it}} is the transition function taking the logistic form with the transition variable for individual \eqn{i}, \eqn{x_{it}} contains the explanatory variables in the linear part, and \eqn{z_{it}} contains the explanatory variables in the nonlinear part, and they can be different.
#'
#' @section Author and Maintainer:
#' Yukai Yang
#'
#' Department of Statistics, Uppsala University
#'
#' \email{yukai.yang@@statistik.uu.se}
#'
#' @section References:
#' González, A., Teräsvirta, T., van Dijk, D. and Yang, Y. (2017) \emph{Panel Smooth Transition Regression Models}
#'
#' @section Function for Initialization:
#' \code{\link{NewPSTR}} initialize the modelling by creating an object of the class PSTR.
#'
#' @section Functions for Model Specification:
#' \code{\link{LinTest}} implements the linearity tests.
#'
#' \code{\link{WCB_LinTest}} implements the wild bootstrap (WB) and the wild cluster bootstrap (WCB) linearity tests.
#'
#' @section Function for Model Estimation:
#' \code{\link{EstPSTR}} implements the estimation of the PSTR model.
#'
#' @section Functions for Model Evaluation:
#' \code{\link{EvalTest}} implements the evaluation tests.
#'
#' \code{\link{WCB_TVTest}} implements the wild bootstrap (WB) and the wild cluster bootstrap (WCB) evaluation test of no time-varying parameters.
#'
#' \code{\link{WCB_HETest}} implements the wild bootstrap (WB) and the wild cluster bootstrap (WCB) evaluation test of no remaining nonlinearity (no remaining heterogeneity).
#'
#' @section Other Functions:
#' \code{\link{version}} shows the version number and some information of the package.
#'
#' \code{\link{print.PSTR}} prints the object of the class PSTR.
#'
#' @docType package
#' @name PSTR
NULL


#' @importFrom stats optim pchisq pf
NULL

#' @import tibble
NULL

#' @importFrom snowfall sfInit sfExport sfSapply sfStop
NULL


#' Create an object of the class PSTR.
#'
#' Create an object of the S3 class PSTR for later usage. This function should be run prior to the other functions in the package. It will return an object which you will use as an input for the other functions. It builds up the basic settings for the Panel Smooth Transition Regression (PSTR) Modelling.
#'
#' Potential transition variables in \code{tvars} will be tested one by one in, for example, \code{LinTest} function.
#'
#' There is no need to specify the number of individuals,  as it will be obtained automatically inside the function given the number of rows and the sample size \code{iT}.
#'
#' \code{NA}s in \code{data} are removed automatically inside the function.
#'
#' @param data a tibble of data. The number of rows of \code{data} must be the sample size \code{iT} times individuals number N.
#' @param dep column number or name of the dependent variable.
#' @param indep a vector of column numbers of names of the independent variables.
#' @param indep_k a vector of column numbers of names of the independent variables in the nonlinear part. If \code{indep_k} is not given (\code{= NULL}), the nonlinear part will be the same as the linear part.
#' @param tvars a vector of column numbers or names of the potential transition variables to be tested.
#' @param im maximal number of switches in the transition function used in the linearity evaluation tests, by default \code{im=1}.
#' @param iT sample size.
#'
#' @return An object of the class PSTR for later usage.
#'
#' The object is a list containing the following components:
#' \item{iT}{the time length of the panel}
#' \item{iN}{the number of individuals}
#' \item{vY}{the vector of the dependent variable}
#' \item{mX}{the matrix of the explanatory variables in the linear part}
#' \item{mK}{the matrix of the explanatory variables in the nonlinear part}
#' \item{mQ}{the matrix of the potential transition variables}
#' \item{im}{the maximal number of switches used in the linearity test}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{LinTest}}
#' @keywords initialization
#'
#' @examples
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala','debta'), iT=14)
#'
#' pstr
#'
#' print(pstr,"summary")
#'
#' @export
NewPSTR <- function(data, dep, indep, indep_k=NULL, tvars, im=1, iT)
{
  ret = list(); class(ret) = "PSTR"
  if(!is_tibble(data)) stop(simpleError("data should be a tibble!"))

  ret$iT = iT
  iNN = dim(data)[1]/iT
  coln = c(t(matrix(1:iNN, iNN, iT)))

  if(length(dep)>1) stop(simpleError("Only one dependent variable!"))
  vY = data[,dep]; ret$vY_name =  names(data[,dep])
  mX = data[,indep]; ret$mX_name = names(data[,indep])
  if(is.null(indep_k)){
    mK = mX; ret$mK_name = ret$mX_name
  }else{
    mK = data[,indep_k]; ret$mK_name = names(data[,indep_k])
  }
  mQ = data[,tvars]; ret$mQ_name = names(data[,tvars])

  ## remove the NAs
  tmp = is.na(vY) | apply(is.na(mX), 1, any) | apply(is.na(mK), 1, any) | apply(is.na(mQ), 1, any)
  tmp = is.na(match(coln,coln[tmp]))
  vY = vY[tmp,]; mX = mX[tmp,]; mK = mK[tmp,]; mQ = mQ[tmp,]
  ret$vY = c(as.matrix(vY))
  ret$mX = as.matrix(mX)
  ret$mK = as.matrix(mK)
  ret$mQ = as.matrix(mQ)
  ret$im = im

  iN = sum(tmp)/iT; ret$iN = iN

  coln = c(t(matrix(1:iN, iN, iT)))
  vYb = NULL; mXb = NULL
  for(nter in 1:iN){
    tmp = coln==nter
    vYb = c(vYb, ret$vY[tmp] - mean(ret$vY[tmp]))
    mXb = rbind(mXb, t(t(mX[tmp,])-apply(t(mX[tmp,]),1,mean)))
  }
  ret$vYb = vYb; ret$mXb = mXb

  return(ret)
}



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
    ttmp = svd(S1); invS1 = ttmp$u %*% diag(1/ttmp$d) %*% t(ttmp$u)
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
    ttmp = svd(S2); invS2 = ttmp$u %*% diag(1/ttmp$d) %*% t(ttmp$u)
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
    ttmp = svd(S1); invS1 = ttmp$u %*% diag(1/ttmp$d) %*% t(ttmp$u)
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
#' # pstr = WCB_LinTest(use=pstr, iB=5000, parallel=T, cpus=50)
#'
#' # a light version for checking on your personal computer.
#' pstr = WCB_LinTest(use=pstr, iB=4, parallel=T, cpus=2)
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
  tmp2 = matrix(vs, length(vs), length(cc))
  tmp2 = t(tmp2) - cc # s - c
  tmp3 = apply(tmp2,2,prod) # prod all

  de1 = tmp1 * tmp3

  ftmp <- function(iter){
    tmp = tmp2; tmp[iter,] = 1
    return(apply(tmp,2,prod))
  }
  tmp4 = sapply(1:length(cc),ftmp) # prod without k

  de1 = cbind(de1, - tmp1 * tmp4 * gamma) # columns are the parameters

  de2 = de1[,1] * (1-2*vg) * tmp3 # d^2 g / d gamma^2

  # d^2 g / d gamma d c
  de2 = cbind(de2, 2*(1-vg) * de1[,2:ncol(de1)] * tmp3 + tmp1 * tmp3 * gamma * tmp4 - tmp1 * tmp4)

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



#' Estimate the PSTR model.
#'
#' This function implements the estimation of the \code{\link{PSTR}} model.
#'
#' The function needs the return value (an object of the class PSTR) from the \code{\link{NewPSTR}}. It copies the object, reuses its contents to estimate the correspdonding PSTR model, and then returns a new object of the class PSTR containing the results from the estimation. The user can choose to save the return value to a new object or simply to overwrite the object returned from \code{NewPSTR}.
#'
#' The PSTR model to be estimated takes the logistic form in nonlinearity. Remember the \eqn{g} function in the model. It takes the form
#' \deqn{g(q_{it} ; \gamma, c) = \left( 1 + \exp \left( - \gamma \prod_{j=1}^{m} (q_{it} - c_j) \right) \right)^{-1}}
#' with \eqn{\gamma > 0} and \eqn{c_1 < c_2 < ... < c_m}. \eqn{\gamma} can be reparametrized as \eqn{\gamma = \exp{\delta}} where \eqn{\delta} is a real number.
#'
#' The user should have obtained the information about which transition variable (\eqn{q_{it}}) to use (from \code{\link{LinTest}} and/or \code{\link{WCB_LinTest}}) in estimation before running the function to estimate the model.
#'
#' The estimation function never change the existing values in the input PSTR object. It adds more values (attributes) into the input object and return.
#'
#' @param use an object of the class PSTR, created by \code{\link{NewPSTR}} function.
#' @param im specifies the number of switches in the transtion function. The default value is 1.
#' @param iq a column number (in \code{mQ}) or variable name specifying the transition variable to use.
#' @param par initial values for the parameters (\eqn{\delta} and \eqn{c}) to be optimized over. It is a vector of length \code{im}+1, where \code{im} is the number of switches.
#' @param vLower a vector or number of the lower offsets determining the lower bounds of the parameters. The lower bounds of the parameters are \code{par - vLower}.
#' @param vUpper a vector or number of the upper offsets determining the upper bounds of the parameters. The upper bounds of the parameters are \code{par + vUpper}.
#' @param method the method to be used in optimization. See the function \code{stats::optim}.
#'
#' @return a new object of the class PSTR containing the results from the estimation.
#'
#' The object is a list containing the components made in \code{\link{NewPSTR}} and the following new components:
#' \item{iq}{specify which transition variable will be used in estimation.}
#' \item{delta}{the estimate of \eqn{\delta}.}
#' \item{c}{the estimates of \eqn{c}.}
#' \item{vg}{the values of the transition function given the estimates of \eqn{\delta} and \eqn{c} and the transition variables \eqn{q_{it}}.}
#' \item{beta}{the estimates of the coefficient parameters.}
#' \item{vU}{the residuals.}
#' \item{s2}{the variance of the residuals.}
#' \item{cov}{the covariance matrix of the estimates which is cluster-dependency and heteroskedasticity consistent.}
#' \item{est}{a vector of all the estimates}
#' \item{se}{a vector of the standard errors of all the estimates which is cluster-dependency and heteroskedasticity consistent.}
#' \item{mbeta}{a vector of the estimates of the parameters in the second extreme regime.}
#' \item{mse}{a vector of the standard errors of the estimates of the parameters in the second extreme regime.}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso \code{\link{NewPSTR}}, \code{\link{LinTest}} and \code{\link{WCB_LinTest}}
#' @keywords estimation
#'
#' @examples
#' \donttest{
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala'), iT=14) # create a new PSTR object
#'
#' # "L-BFGS-B" is used by default
#' pstr = EstPSTR(use=pstr, im=1, iq=1, par=c(1.6,.5), vLower=4, vUpper=4)
#' # You can also choose the method yourself.
#' pstr = EstPSTR(use=pstr, im=1, iq=1, par=c(1.6,.5), method='CG')
#'
#' print(pstr, "estimates", digits=6)
#' }
#'
#' @export
EstPSTR <- function(use, im=1, iq, par, vLower=2, vUpper=2, method='L-BFGS-B')
{
  if(class(use)!="PSTR")
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  ret = use
  iT = use$iT; iN = use$iN

  # get the data here
  vY = use$vY; vYb = use$vYb
  mX = use$mX; mXb = use$mXb
  mK = use$mK
  ik = ncol(mK)

  vQ = use$mQ[,iq]
  mQ = t(matrix(vQ,iT*iN,im))

  ftmp <- function(vx) return(vx - mean(vx))

  ResiduleSumSquare <- function(vp){
    # vp[1] = log(gamma) or delta
    vg = fTF(vx=mQ,gamma=exp(vp[1]),vc=vp[2:length(vp)])
    mXX = mK * vg
    aXX = array(c(mXX), dim=c(iT,iN,ik))
    mXXb = cbind(mXb, matrix(c(apply(aXX,c(2,3),ftmp)), iT*iN, ik))
    tmp = chol2inv(chol(t(mXXb)%*%mXXb)) %*% t(mXXb) %*% vYb
    vE = c(vYb-mXXb%*%tmp)
    return(sum(vE*vE)/iT/iN)
  }

  if(method=='L-BFGS-B') opt = optim(par=par,fn=ResiduleSumSquare,method="L-BFGS-B",
    lower=par-vLower,upper=par+vUpper)
  else opt = optim(par=par,fn=ResiduleSumSquare,method=method)

  # return value
  ret$iq=iq
  ret$delta = opt$par[1]; ret$gamma = exp(ret$delta)
  ret$c = opt$par[2:length(opt$par)]

  vg = fTF(vx=mQ,gamma=ret$gamma,vc=ret$c) # g_it
  ret$vg = vg
  mXX = mK * vg # x_it * g_it

  aXX = array(c(mXX), dim=c(iT,iN,ik))
  mXXb = cbind(mXb, matrix(c(apply(aXX,c(2,3),ftmp)), iT*iN, ik)) # mean adjust
  tmp = chol2inv(chol(t(mXXb)%*%mXXb)) %*% t(mXXb) %*% vYb # beta
  ret$beta = c(tmp); names(ret$beta) = c(paste0(use$mX_name,'_0'), paste0(use$mK_name,'_1'))

  mXX = cbind(mX, mXX) # (x_it, x_it*g_it)
  ret$mXX = mXX

  ret$vU = c(apply(matrix(c(vY-mXX%*%tmp),iT,iN),2,ftmp))
  ret$s2 = c(ret$vU %*% ret$vU) / (iT*iN)

  # computing standard errors
  tmp = Der2GFunc(vg=vg,vs=vQ,vp=c(ret$gamma,ret$c))
  de1 = tmp$de1; de2 = tmp$de2
  beta1 = ret$beta[(ncol(mX)+1):length(ret$beta)]

  dedp = -mXXb
  d2edp2 = array(0,dim=c(iT*iN,length(ret$beta)+1+im,length(ret$beta)+1+im))
  tmp = 1
  for(iter in 1:ncol(de1)){
    mKK = mK * de1[,iter]
    aKK = array(c(mKK),dim=c(iT,iN,ik))
    mKK = matrix(c(apply(aKK,c(2,3),ftmp)), iT*iN, ik)
    dedp = cbind(dedp, -mKK %*% beta1)

    d2edp2[,(ncol(mX)+1):length(ret$beta),length(ret$beta)+iter] = -mKK
    d2edp2[,length(ret$beta)+iter,(ncol(mX)+1):length(ret$beta)] = -mKK

    for(jter in iter:ncol(de1)){
      mKK = mK * de2[,tmp]
      aKK = array(c(mKK),dim=c(iT,iN,ik))
      mKK = matrix(c(apply(aKK,c(2,3),ftmp)), iT*iN, ik)

      d2edp2[,length(ret$beta)+iter,length(ret$beta)+jter] = -mKK %*% beta1
      d2edp2[,length(ret$beta)+jter,length(ret$beta)+iter] = -mKK %*% beta1
      tmp = tmp+1
    }
  }

  mh = 2*ret$vU*dedp
  ah = array(c(mh),dim=c(iT,iN,ncol(dedp)))
  hi = matrix(c(apply(ah,c(2,3),sum)), iN, ncol(dedp))
  mB = 0
  for(iter in 1:iN){
    mB = mB + hi[iter,] %*% t(hi[iter,])
  }

  invA = 0
  for(iter in 1:(iT*iN))
    invA = invA + (dedp[iter,]%*%t(dedp[iter,]) + d2edp2[iter,,]*ret$vU[iter])*2
  invA = solve(invA)
  # done

  ret$cov = invA %*% mB %*% t(invA)
  ret$se = sqrt(diag(ret$cov))
  names(ret$se) = c(names(ret$beta),'gamma',paste0('c_',1:im))

  ret$est = c(ret$beta, ret$gamma, ret$c)
  names(ret$est) = names(ret$se)

  mM = NULL; mname = NULL
  mTmp = diag(length(ret$mX_name))
  for(iter in 1:length(ret$mX_name)){
    tmp = ret$mX_name[iter] == ret$mK_name
    if(any(tmp)){
      mM = rbind(mM, c(mTmp[iter,], tmp))
      mname = c(mname, ret$mX_name[iter])
    }
  }

  if(!is.null(mM)){
    mM = cbind(mM, matrix(0,nrow(mM),1+im))
    ret$mbeta = c(mM %*% ret$est)
    names(ret$mbeta) = mname
    ret$mse = sqrt(diag(mM %*% ret$cov %*% t(mM)))
    names(ret$mse) = mname
  }

  return(ret)
}




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
#' The functions need the return value (an object of the class PSTR) from the \code{\link{EstPSTR}}. The model should be estimated before conducting the evaluation tests. They copy the object, reuse its contents, especially the estimates, to produce the evaluation test results, and then return a new object of the class PSTR. The user can choose to save the return value to a new object or simply to overwrite the object returned from \code{EstPSTR}. See the example below.
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
#' pstr = EstPSTR(use=pstr, im=1, iq=1, par=c(1.6,.5), method='CG')
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
#' pstr = WCB_TVTest(use=pstr, iB=4, parallel=T, cpus=2)
#'
#' # pstr$mQ[,1] is the transition variable stored in the object
#' # You can also try other transition variables.
#' pstr = WCB_HETest(use=pstr, vq=pstr$mQ[,1], iB=4, parallel=T, cpus=2)
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
  if(class(use)!="PSTR")
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  if(is.null(use$beta))
    stop(simpleError("Estimate the model first!"))
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
  if(class(use)!="PSTR")
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  if(is.null(use$beta))
    stop(simpleError("Estimate the model first!"))
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
    EST = EstPSTR(use=ruse,im=1,iq=ruse$iq,par=c(use$delta,use$c),vLower=1,vUpper=1)
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
    EST = EstPSTR(use=ruse,im=1,iq=ruse$iq,par=c(use$delta,use$c),vLower=1,vUpper=1)
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
    #sfExport(list=c('sLMTEST','EstPSTR','fTF','DerGFunc'))
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
  if(class(use)!="PSTR")
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  if(is.null(use$beta))
    stop(simpleError("Estimate the model first!"))
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
    EST = EstPSTR(use=ruse,im=1,iq=ruse$iq,par=c(use$delta,use$c),vLower=1,vUpper=1)
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
    EST = EstPSTR(use=ruse,im=1,iq=ruse$iq,par=c(use$delta,use$c),vLower=1,vUpper=1)
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
    #sfExport(list=c('iT','iN','vU','eY','mXb','mX','sLMTEST'))
    #sfExport(list=c('sLMTEST','EstPSTR','fTF','DerGFunc'))
    qLM1 = sfSapply(1:iB,ftmp_wb)
    qLM2 = sfSapply(1:iB,ftmp_wcb)
    sfStop()

    ret$wcb_ht = rbind(ret$wcb_ht,c(LM, mean(LM<=qLM1), mean(LM<=qLM2)))
  }

  ret$wcb_ht = matrix(ret$wcb_ht, nrow=im)

  return(ret)
}
