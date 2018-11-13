#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## Sep 2017
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
#' @param par initial values for the parameters \eqn{\gamma} or \eqn{\delta}, and \eqn{c} to be optimized over. It is a vector of length \code{im}+1, where \code{im} is the number of switches. When missing, the function will choose the initial values automatically, and \code{useDelta=TRUE}.
#' @param useDelta whether delta is used in par in the estimation. Note that if \code{par} is missing, this argument will be ignored.
#' @param vLower a vector or number of the lower offsets determining the lower bounds of the parameters. The lower bounds of the parameters are \code{par - vLower}.
#' @param vUpper a vector or number of the upper offsets determining the upper bounds of the parameters. The upper bounds of the parameters are \code{par + vUpper}.
#' @param method the method to be used in optimization. See the function \code{stats::optim}.
#'
#' @return a new object of the class PSTR containing the results from the estimation.
#'
#' The object is a list containing the components made in \code{\link{NewPSTR}} and the following new components:
#' \item{iq}{specify which transition variable will be used in estimation. The default value \code{NULL} implies a linear panel regression model.}
#' \item{delta}{the estimate of \eqn{\delta}.}
#' \item{c}{the estimates of \eqn{c}.}
#' \item{vg}{the values of the transition function given the estimates of \eqn{\delta} and \eqn{c} and the transition variables \eqn{q_{it}}.}
#' \item{beta}{the estimates of the coefficient parameters.}
#' \item{vU}{the residuals.}
#' \item{vM}{a vector of the estimated time-invariant individual effect.}
#' \item{s2}{the variance of the residuals.}
#' \item{cov}{the covariance matrix of the estimates which is cluster-dependency and heteroskedasticity consistent.}
#' \item{est}{a vector of all the estimates}
#' \item{se}{a vector of the standard errors of all the estimates which is cluster-dependency and heteroskedasticity consistent.}
#' \item{mbeta}{a vector of the estimates of the parameters in the second extreme regime.}
#' \item{mse}{a vector of the standard errors of the estimates of the parameters in the second extreme regime.}
#' \item{convergence}{an integer code showing the convergence, see \code{optim}.}
#' \item{par}{a vector of the initial values used in the optimization. Note that the first element is always delta, no matter whether gamma is used as input.}
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
#' # estimate a linear panel regression model
#' pstr = EstPSTR(use=pstr)
#' print(pstr, "estimates", digits=6)
#'
#' # "L-BFGS-B" is used by default
#' pstr = EstPSTR(use=pstr, im=1, iq=1, useDelta=TRUE, par=c(1.6,.5), vLower=4, vUpper=4)
#' # You can also choose the method yourself.
#' pstr = EstPSTR(use=pstr, im=1, iq=1, useDelta=TRUE, par=c(1.6,.5), method='CG')
#'
#' print(pstr, "estimates", digits=6)
#' 
#' # The estimation of a linear panel regression model with fix effects is also implemented.
#' pstr0 = EstPSTR(use=pstr)
#' 
#' print(pstr0,"estimates")
#' }
#'
#' @export
EstPSTR <- function(use, im=1, iq=NULL, par=NULL, useDelta=FALSE, vLower=2, vUpper=2, method='L-BFGS-B')
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
  
  ftmp <- function(vx) return(vx - mean(vx))
  
  ret$imm = im # used in estimation
  
  if(is.numeric(iq)) ret$iq=iq
  else ret$iq=which(use$mQ_name==iq)
  
  if(length(ret$iq)>1) stop(simpleError("Sorry! We only support the one transition variable case."))
  if(length(iq)<1) stop(simpleError("Fail to find the transition variable."))
  
  if(!is.null(iq)){ 
    if(im < 1) stop(simpleError("The number of switches is invalid."))
    
    vQ = use$mQ[,iq]
    mQ = t(matrix(vQ,iT*iN,im))  
    
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
    
    if(is.null(par)){
      useDelta = T
      tmp = unname(quantile(vQ, (1:im) / (im+1)))
      par = c(log(8/min(diff(c(0,tmp)))), tmp)
    }
    
    if(!useDelta) par[1] = log(par[1])
    
    ret$par = par
    
    if(method=='L-BFGS-B') opt = optim(par=par,fn=ResiduleSumSquare,method="L-BFGS-B",lower=par-vLower,upper=par+vUpper)
    else opt = optim(par=par,fn=ResiduleSumSquare,method=method)
    
    # return value 
    ret$delta = opt$par[1]; ret$gamma = exp(ret$delta)
    ret$c = opt$par[2:length(opt$par)]
    ret$convergence = opt$convergence
    
    vg = fTF(vx=mQ,gamma=ret$gamma,vc=ret$c) # g_it
    ret$vg = vg
    mXX = mK * vg # x_it * g_it 
    aXX = array(c(mXX), dim=c(iT,iN,ik))
    mXXb = cbind(mXb, matrix(c(apply(aXX,c(2,3),ftmp)), iT*iN, ik)) # mean adjust
    tmp = chol2inv(chol(t(mXXb)%*%mXXb)) %*% t(mXXb) %*% vYb # beta
    ret$beta = c(tmp); names(ret$beta) = c(paste0(use$mX_name,'_0'), paste0(use$mK_name,'_1'))
    
    mXX = cbind(mX, mXX) # (x_it, x_it*g_it)
    ret$mXX = mXX
    
    mtmp = matrix(c(vY-mXX%*%tmp),iT,iN)
    ret$vM = c(apply(mtmp,2,mean))
    ret$vU = c(apply(mtmp,2,ftmp))
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
    ttmp = try(solve(invA), silent=T)
    if(class(ttmp)=='try-error'){
      ttmp = svd(invA); invA = ttmp$u %*% diag(1/ttmp$d) %*% t(ttmp$u)
    }else invA = ttmp
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
 
  }else{
    
    tmp = chol2inv(chol(t(mXb)%*%mXb)) %*% t(mXb) %*% vYb
    #vE = c(vYb-mXb%*%tmp)
    
    # return value 
    ret$beta = c(tmp); names(ret$beta) = use$mX_name
    
    mtmp = matrix(c(vY-mX%*%tmp),iT,iN)
    ret$vM = c(apply(mtmp,2,mean))
    ret$vU = c(apply(mtmp,2,ftmp))
    ret$s2 = c(ret$vU %*% ret$vU) / (iT*iN)
    
    # computing standard errors
    dedp = -mXb
    
    mh = 2*ret$vU*dedp
    ah = array(c(mh),dim=c(iT,iN,ncol(dedp)))
    hi = matrix(c(apply(ah,c(2,3),sum)), iN, ncol(dedp))
    mB = 0
    for(iter in 1:iN){
      mB = mB + hi[iter,] %*% t(hi[iter,])
    }
    
    invA = 0
    for(iter in 1:(iT*iN))
      invA = invA + dedp[iter,]%*%t(dedp[iter,])*2
    ttmp = try(solve(invA), silent=T)
    if(class(ttmp)=='try-error'){
      ttmp = svd(invA); invA = ttmp$u %*% diag(1/ttmp$d) %*% t(ttmp$u)
    }else invA = ttmp
    # done
    
    ret$cov = invA %*% mB %*% t(invA)
    ret$se = sqrt(diag(ret$cov))
    names(ret$se) = names(ret$beta)
    
    ret$est = ret$beta
    names(ret$est) = names(ret$se)
    
  }
  
  return(ret)
}