#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## Aug 2023
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
#' The transition function \eqn{g_{it}} takes the logistic form
#' \deqn{g(q_{it} ; \gamma, c) = \left( 1 + \exp \left( - \gamma \prod_{j=1}^{m} (q_{it} - c_j) \right) \right)^{-1}}
#' with \eqn{\gamma > 0} and \eqn{c_1 < c_2 < ... < c_m}. \eqn{\gamma} can be reparametrized as \eqn{\gamma = \exp{\delta}} where \eqn{\delta} is a real number.
#'
#'
#' @section Author and Maintainer:
#' Yukai Yang
#'
#' Department of Statistics, Uppsala University
#'
#' \email{yukai.yang@@statistik.uu.se}
#'
#' @section References:
#' González, A., Teräsvirta, T., van Dijk, D. and Yang, Y. (2005) "\href{http://swopec.hhs.se/hastef/papers/hastef0604.pdf}{Panel Smooth Transition Regression Models}", SSE/EFI Working Paper Series in Economics and Finance 604, Stockholm School of Economics, revised 11 Oct 2017.
#'
#' @section Function for Initialization:
#' \code{\link{NewPSTR}} initialize the modelling by creating an object of the class PSTR.
#'
#' @section Functions for Model Specification:
#' \code{PSTR$LinTest} implements the linearity tests.
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
#' \code{\link{plot_transition}} plots the transition function of an estimated PSTR model.
#'
#' \code{\link{plot_response}} plots curve or surfaces of the expected reponse agaist the corresponding variable.
#'
#' \code{\link{plot_target}} plots the surface of the target function for the nonlinear least square estimation.
#'
#' @section  Data:
#' \code{\link{Hansen99}} a balanced panel of 565 US firms observed for the years 1973–1987.
#'
#' \code{\link{sunspot}} transformed Wolf annual sunspot numbers for the years 1710-1979.
#'
#' @docType package
#' @name PSTR
NULL


#' @importFrom stats optim pchisq pf quantile
NULL

#' @importFrom R6 R6Class
NULL

#' @importFrom knitr kable
NULL

#' @importFrom cli cli_alert_success cli_h1 cli_h2 cli_h3 cli_alert_info
NULL

#' @importFrom magrittr %>%
NULL

#' @import tibble
NULL

#' @importFrom ggplot2 ggplot aes geom_point geom_rug geom_rect geom_line geom_hline facet_grid vars labs scale_x_log10 xlim ylim
NULL

#' @importFrom plotly add_surface add_trace plot_ly layout
NULL

#' @importFrom snowfall sfInit sfExport sfSapply sfStop
NULL


#' Create an object of the class PSTR.
#'
#' Create an object of the R6 class PSTR for later usage. This function should be run prior to the other functions in the package. It will return an object which you will use as an input for the other functions. It builds up the basic settings for the Panel Smooth Transition Regression (PSTR) Modelling.
#'
#' Potential transition variables in \code{tvars} will be tested one by one in, for example, \code{LinTest} function.
#'
#' There is no need to specify the number of individuals,  as it will be obtained automatically inside the function given the number of rows and the sample size \code{iT}.
#'
#' \code{NA}s in \code{data} are removed automatically inside the function.
#'
#' @param data a tibble of data. The number of rows of \code{data} must be the sample size \code{iT} times individuals number N.
#' @param dep column number or name of the dependent variable. Note that this must be specified.
#' @param indep a vector of column numbers of names of the independent variables. Note that this must be specified.
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
  # checking
  if(!is_tibble(data)) stop(simpleError("data should be a tibble!"))
  if(length(dep)>1) stop(simpleError("Only one dependent variable!"))
  if(length(indep)<1) stop(simpleError("There is no independent variable!"))
  if(length(tvars)<1) stop(simpleError("Please specify the candidates of the transition variables 'tvars'."))
  
  #data, dep, indep, indep_k, tvars, im, iT
  ret = PSTR$new(data=data, dep=dep, indep=indep, indep_k=indep_k, tvars=tvars, im=im, iT=iT)
  cli::cli_alert_success("The PSTR model is ready.")
  return(ret)
}


# make the PSTR R6 class
PSTR <- R6::R6Class(
  "PSTR",
  public = list(
    initialize = function(data, dep, indep, indep_k, tvars, im, iT) {
      private$iT = iT
      
      iNN = nrow(data)/iT; coln = c(t(matrix(1:iNN, iNN, iT)))
      
      vY = data[,dep]; private$vY_name =  names(data[,dep])
      mX = data[,indep]; private$mX_name = names(data[,indep])
      
      if(is.null(indep_k)){
        mK = mX; private$mK_name = private$mX_name
      }else{
        mK = data[,indep_k]; private$mK_name = names(data[,indep_k])
      }
      
      mQ = data[,tvars]; private$mQ_name = names(data[,tvars])
      
      # remove the NAs
      tmp = is.na(vY) | apply(is.na(mX), 1, any) | apply(is.na(mK), 1, any) | apply(is.na(mQ), 1, any)
      tmp = is.na(match(coln,coln[tmp]))
      vY = vY[tmp,]; mX = mX[tmp,]; mK = mK[tmp,]; mQ = mQ[tmp,]
      private$vY = c(as.matrix(vY))
      private$mX = as.matrix(mX)
      private$mK = as.matrix(mK)
      private$mQ = as.matrix(mQ)
      private$im = im
      
      iN = sum(tmp)/iT; private$iN = iN
      
      coln = c(t(matrix(1:iN, iN, iT)))
      vYb = NULL; mXb = NULL
      for(nter in 1:iN){
        tmp = coln==nter
        vYb = c(vYb, private$vY[tmp] - mean(private$vY[tmp]))
        mXb = rbind(mXb, t(t(mX[tmp,])-apply(t(mX[tmp,]),1,mean)))
      }
      
      private$vYb = vYb; private$mXb = mXb
    },
    getTest = function(){private$test},
    getSqTest = function(){private$sqtest}
  ),
  private = list(
    iT=NULL, vY_name=NULL, mX_name=NULL, mK_name=NULL, mQ_name=NULL,
    vY=NULL, mX=NULL, mK=NULL, mQ=NULL, im=NULL, iN=NULL,
    vYb=NULL, mXb=NULL,
    # assigned by LinTest
    test=list(), sqtest=list()
  )
)
