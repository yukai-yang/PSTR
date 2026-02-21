#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
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
#' \code{\link{plot_transition}} plots the transition function of an estimated PSTR model.
#' 
#' \code{\link{plot_coefficients}} plots coefficients, standard errors, and p-values against the transition variable.
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
#' @keywords internal
"_PACKAGE"


#' @importFrom stats optim pchisq pf quantile
NULL

#' @importFrom R6 R6Class
NULL

#' @importFrom knitr kable
NULL

#' @importFrom cli cli_alert_success cli_h1 cli_h2 cli_h3 cli_alert_info cli_alert_warning cli_alert_danger
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


#' Create a PSTR model object
#'
#' Create an R6 object of class \code{"PSTR"} to be used as the main container for
#' Panel Smooth Transition Regression (PSTR) modelling in this package.
#' You typically call \code{NewPSTR()} once, and then pass the returned object to
#' specification, estimation and evaluation functions.
#'
#' The candidate transition variables in \code{tvars} will be stored in the object
#' and can be tested one by one by functions such as \code{\link{LinTest}}.
#'
#' Missing values in the dependent variable, linear regressors, non-linear regressors,
#' or transition variables are removed internally (row-wise).
#' The number of individuals \eqn{N} is inferred from \code{nrow(data)} and \code{iT}
#' after removing missing values.
#'
#' @param data A tibble containing the panel in long format. The number of rows must be
#'   \code{iT * N} for some integer \code{N}. Rows are assumed to be ordered by time within
#'   individual, consistently with the package conventions.
#' @param dep A single column index or a single column name specifying the dependent variable.
#' @param indep A vector of column indices or column names specifying the regressors in the
#'   linear part.
#' @param indep_k Optional. A vector of column indices or column names specifying the regressors
#'   in the non-linear part. If \code{NULL}, the non-linear part is set equal to the linear part.
#' @param tvars A vector of column indices or column names specifying the candidate transition
#'   variables.
#' @param im Integer. The maximal number of switches used in linearity-related tests.
#'   Default is \code{1}.
#' @param iT Integer. The time dimension (number of time observations per individual).
#'
#' @return An R6 object of class \code{"PSTR"}.
#'
#' @seealso \code{\link{LinTest}}, \code{\link{WCB_LinTest}}, \code{\link{EstPSTR}},
#'   \code{\link{EvalTest}}, \code{\link{WCB_TVTest}}, \code{\link{WCB_HETest}}.
#'
#' @keywords initialization
#'
#' @examples
#' \donttest{
#' pstr <- NewPSTR(
#'   Hansen99,
#'   dep = "inva",
#'   indep = 4:20,
#'   indep_k = c("vala", "debta", "cfa", "sales"),
#'   tvars = c("vala", "debta"),
#'   iT = 14
#' )
#'
#' # print summary (your R6 print method)
#' pstr
#' print(pstr, mode = "summary")
#'
#' # after running tests/estimation, you can print other sections
#' # print(pstr, mode = "tests")
#' # print(pstr, mode = "estimates")
#' # print(pstr, mode = "evaluation")
#' }
#'
#' @export
NewPSTR <- function(data, dep, indep, indep_k=NULL, tvars, im=1, iT)
{
  # checking
  if(!is_tibble(data)){cli::cli_alert_danger("data should be a tibble!"); return(0)}
  if(length(dep)>1){cli::cli_alert_danger("Only one dependent variable!"); return(0)}
  if(length(indep)<1){cli::cli_alert_danger("There is no independent variable!"); return(0)}
  if(length(tvars)<1){cli::cli_alert_danger("Please specify the candidates of the transition variables 'tvars'."); return(0)}
  
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
    
    # assigned by LinTest
    test=NULL, sqtest=NULL,
    # assigned by WCB_LinTest
    wcb_test=NULL, wcb_sqtest=NULL,
    
    # estimation results
    delta=NULL, gamma=NULL, c=NULL,
    vg=NULL, beta=NULL,
    vM=NULL, vU=NULL, s2=NULL,
    cov=NULL, se=NULL, est=NULL,
    mbeta=NULL, mse=NULL,
    
    # assigned by EvalTest
    tv=NULL, ht=NULL,
    
    # assigned by WCB_TVTest
    wcb_tv=NULL,
    
    # assigned by WCB_HETest
    wcb_ht = NULL
  ),
  private = list(
    iT=NULL, vY_name=NULL, mX_name=NULL, mK_name=NULL, mQ_name=NULL,
    vY=NULL, mX=NULL, mK=NULL, mQ=NULL, im=NULL, iN=NULL,
    vYb=NULL, mXb=NULL,
    
    # estimation results
    imm=NULL, iq=NULL, par=NULL, convergence=NULL,
    mXX=NULL
  )
)
