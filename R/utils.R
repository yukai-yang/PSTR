#################################################################################
## package name: PSTR.R
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## Aug 2017
#################################################################################

#################################################################################
## utility functions
#################################################################################



# simple cat
cat0 <- function(...)
{
  words = list(...)
  for(tmp in words) cat(tmp)
  cat("\n")
}


#' Show the version number of some information.
#'
#' This function shows the version number and some information of the package.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @keywords utils
#'
#' @export
version <- function(){
  cat0("#########################################################################")
  cat0("## package name: PSTR")
  cat0("## author: Yukai Yang")
  cat0("## Department of Statistics")
  cat0("## Uppsala University")
  cat0("## yukai.yang@statistik.uu.se")
  cat0("## Version 1.0.0, Sep. 2017")
  cat0("#########################################################################")
}



#' Print the object of the class PSTR.
#'
#' This function prints the object of the class PSTR.
#'
#' @param x an object of the class PSTR returned from some functions in the package. See below "See Also" for a list of these functions.
#' @param mode a vector of character strings specifying which results to print. It takes the values c('summary', 'tests', 'estimates', 'evaluation'). By default 'summary'.
#' @param digits integer indicating the number of decimal places (for the \code{round} function inside) to be used. Negative values are allowed (see \code{round}).
#' @param ... further arguments passed to or from other methods. Ignored here.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso Functions which return an object of the class PSTR:
#'
#' \code{\link{NewPSTR}}, \code{\link{LinTest}}, \code{\link{WCB_LinTest}}, \code{\link{EstPSTR}}, \code{\link{EvalTest}}, \code{\link{WCB_TVTest}} and \code{\link{WCB_HETest}}
#' @keywords utils
#'
#' @examples
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala','debta','cfa','sales'), iT=14)
#' print(pstr)
#' print(pstr, mode='summary',digits=2)
#' @export
print.PSTR <- function(x, mode="summary", digits=4, ...)
{
  cat0("#########################################################################")
  cat0("## package name: PSTR")
  cat0("## Version 1.0.0, Sep. 2017")

  tmp = NULL
  for(iter in 1:length(mode)){
    tmp = c(tmp, grep(mode[iter], c("summary","tests","estimates","evaluation")))
  }
  tmp = unique(tmp)

  if(1 %in% tmp) print_summary(x)

  if(2 %in% tmp) print_tests(x,digits)

  if(3 %in% tmp) print_estimates(x,digits)

  if(4 %in% tmp) print_evaluation(x,digits)

  if(length(tmp)==0){
    cat0("***********************************************************************")
    cat0("The argument 'mode' only accepts the values:")
    cat0("  'summary'(default), 'tests', 'estimates' or 'evaluation'.")
    cat0("Incomplete words, such as 'su' or 'mm' for 'summary', are allowed.")
  }

  cat0("***********************************************************************")
  cat0("#########################################################################")
}


print_summary <- function(obj)
{
  cat0("#########################################################################")
  cat0("***********************************************************************")
  cat0("Summary of the model:")

  cat0("-----------------------------------------------------------------------")
  cat0("  time horizon sample size = ",obj$iT,",  number of individuals = ",obj$iN)

  cat0("-----------------------------------------------------------------------")
  cat0("Dependent variable:  ",obj$vY_name)
  cat0("-----------------------------------------------------------------------")
  cat0("Explanatory variables in the linear part:")
  cat0("  ",obj$mX_name)
  cat0("-----------------------------------------------------------------------")
  cat0("Explanatory variables in the non-linear part:")
  cat0("  ",obj$mK_name)
  cat0("-----------------------------------------------------------------------")
  cat0("Potential transition variable(s) to be tested:")
  cat0("  ",obj$mQ_name)
}


print_tests <- function(obj,digits)
{
  cat0("#########################################################################")
  cat0("***********************************************************************")
  cat0("Results of the linearity (homogeneity) tests:")

  if(!is.null(obj$test)){
    for(iter in 1:length(obj$test)){
      im = length(obj$test[[iter]])

      cat0("-----------------------------------------------------------------------")
      cat0("LM tests based on transition variable '",obj$mQ_name[iter],"'")

      tmp = NULL
      for(jter in 1:im){
        ttmp = obj$test[[iter]][[jter]]
        tmp = rbind(tmp, c(jter, ttmp$LM1_X, ttmp$PV1_X, ttmp$LM1_F, ttmp$PV1_F,
                           ttmp$LM2_X, ttmp$PV2_X, ttmp$LM2_F, ttmp$PV2_F) )
      }
      tmp = matrix(tmp, nrow=im)
      rownames(tmp) = rep(" ",im)
      colnames(tmp) = c("m", "LM_X", "PV", "LM_F", "PV", "HAC_X", "PV", "HAC_F", "PV")

      if(!is.null(obj$wcb_test)){
        ttmp = obj$wcb_test[[iter]]
        ttmp = matrix(ttmp, nrow=im)[,2:3]
        colnames(ttmp) = c("WB_PV", "WCB_PV")
        tmp = cbind(tmp, ttmp)
      }

      print(signif(tmp,digits))

    }
  }

  cat0("***********************************************************************")
  cat0("Sequence of homogeneity tests for selecting number of switches 'm':")

  if(!is.null(obj$sqtest)){
    for(iter in 1:length(obj$sqtest)){
      im = length(obj$sqtest[[iter]])

      cat0("-----------------------------------------------------------------------")
      cat0("LM tests based on transition variable '",obj$mQ_name[iter],"'")

      tmp = NULL
      for(jter in 1:im){
        ttmp = obj$sqtest[[iter]][[jter]]
        tmp = rbind(tmp, c(jter, ttmp$LM1_X, ttmp$PV1_X, ttmp$LM1_F, ttmp$PV1_F,
                           ttmp$LM2_X, ttmp$PV2_X, ttmp$LM2_F, ttmp$PV2_F) )
      }
      tmp = matrix(tmp, nrow=im)
      rownames(tmp) = rep(" ",im)
      colnames(tmp) = c("m", "LM_X", "PV", "LM_F", "PV", "HAC_X", "PV", "HAC_F", "PV")

      if(!is.null(obj$wcb_sqtest)){
        ttmp = obj$wcb_sqtest[[iter]]
        ttmp = matrix(ttmp, nrow=im)[,2:3]
        colnames(ttmp) = c("WB_PV", "WCB_PV")
        tmp = cbind(tmp, ttmp)
      }

      print(signif(tmp,digits))

    }
  }
}


print_estimates <- function(obj,digits)
{
  cat0("#########################################################################")
  cat0("***********************************************************************")

  cat0("Results of the PSTR estimation:")

  if(!is.null(obj$iq)){
    cat0("-----------------------------------------------------------------------")
    cat0("Transition variable '",obj$mQ_name[obj$iq],"' is used in the estimation.")
    cat0("-----------------------------------------------------------------------")
    cat0("Parameter estimates in the linear part (first extreme regime) are")
    tmp = rbind(obj$beta[1:length(obj$mX_name)],obj$se[1:length(obj$mX_name)])
    rownames(tmp) = c('Est','s.e.')
    print(signif(tmp,digits))
    cat0("-----------------------------------------------------------------------")
    cat0("Parameter estimates in the non-linear part are")
    tmp = rbind(obj$beta[(length(obj$mX_name)+1):length(obj$beta)],obj$se[(length(obj$mX_name)+1):length(obj$beta)])
    rownames(tmp) = c('Est','s.e.')
    print(signif(tmp,digits))
    cat0("-----------------------------------------------------------------------")
    cat0("Parameter estimates in the second extreme regime are")
    tmp = rbind(obj$mbeta,obj$mse)
    rownames(tmp) = c('Est','s.e.')
    colnames(tmp) = paste0(colnames(tmp),'_{0+1}')
    print(signif(tmp,digits))
    cat0("-----------------------------------------------------------------------")
    cat0("Non-linear parameter estimates are")
    tmp = rbind(obj$est[(length(obj$beta)+1):length(obj$est)],obj$se[(length(obj$beta)+1):length(obj$se)])
    rownames(tmp) = c('Est','s.e.')
    print(signif(tmp,digits))
    cat0("-----------------------------------------------------------------------")
    cat0("Estimated standard deviation of the residuals is ",signif(sqrt(obj$s2),digits))
  }
}


print_evaluation <- function(obj,digits)
{
  cat0("#########################################################################")
  cat0("***********************************************************************")

  cat0("Results of the evaluation tests:")

  if(!is.null(obj$tv)){
    cat0("-----------------------------------------------------------------------")
    cat0("Parameter constancy test")
    im = length(obj$tv)

    tmp = NULL
    for(jter in 1:im){
      ttmp = obj$tv[[jter]]
      tmp = rbind(tmp, c(jter, ttmp$LM1_X, ttmp$PV1_X, ttmp$LM1_F, ttmp$PV1_F,
                         ttmp$LM2_X, ttmp$PV2_X, ttmp$LM2_F, ttmp$PV2_F) )
    }
    tmp = matrix(tmp, nrow=im)
    rownames(tmp) = rep(" ",im)
    colnames(tmp) = c("m", "LM_X", "PV", "LM_F", "PV", "HAC_X", "PV", "HAC_F", "PV")

    print(signif(tmp,digits))
  }

  if(!is.null(obj$wcb_tv)){
    cat0("-----------------------------------------------------------------------")
    cat0("WB and WCB parameter constancy test")
    tmp = obj$wcb_tv[,2:3,drop=F]
    im = nrow(tmp)
    tmp = cbind(1:im, tmp)
    rownames(tmp) = rep(" ",im)
    colnames(tmp) = c("m","WB_PV", "WCB_PV")

    print(signif(tmp,digits))
  }

  if(!is.null(obj$ht)){
    cat0("-----------------------------------------------------------------------")
    cat0("No remaining nonliearity (heterogeneity) test")
    im = length(obj$ht)

    tmp = NULL
    for(jter in 1:im){
      ttmp = obj$ht[[jter]]
      tmp = rbind(tmp, c(jter, ttmp$LM1_X, ttmp$PV1_X, ttmp$LM1_F, ttmp$PV1_F,
                         ttmp$LM2_X, ttmp$PV2_X, ttmp$LM2_F, ttmp$PV2_F) )
    }
    tmp = matrix(tmp, nrow=im)
    rownames(tmp) = rep(" ",im)
    colnames(tmp) = c("m", "LM_X", "PV", "LM_F", "PV", "HAC_X", "PV", "HAC_F", "PV")

    print(signif(tmp,digits))
  }

  if(!is.null(obj$wcb_ht)){
    cat0("-----------------------------------------------------------------------")
    cat0("WB and WCB no remaining nonliearity (heterogeneity) test")
    tmp = obj$wcb_ht[,2:3,drop=F]
    im = nrow(tmp)
    tmp = cbind(1:im, tmp)
    rownames(tmp) = rep(" ",im)
    colnames(tmp) = c("m","WB_PV", "WCB_PV")

    print(signif(tmp,digits))
  }

}
