#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## May 2018
#################################################################################

#################################################################################
## utility functions
#################################################################################

vnum = "1.2.5"
packname = "(Orange Panel)"

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
  cat0("PSTR version ", vnum, " ",packname)
}


#' Print the object of the class PSTR.
#'
#' This function prints the object of the class PSTR.
#'
#' @param x an object of the class PSTR returned from some functions in the package. See below "See Also" for a list of these functions.
#' @param mode a vector of character strings specifying which results to print. It takes the values c('summary', 'tests', 'estimates', 'evaluation'). By default 'su' and 'e' which means all.
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
print.PSTR <- function(x, mode=c("su","e"), digits=4, ...)
{
  cat0(paste0(rep("#",getOption("width")),collapse=''))
  cat0("## PSTR ", vnum, " ",packname)

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
    cat0(paste0(rep("*",getOption("width")),collapse=''))
    cat0("The argument 'mode' only accepts the values:")
    cat0("  'summary', 'tests', 'estimates' or 'evaluation'.")
    cat0("Incomplete words, such as 'su' or 'mm' for 'summary', are allowed.")
  }

  cat0(paste0(rep("*",getOption("width")),collapse=''))
  cat0(paste0(rep("#",getOption("width")),collapse=''))
}


print_summary <- function(obj)
{
  cat0(paste0(rep("#",getOption("width")),collapse=''))
  cat0(paste0(rep("*",getOption("width")),collapse=''))
  cat0("Summary of the model:")

  cat0(paste0(rep("-",getOption("width")),collapse=''))
  cat0("  time horizon sample size = ",obj$iT,",  number of individuals = ",obj$iN)

  cat0(paste0(rep("-",getOption("width")),collapse=''))
  cat0("Dependent variable:  ",obj$vY_name)
  cat0(paste0(rep("-",getOption("width")),collapse=''))
  cat0("Explanatory variables in the linear part:")
  cat0("  ",obj$mX_name)
  cat0(paste0(rep("-",getOption("width")),collapse=''))
  cat0("Explanatory variables in the non-linear part:")
  cat0("  ",obj$mK_name)
  cat0(paste0(rep("-",getOption("width")),collapse=''))
  cat0("Potential transition variable(s) to be tested:")
  cat0("  ",obj$mQ_name)
}


print_tests <- function(obj,digits)
{
  cat0(paste0(rep("#",getOption("width")),collapse=''))
  cat0(paste0(rep("*",getOption("width")),collapse=''))
  cat0("Results of the linearity (homogeneity) tests:")

  im = obj$im

  if(!is.null(obj$test)){
    for(iter in 1:length(obj$test)){

      cat0(paste0(rep("-",getOption("width")),collapse=''))
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
        ttmp = obj$wcb_test[[iter]][,2:3,drop=F]
        colnames(ttmp) = c("WB_PV", "WCB_PV")
        tmp = cbind(tmp, ttmp)
      }

      print(signif(tmp,digits))
    }
  }else{
    if(!is.null(obj$wcb_test)){
      for(iter in 1:length(obj$wcb_test)){

        cat0(paste0(rep("-",getOption("width")),collapse=''))
        cat0("LM tests based on transition variable '",obj$mQ_name[iter],"'")

        ttmp = cbind(1:im, obj$wcb_test[[iter]][,2:3,drop=F])
        rownames(ttmp) = rep(" ",im)
        colnames(ttmp) = c("m","WB_PV", "WCB_PV")
      }

      print(signif(ttmp,digits))
    }
  }

  cat0(paste0(rep("*",getOption("width")),collapse=''))
  cat0("Sequence of homogeneity tests for selecting number of switches 'm':")

  if(!is.null(obj$sqtest)){
    for(iter in 1:length(obj$sqtest)){

      cat0(paste0(rep("-",getOption("width")),collapse=''))
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
        ttmp = obj$wcb_sqtest[[iter]][,2:3,drop=F]
        colnames(ttmp) = c("WB_PV", "WCB_PV")
        tmp = cbind(tmp, ttmp)
      }

      print(signif(tmp,digits))
    }
  }
}


print_estimates <- function(obj,digits)
{
  if(!is.null(obj$iq)){
    cat0(paste0(rep("#",getOption("width")),collapse=''))
    cat0(paste0(rep("*",getOption("width")),collapse=''))
    cat0("Results of the PSTR estimation:")
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    cat0("Transition variable '",obj$mQ_name[obj$iq],"' is used in the estimation.")
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    cat0("Parameter estimates in the linear part (first extreme regime) are")
    tmp = rbind(obj$beta[1:length(obj$mX_name)],obj$se[1:length(obj$mX_name)])
    rownames(tmp) = c('Est','s.e.')
    print(signif(tmp,digits))
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    cat0("Parameter estimates in the non-linear part are")
    tmp = rbind(obj$beta[(length(obj$mX_name)+1):length(obj$beta)],obj$se[(length(obj$mX_name)+1):length(obj$beta)])
    rownames(tmp) = c('Est','s.e.')
    print(signif(tmp,digits))
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    cat0("Parameter estimates in the second extreme regime are")
    tmp = rbind(obj$mbeta,obj$mse)
    rownames(tmp) = c('Est','s.e.')
    colnames(tmp) = paste0(colnames(tmp),'_{0+1}')
    print(signif(tmp,digits))
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    cat0("Non-linear parameter estimates are")
    tmp = rbind(obj$est[(length(obj$beta)+1):length(obj$est)],obj$se[(length(obj$beta)+1):length(obj$se)])
    rownames(tmp) = c('Est','s.e.')
    print(signif(tmp,digits))
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    cat0("Estimated standard deviation of the residuals is ",signif(sqrt(obj$s2),digits))
  }else{
    if(!is.null(obj$est)){
      cat0(paste0(rep("#",getOption("width")),collapse=''))
      cat0(paste0(rep("*",getOption("width")),collapse=''))
      cat0("A linear panel regression with fixed effects is estimated.")
      cat0(paste0(rep("-",getOption("width")),collapse=''))
      cat0("Parameter estimates are")
      tmp = rbind(obj$est,obj$se)
      rownames(tmp) = c('Est','s.e.')
      print(signif(tmp,digits))
      cat0(paste0(rep("-",getOption("width")),collapse=''))
      cat0("Estimated standard deviation of the residuals is ",signif(sqrt(obj$s2),digits))
    }
  }
}


print_evaluation <- function(obj,digits)
{
  if(!is.null(obj$tv) || !is.null(obj$wcb_tv) || !is.null(obj$ht) || !is.null(obj$wcb_ht)){
    cat0(paste0(rep("#",getOption("width")),collapse=''))
    cat0(paste0(rep("*",getOption("width")),collapse=''))

    cat0("Results of the evaluation tests:")
  }

  if(!is.null(obj$tv)){
    cat0(paste0(rep("-",getOption("width")),collapse=''))
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
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    cat0("WB and WCB parameter constancy test")
    tmp = obj$wcb_tv[,2:3,drop=F]
    im = nrow(tmp)
    tmp = cbind(1:im, tmp)
    rownames(tmp) = rep(" ",im)
    colnames(tmp) = c("m","WB_PV", "WCB_PV")

    print(signif(tmp,digits))
  }

  if(!is.null(obj$ht)){
    cat0(paste0(rep("-",getOption("width")),collapse=''))
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
    cat0(paste0(rep("-",getOption("width")),collapse=''))
    cat0("WB and WCB no remaining nonliearity (heterogeneity) test")
    tmp = obj$wcb_ht[,2:3,drop=F]
    im = nrow(tmp)
    tmp = cbind(1:im, tmp)
    rownames(tmp) = rep(" ",im)
    colnames(tmp) = c("m","WB_PV", "WCB_PV")

    print(signif(tmp,digits))
  }

}


#' Plot the transition function of the estimated PSTR model.
#'
#' This function plots the transition function of the estimated PSTR model.
#'
#' The funciton uses some functions in the ggplot2 package and aims to give a quick plot of the transtion function.
#' The user can customize the title, subtitle, caption, x and y labels, for details, read the help file for the \code{labs} function in ggplot2.
#'
#' @param obj an object of the class PSTR returned from some functions in the package. Note that the corresponding PSTR model must be estimated first.
#' @param size the size of the circle.
#' @param color the color of the circle.
#' @param xlim a numeric vector of dimension 2 specifying the limits of x-axis.
#' @param ylim a numeric vector of dimension 2 specifying the limits of y-axis.
#' @param fill the color used to fill the area on the transition curve with observations, \code{NULL} for not fill, see ggplot2.
#' @param alpha a number controlling the transparency of the points and filled area, \code{NULL} for default, see ggplot2.
#'
#' @return A ggplot object. The user can plot it simply by print the object.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso Functions which return an object of the class PSTR and can be input into this function
#'
#' \code{\link{EstPSTR}}
#' @keywords utils
#'
#' @examples
#' \donttest{
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala'), iT=14) # create a new PSTR object
#'
#' # estimate the PSTR model
#' pstr = EstPSTR(use=pstr, im=1, iq=1, useDelta=TRUE, par=c(.63,0), method='CG')
#'
#' # plot the transition function
#'
#' ret = plot_transition(pstr)
#' # plot by running
#' ret
#'
#' ret = plot_transition(pstr, fill='blue', xlim=c(-2,20), color = "dodgerblue4", size = 2, alpha=.3)
#' ret
#' }
#'
#' @export
plot_transition <- function(obj, size=1.5, color="blue", xlim=NULL, ylim=NULL, fill=NULL, alpha=NULL)
{
  if(class(obj)!="PSTR")
    stop(simpleError("The argument 'obj' is not an object of class 'PSTR'."))
  if(is.null(obj$vg)) stop(simpleError("The PSTR model is not estimated yet."))
  
  qq=obj$mQ[,obj$iq]

  ###### new
  if(is.null(xlim)){ xlim = c(min(min(qq), obj$c - log(1/0.002472623 - 1)/obj$gamma),
                             max(max(qq), obj$c - log(1/0.9975274 - 1)/obj$gamma)) 
  }else{
    if(length(xlim)!=2) stop(simpleError("xlim must be a 2-vector."))
    if(!is.numeric(xlim)) stop(simpleError("xlim must be numeric."))
  }
  
  if(is.null(ylim)){ ylim = c(0, 1)
  }else{
    if(length(ylim)!=2) stop(simpleError("ylim must be a 2-vector."))
    if(!is.numeric(ylim)) stop(simpleError("ylim must be numeric."))
  }
  
  if(is.null(alpha)) alpha = .2
  
  vx = seq(xlim[1],xlim[2],length.out=1001)
  mx = t(matrix(vx,length(vx),obj$imm))
  vy = fTF(mx, obj$gamma, obj$c)
  
  ret = ggplot() + ggplot2::xlim(xlim) + ggplot2::ylim(ylim)
  if(!is.null(fill)) ret = ret + geom_rect(aes(xmin=min(qq),ymin=0,xmax=max(qq),ymax=1), alpha=alpha/2, fill=fill) 
  ret = ret + geom_line(aes(x=vx,y=vy),color='red') +
    geom_rug(aes(x=qq, y=obj$vg), sides = "b", color=color) +
    geom_point(aes(x=qq, y= obj$vg), size=size, stroke=T, alpha=alpha, color=color) +
    labs(y="transition function", x=obj$mQ_name[obj$iq])
  ######

  return(ret)
}


#' Curve or surfaces of the expected reponse agaist the corresponding variable.
#'
#' This function plots the curve or the surfaces of the expected reponse agaist the corresponding variable (and the transition variable if surface).
#'
#' The expected response is the expected value of the dependent variable minus the individual effect and all the other variables times their estimated coefficients.
#' That is, if the variable is \eqn{z_{k,it}} in both \eqn{x_{it}} and \eqn{z_{it}},
#' then the function plots the surface of
#' \deqn{y_{it} - \mu_i - \beta_{-k,0}' x_{-k,it} + \beta_{-k,1}' z_{-k,it} g_{it} - u_{it}}
#' or simply
#' \deqn{(\beta_{k,0} + \beta_{k,1}g_{it}) \cdot z_{k,it}}
#' where \eqn{-k} means with the \eqn{k}th element removed,
#' against \eqn{z_{k,it}} and \eqn{q_{it}} if \eqn{z_{k,it} \neq q_{it}}.
#'
#' If \eqn{z_{k,it} = q_{it}}, then the function plots the curve of the expected response defined above against \eqn{z_{k,it}}.
#'
#' More than one variable can be put in \code{vars}.
#' If \code{vars} contains the transition variable and the transition variable belongs to the nonlinear part,
#' the function will plot a curve of the effect-adjusted expected response and the transition variable,
#' otherwise, the function will plot a 3-D surface of the effect-adjusted expected response against a chosen variable in the nonlinear part and the transition variable.
#'
#' \code{length.out} takes a vector or a scalar.
#' The vector must be two dimensional specifying numbers of points in the grid built for the surface.
#' The first element of the vector corresponds to the variables, and the second to the transition variable.
#' If it is a scalar, then grid has the same number of points for the variables and the transition varible.
#'
#' The return value is a list of the same length as \code{vars}, whose elements are plottable objects.
#'
#' @param obj an object of the class PSTR returned from some functions in the package. Note that the corresponding PSTR model must be estimated first.
#' @param vars a vector of column numbers or names (character strings) specifying which variables in the nonlinear part to use.
#' @param log_scale a 2-dim vector or scalar specifying whether to take log scale for the variables and the transition variable.
#' @param length.out a 2-dim vector or scalar of desired length (number of points) for the parameters. 20 by default.
#' @param color the color of the line.
#' @param size the size of the line.
#'
#' @return A list of plottable objects from the \code{ggplot2} (for curve) and/or \code{plotly} (for surface) package.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso Functions which return an object of the class PSTR and can be input into this function
#'
#' \code{\link{EstPSTR}}
#' @keywords utils
#'
#' @examples
#' \donttest{
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala','debta','cfa','sales'), iT=14) # create a new PSTR object
#'
#' # estimate the PSTR model first
#' pstr = EstPSTR(use=pstr, im=1, iq=1, useDelta=TRUE, par=c(.63,0), method='CG')
#'
#' # plot the curve and surfaces
#' ret = plot_response(obj=pstr, vars=1:4, log_scale = c(FALSE,TRUE), length.out=40)
#' attributes(ret)
#' ret$vala
#' ret$debta
#' }
#'
#' @export
plot_response <- function(obj, vars, log_scale=FALSE, length.out=20, color="blue", size=1.5)
{
  if(class(obj)!="PSTR")
    stop(simpleError("The argument 'obj' is not an object of class 'PSTR'"))
  if(is.null(obj$vg)) stop(simpleError("The PSTR model is not estimated yet."))

  if(length(length.out)==1) length.out = rep(length.out,2)
  if(length(log_scale)==1) log_scale = rep(log_scale,2)

  tvar = obj$mQ[,obj$iq]

  vy = seq(from=min(tvar), to=max(tvar), length.out=length.out[2])
  vg = fTF(vx=t(matrix(vy,length(vy),obj$imm)),gamma=obj$gamma,vc=obj$c)
  tvarname = obj$mQ_name[obj$iq]

  vyy = rep(obj$c, length.out[1])
  vgg = rep(.5, length.out[1])

  ftmp <- function(vu) (phi0 + phi1*vu[2])*vu[1]

  ret = list()
  for(vter in vars){
    vK = try(obj$mK[,vter],silent=T)
    if(class(vK)=='try-error' || length(vK)==0) next

    varname = obj$mK_name[vter]
    phi0 = obj$beta[paste0(varname,'_0')]
    phi1 = obj$beta[paste0(varname,'_1')]

    if(varname != tvarname){
      vx = seq(from=min(vK), to=max(vK), length.out=length.out[1])
      mz = t(matrix(apply(expand.grid(vx, vg),1,ftmp), nrow=length(vx)))
      vzz = c(apply(cbind(vx, vgg), 1, ftmp))

      tmpp = list(xaxis=list(title=paste0(varname,"_x")),
                  yaxis=list(title=paste0(tvarname,"_y")),zaxis=list(title="response"))
      if(log_scale[1]) tmpp$xaxis$type = "log"
      if(log_scale[2]) tmpp$yaxis$type = "log"

      tmp = add_surface(plot_ly(x=vx, y=vy, z=mz))

      tmp = tmp %>% layout(scene=tmpp)


      eval(parse(text=paste0("ret$",varname," = tmp")))
    }else{
      vz = c(apply(cbind(vy, vg),1,ftmp))

      tmp = ggplot(tibble(vy=vy,vz=vz), aes(x=vy,y=vz)) +
        labs(y="response", x=varname) + geom_line(color=color,size=size)
      if(log_scale[2]) tmp = tmp + scale_x_log10()

      eval(parse(text=paste0("ret$",varname," = tmp")))
    }

  }

  return(ret)
}



#' Plot the surface of the target function for the nonlinear least square estimation.
#'
#' This function plots the surface of the target function for the nonlinear least square estimation.
#' It is useful for finding the suitable initial value for the estimation.
#'
#' The funciton uses the \code{plotly} package to plot the 3-D surface of the target function for the nonlinear least square estimation.
#'
#' The function takes the PSTR object as one of the inputs. The user needs to give the number of switches \code{im}, and the transition variable \code{iq},
#' such that the target function values can be computed.
#'
#' The number of parameters to estimate in the nonlinear least square estimation is \code{1+im}, that is, one smoothness parameter and the \code{im} switching locations.
#' However, the 3-D plot is based on only two changing parameters with the others (if more than two parameters) constant. Thus, the user needs to input a vector \code{par},
#' which gives the values of the other parameters. Note that \code{par} should still be of length \code{1+im} with the order \eqn{\delta} (always use delta in this function),
#' \eqn{c_1}, ..., \eqn{c_m}.
#'
#' The user should give the vector \code{basedon} of length two, that shows which two parameters will be used to build the grid.
#' \code{basedon} gives the positions of the two parameters in \code{par}. Thus, the values in the positions \code{basedon} in \code{par} will not be used.
#'
#' \code{from}, \code{to} and \code{length.out} serve to build the grid for the two parameters.
#' These arguments must be of length two for the two parameters, respectively.
#' See the \code{seq} function for the details.
#'
#'
#' @param obj an object of the class PSTR returned from some functions in the package.
#' @param im specifies the number of switches in the transtion function. The default value is 1.
#' @param iq a column number (in \code{mQ}) or variable name specifying the transition variable to use.
#' @param par a vector of the values of the parameters. NULL by default, then it will be made automatically.
#' @param basedon an integer vector of length 2 specify which two parameters to use to build the grid.
#' @param from a vector of length 2 of the starting (minimal) values of the parameters.
#' @param to a vector of length 2 of the end (maximal) values of the parameters.
#' @param length.out a 2-dim vector or scalar of desired length (number of points) for the parameters. 40 by default.
#'
#' @return A plottable object from the \code{plotly} package.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso Functions which return an object of the class PSTR and can be input into this function
#'
#' \code{\link{NewPSTR}}, \code{\link{LinTest}}, \code{\link{WCB_LinTest}}, \code{\link{EstPSTR}}, \code{\link{EvalTest}}, \code{\link{WCB_TVTest}} and \code{\link{WCB_HETest}}
#'
#' @keywords utils
#'
#' @examples
#' \donttest{
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala'), iT=14) # create a new PSTR object
#'
#' # build the grid based on the first two parameters
#' ret = plot_target(obj=pstr,iq=1,basedon=c(1,2),from=c(log(1),6),
#'   to=c(log(18),10),length.out=c(40,40))
#' }
#' 
#' @export
plot_target <- function(obj,im=1,iq=NULL,par=NULL,basedon=c(1,2),from,to,length.out=40)
{
  if(class(obj)!="PSTR")
    stop(simpleError("The argument 'obj' is not an object of class 'PSTR'"))
  ret = NULL
  iT = obj$iT; iN = obj$iN

  # get the data here
  vY = obj$vY; vYb = obj$vYb
  mX = obj$mX; mXb = obj$mXb
  mK = obj$mK
  ik = ncol(mK)

  ftmp <- function(vx) return(vx - mean(vx))

  if(!is.null(iq)){
    if(im < 1) stop(simpleError("The number of switches is invalid."))

    pnames = paste0('c_',basedon-1)
    pnames[pnames=='c_0'] = "delta"

    vQ = obj$mQ[,iq]
    mQ = t(matrix(vQ,iT*iN,im))

    if(is.null(par)){
      tmp = unname(quantile(vQ, (1:im) / (im+1)))
      par = c(log(8/min(diff(c(0,tmp)))), tmp)
    }

    if(length(length.out)==1) length.out = rep(length.out,2)

    ret$x = seq(from=from[1], to=to[1], length.out=length.out[1])
    ret$y = seq(from=from[2], to=to[2], length.out=length.out[2])

    ret$com = expand.grid(ret$x, ret$y)

    ResiduleSumSquare <- function(vpp){
      vp = par
      vp[basedon] = vpp
      vg = fTF(vx=mQ,gamma=exp(vp[1]),vc=vp[2:length(vp)])
      mXX = mK * vg
      aXX = array(c(mXX), dim=c(iT,iN,ik))
      mXXb = cbind(mXb, matrix(c(apply(aXX,c(2,3),ftmp)), iT*iN, ik))
      tmp = chol2inv(chol(t(mXXb)%*%mXXb)) %*% t(mXXb) %*% vYb
      vE = c(vYb-mXXb%*%tmp)
      return(sum(vE*vE)/iT/iN)
    }

    ret$val = apply(ret$com,1,ResiduleSumSquare)
    ret$val = t(matrix(ret$val, nrow=length(ret$x)))

    ret = add_surface(plot_ly(x=ret$x, y=ret$y, z=ret$val))

    tmpp = list(xaxis=list(title=pnames[1]),
                yaxis=list(title=pnames[2]),zaxis=list(title="target"))
    ret = ret %>% layout(scene=tmpp)

    return(ret)
  }
  else stop(simpleError("Transition variable missing! Please specify iq."))

}

#' Plot the coefficients, the standard errors and the p-values against the transition variable.
#'
#' This function plots the curves of the coefficients, the standard errors and the p-values against the transition variable.
#'
#' The curves of the coefficients, the standard errors and the p-values against the transition variable are functions
#' \deqn{f_1(x) = \beta_{0j} + \beta_{1j}g(x ; \gamma, c)}
#' \deqn{f_2(x) = se(f_1(x))}
#' \deqn{f_3(x) = 1 - Prob\{ X < [f_1(x)/f_2(x)]^2 \} }
#' where \eqn{x} is a variable taking the position of the transition variable,
#' \eqn{se} stands for the cluster-robust and heteroskedasticity-consistent standard error of the estimate \eqn{f_1(x)} at \eqn{x},
#' \eqn{X} is a random variable following chi-square distribution with degrees of freedom one.
#' 
#' More than one variable can be put in \code{vars}.
#'
#' The return value is a list of the same length as \code{vars}, whose elements are plottable objects.
#'
#' @param obj an object of the class PSTR returned from some functions in the package. Note that the corresponding PSTR model must be estimated first.
#' @param vars a vector of column numbers or names (character strings) specifying which variables in the nonlinear part to use.
#' @param length.out a scalar of desired length (number of points) for building the x-axis. 100 by default.
#' @param color the color of the lines.
#' @param size the size of the lines.
#'
#' @return A list of plottable objects from the \code{ggplot2} package.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @seealso Functions which return an object of the class PSTR can be input into this function
#'
#' \code{\link{EstPSTR}}
#' 
#' @keywords utils
#'
#' @examples
#' \donttest{
#' pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#'     tvars=c('vala','debta','cfa','sales'), iT=14) # create a new PSTR object
#'
#' # estimate the PSTR model first
#' pstr = EstPSTR(use=pstr, im=1, iq=1, useDelta=TRUE, par=c(.63,0), method='CG')
#'
#' # plot the curve and surfaces
#' ret = plot_coefficients(pstr, vars=1:4, length.out=100, color="dodgerblue4", size=2)
#' ret[[1]]
#' ret[[1]] + ggplot2::scale_x_log10()
#' }
#'
#' @export
plot_coefficients <- function(obj, vars, length.out=100, color="blue", size=1.5)
{
  if(class(obj)!="PSTR")
    stop(simpleError("The argument 'obj' is not an object of class 'PSTR'"))
  if(is.null(obj$vg)) stop(simpleError("The PSTR model is not estimated yet."))
  
  tvar = obj$mQ[,obj$iq]
  tvarname = obj$mQ_name[obj$iq]
  
  pp = seq(from=min(tvar), to=max(tvar), length.out=length.out)
  mx = t(matrix(pp,length(pp),obj$imm))
  ratio = fTF(mx, obj$gamma, obj$c)
  
  tnames = names(obj$est)
  
  ret = list()
  for(vter in vars){
    vK = try(obj$mK[,vter],silent=T)
    if(class(vK)=='try-error' || length(vK)==0) next
    
    if(is.numeric(vter)) tchar = obj$mK_name[vter]
    else tchar = vter
    
    tmp = rep(0, length(obj$est))
    tmp[match(paste0(tchar,"_0"), tnames)] = 1
    kter = match(paste0(tchar,"_1"), tnames)
    
    bb = NULL; se = NULL
    for(jter in ratio){
      tmp[kter] = jter
      bb = c(bb, crossprod(tmp, obj$est))
      se = c(se, sqrt(t(tmp)%*%obj$cov%*%tmp))
    }
    pv = 1-pchisq((bb/se)**2,df=1)
    
    dtmp = data.frame(xx=rep(pp,3), yy=c(bb, se, pv),
                  gg=factor(c( rep("\u03b2", length(bb)), rep("s.e.", length(se)),
                               rep("p-val", length(pv))), levels=c("\u03b2", "s.e.", "p-val")) ) 
    
    tmp = with(dtmp, ggplot(dtmp, aes(x=xx,y=yy))) + geom_line(col=color,size=size) +
      with(dtmp,facet_grid(rows = vars(gg), scales = "free")) +
      geom_hline(data=with(dtmp, subset(dtmp, gg=='p-val')), aes(yintercept=.05), col='red', linetype=2) +
      labs(x=tvarname, y='', title=paste("coefficient",tchar))
    
    eval(parse(text=paste0("ret$",tchar," = tmp")))
    
  }
  
  return(ret)
}
