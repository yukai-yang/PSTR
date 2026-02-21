#################################################################################
## package name: PSTR
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## May 2018
#################################################################################

#################################################################################
## utility functions
#################################################################################

vnum = "2.0.0"
packname = "(Green Panel)"

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


svd_pinv <- function(A, tol = 1e-10) {
  s <- svd(A)
  d_inv <- ifelse(s$d > tol * max(s$d), 1 / s$d, 0)
  s$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(s$u)
}


# Print the object of the class PSTR.
#
# This function prints the object of the class PSTR.
#
# @param format argument passed to knitr::kable.
# @param mode a vector of character strings specifying which results to print. It takes the values c('summary', 'tests', 'estimates', 'evaluation'). By default 'su' and 'e' which means all.
# @param digits integer indicating the number of decimal places (for the \code{round} function inside) to be used. Negative values are allowed (see \code{round}).
# @param ... further arguments passed to knitr::kable.
#
# @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
# @seealso Functions which return an object of the class PSTR:
#
# \code{\link{NewPSTR}}, \code{\link{LinTest}}, \code{\link{WCB_LinTest}}, \code{\link{EstPSTR}}, \code{\link{EvalTest}}, \code{\link{WCB_TVTest}} and \code{\link{WCB_HETest}}
# @keywords utils
#
# @examples
# pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
#     tvars=c('vala','debta','cfa','sales'), iT=14)
#
# pstr
# print(pstr, mode="tests", format="simple")
# print(pstr, mode="tests", format="pipe",caption="The test results")
# print(pstr, mode="tests", format = "latex", caption = "The test results")
#
PSTR$set("public", "print", function(format="simple", mode=c("summary"), digits=4, ...){
  cli::cli_h1(paste0("R package PSTR ", vnum, " ",packname))
  
  tmp = NULL
  for(iter in 1:length(mode)){
    tmp = c(tmp, grep(mode[iter], c("summary","tests","estimates","evaluation")))
  }
  tmp = unique(tmp)
  
  if(1 %in% tmp){ private$print_summary(format, ...); cat("\n")}
  
  print_message = character()
  
  if(2 %in% tmp){
    private$print_tests(format, digits, ...); cat("\n")
  }else{
    if(!is.null(self$test) || !is.null(self$wcb_test)){
      #code = '`print(obj, mode="tests")`'
      #cli::cli_alert_info("The results of the linearity tests are ready, run {code} to show the results."); cat("\n")
      print_message = c(print_message, 'The specification results are ready, run `print(obj, mode="tests")` to show the results.')
    }
  }
  
  if(3 %in% tmp){
    private$print_estimates(format, digits, ...); cat("\n")
  }else{
    if(!is.null(self$est)){
      #code = '`print(obj, mode="estimates")`'
      #cli::cli_alert_info("The estimation results are ready, run {code} to show the results."); cat("\n")
      print_message = c(print_message, 'The estimation results are ready, run `print(obj, mode="estimates")` to show the results.')
    }
  }
    
  if(4 %in% tmp){
    private$print_evaluation(format, digits, ...); cat("\n")
  }else{ # tv, ht, wcb_tv, wcb_ht
    if(!is.null(self$tv) || !is.null(self$ht) || !is.null(self$wcb_tv) || !is.null(self$wcb_ht)){
      #code = '`print(obj, mode="evaluation")`'
      #cli::cli_alert_info("The evaluation results are ready, run {code} to show the results."); cat("\n")
      print_message = c(print_message, 'The evaluation results are ready, run `print(obj, mode="evaluation")` to show the results.')
    }
  }
  
  if(length(tmp)==0){
    cli::cli_alert_info("The argument 'mode' only accepts the values:")
    cat0("  'summary', 'tests', 'estimates' or 'evaluation'.")
    cat0("  Incomplete words, such as 'su' or 'mm' for 'summary', are allowed.")
  }else{
    if (length(print_message) > 0) {
      cli::cli_bullets(setNames(print_message, rep("i", length(print_message))))
    }
  }
  
})


PSTR$set("private", "print_summary", function(...) {
  cli::cli_h2("Summary of the model")
  
  cli::cli_text(
    "The long format panel is {private$iT} × {private$iN} (time × individual)."
  )
  
  cli::cli_rule()
  
  cli::cli_h3("Dependent variable")
  cli::cli_text("{private$vY_name}")
  
  cli::cli_h3("Explanatory variables in the linear part ({length(private$mX_name)})")
  cli::cli_text("{paste(private$mX_name, collapse = ', ')}")
  
  cli::cli_h3("Explanatory variables in the non-linear part ({length(private$mK_name)})")
  cli::cli_text("{paste(private$mK_name, collapse = ', ')}")
  
  cli::cli_h3("Potential transition variable(s) to be tested ({length(private$mQ_name)})")
  cli::cli_text("{paste(private$mQ_name, collapse = ', ')}")
  
  invisible(self)
})


PSTR$set("private", "print_tests", function(format, digits, ...) {
  
  im <- private$im
  
  # helper: build two-row-per-m table for one transition variable
  build_tab <- function(test_list_for_one_q, wcb_mat_for_one_q = NULL, im) {
    
    out <- NULL
    
    for (m in 1:im) {
      
      ttmp <- test_list_for_one_q[[m]]
      
      # 1) stat row
      r1 <- c(
        m,
        ttmp$LM1_X,
        ttmp$LM1_F,
        ttmp$LM2_X,
        ttmp$LM2_F
      )
      
      # 2) p-val row
      r2 <- c(
        "p-val",
        ttmp$PV1_X,
        ttmp$PV1_F,
        ttmp$PV2_X,
        ttmp$PV2_F
      )
      
      # optional: add WB/WCB p-values (already p-values)
      if (!is.null(wcb_mat_for_one_q)) {
        wb  <- wcb_mat_for_one_q[m, 2]
        wcb <- wcb_mat_for_one_q[m, 3]
        r1 <- c(r1, wb, wcb)
        r2 <- c(r2, wb, wcb)
      }
      
      out <- rbind(out, r1, r2)
    }
    
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    
    if (is.null(wcb_mat_for_one_q)) {
      colnames(out) <- c("m", "LM_X", "LM_F", "HAC_X", "HAC_F")
    } else {
      colnames(out) <- c("m", "LM_X", "LM_F", "HAC_X", "HAC_F", "WB_PV", "WCB_PV")
    }
    
    # numeric conversion for numeric columns
    num_cols <- setdiff(colnames(out), "m")
    out[num_cols] <- lapply(out[num_cols], as.numeric)
    
    out
  }
  
  printed_any <- FALSE
  
  # ---------------------------
  # Main linearity tests
  # ---------------------------
  if (!is.null(self$test) || !is.null(self$wcb_test)) {
    
    printed_any <- TRUE
    cli::cli_h2("Results of the linearity (homogeneity) tests")
    
    # how many transition variables are we printing?
    nQ <- max(length(self$test), length(self$wcb_test))
    
    for (iter in seq_len(nQ)) {
      
      cli::cli_h3(
        paste0("LM tests based on transition variable '",
               private$mQ_name[iter], "'")
      )
      
      # main tests exist?
      if (!is.null(self$test) && iter <= length(self$test)) {
        wcb_here <- NULL
        if (!is.null(self$wcb_test) && iter <= length(self$wcb_test)) {
          wcb_here <- self$wcb_test[[iter]]
        }
        tab <- build_tab(self$test[[iter]], wcb_here, im)
        print(knitr::kable(tab, format = format, digits = digits, row.names = FALSE, ...))
        
      } else if (!is.null(self$wcb_test) && iter <= length(self$wcb_test)) {
        
        # edge case: only wcb exists but no asymptotic tests
        # show a compact table (still two rows per m so your layout remains consistent)
        w <- self$wcb_test[[iter]]
        out <- NULL
        for (m in 1:im) {
          wb  <- w[m, 2]
          wcb <- w[m, 3]
          out <- rbind(out,
                       c(m,     NA, NA, NA, NA, wb, wcb),
                       c("p-val", NA, NA, NA, NA, wb, wcb))
        }
        out <- as.data.frame(out, stringsAsFactors = FALSE)
        colnames(out) <- c("m", "LM_X", "LM_F", "HAC_X", "HAC_F", "WB_PV", "WCB_PV")
        num_cols <- setdiff(colnames(out), "m")
        out[num_cols] <- lapply(out[num_cols], as.numeric)
        print(knitr::kable(out, format = format, digits = digits, row.names = FALSE, ...))
      }
    }
  }
  
  # ---------------------------
  # Sequence tests
  # ---------------------------
  if (!is.null(self$sqtest) || !is.null(self$wcb_sqtest)) {
    
    printed_any <- TRUE
    cli::cli_h2("Sequence of homogeneity tests for selecting number of switches 'm'")
    
    nQ <- max(length(self$sqtest), length(self$wcb_sqtest))
    
    for (iter in seq_len(nQ)) {
      
      cli::cli_h3(
        paste0("LM tests based on transition variable '",
               private$mQ_name[iter], "'")
      )
      
      if (!is.null(self$sqtest) && iter <= length(self$sqtest)) {
        wcb_here <- NULL
        if (!is.null(self$wcb_sqtest) && iter <= length(self$wcb_sqtest)) {
          wcb_here <- self$wcb_sqtest[[iter]]
        }
        tab <- build_tab(self$sqtest[[iter]], wcb_here, im)
        print(knitr::kable(tab, format = format, digits = digits, row.names = FALSE, ...))
        
      } else if (!is.null(self$wcb_sqtest) && iter <= length(self$wcb_sqtest)) {
        
        w <- self$wcb_sqtest[[iter]]
        out <- NULL
        for (m in 1:im) {
          wb  <- w[m, 2]
          wcb <- w[m, 3]
          out <- rbind(out,
                       c(m,     NA, NA, NA, NA, wb, wcb),
                       c("p-val", NA, NA, NA, NA, wb, wcb))
        }
        out <- as.data.frame(out, stringsAsFactors = FALSE)
        colnames(out) <- c("m", "LM_X", "LM_F", "HAC_X", "HAC_F", "WB_PV", "WCB_PV")
        num_cols <- setdiff(colnames(out), "m")
        out[num_cols] <- lapply(out[num_cols], as.numeric)
        print(knitr::kable(out, format = format, digits = digits, row.names = FALSE, ...))
      }
    }
  }
  
  if (!printed_any) {
    code <- "`PSTR::LinTest()`"
    cli::cli_alert_warning("The linearity tests have not been conducted yet, run {code}.")
  }
  
  invisible(self)
})


PSTR$set("private", "print_estimates", function(format, digits, ..., max_cols = NULL) {
  
  # nothing estimated yet
  if (is.null(self$est) || is.null(self$se)) {
    code <- "`PSTR::EstPSTR()`"
    cli::cli_alert_warning("The model has not been estimated yet, run {code}.")
    return(invisible(self))
  }
  
  # helper: chunked wide table with Est / s.e. / t-ratio
  print_chunked_coef_table <- function(est, se, title = NULL) {
    
    if (!is.null(title)) cli::cli_h3(title)
    
    nm <- names(est)
    if (is.null(nm) || anyNA(nm) || any(nm == "")) {
      nm <- paste0("p", seq_along(est))
    }
    
    # t-ratio
    tr <- est / se
    
    # choose columns per chunk
    w <- cli::console_width()
    
    # crude but stable width estimate per parameter column
    name_w <- max(nchar(nm), na.rm = TRUE)
    num_w  <- max(10L, digits + 6L)  # sign + integer + dot + decimals
    col_w  <- max(name_w, num_w) + 2L
    
    # room for row names + separators; keep a safety margin
    if (is.null(max_cols)) {
      max_cols <- floor((w - 12L) / col_w)
      max_cols <- max(1L, min(max_cols, length(est)))
    } else {
      max_cols <- max(1L, min(as.integer(max_cols), length(est)))
    }
    
    idx_list <- split(seq_along(est), ceiling(seq_along(est) / max_cols))
    
    for (idx in idx_list) {
      tab <- rbind(
        Est     = est[idx],
        `s.e.`  = se[idx],
        `t-ratio` = tr[idx]
      )
      tab <- signif(tab, digits)
      
      # keep original parameter names as column names
      colnames(tab) <- nm[idx]
      
      print(knitr::kable(tab, format = format, ...))
      cat("\n")
    }
  }
  
  # nonlinear PSTR estimated
  if (!is.null(private$iq)) {
    
    cli::cli_h2("Results of the PSTR estimation:")
    cli::cli_alert_info("Transition variable '{private$mQ_name[private$iq]}' is used in the estimation.")
    
    kx <- length(private$mX_name)
    
    # linear part (beta_0)
    print_chunked_coef_table(
      est = self$beta[1:kx],
      se  = self$se[1:kx],
      title = "Parameter estimates in the linear part (first extreme regime)"
    )
    
    # non-linear part (beta_1)
    print_chunked_coef_table(
      est = self$beta[(kx + 1):length(self$beta)],
      se  = self$se[(kx + 1):length(self$beta)],
      title = "Parameter estimates in the non-linear part"
    )
    
    # second extreme regime (beta_0 + beta_1) if available
    if (!is.null(self$mbeta) && !is.null(self$mse)) {
      est2 <- self$mbeta
      se2  <- self$mse
      nm2  <- names(est2)
      if (!is.null(nm2)) names(est2) <- paste0(nm2, "_{0+1}")
      
      print_chunked_coef_table(
        est = est2,
        se  = se2,
        title = "Parameter estimates in the second extreme regime"
      )
    }
    
    # nonlinear parameters (gamma and c's)
    est3 <- self$est[(length(self$beta) + 1):length(self$est)]
    se3  <- self$se[(length(self$beta) + 1):length(self$se)]
    
    print_chunked_coef_table(
      est = est3,
      se  = se3,
      title = "Non-linear parameter estimates"
    )
    
    cli::cli_alert_info(
      "Estimated standard deviation of the residuals is {signif(sqrt(self$s2), digits)}."
    )
    
  } else {
    
    # linear FE estimated
    cli::cli_h2("A linear panel regression with fixed effects is estimated.")
    
    print_chunked_coef_table(
      est = self$est,
      se  = self$se,
      title = "Parameter estimates"
    )
    
    cli::cli_alert_info(
      "Estimated standard deviation of the residuals is {signif(sqrt(self$s2), digits)}."
    )
  }
  
  invisible(self)
})


PSTR$set("private", "print_evaluation", function(format, digits, ...) {
  
  has_eval <- (!is.null(self$tv) || !is.null(self$wcb_tv) ||
                 !is.null(self$ht) || !is.null(self$wcb_ht))
  
  if (!has_eval) {
    code <- "`PSTR::EvalTest()`"
    code2 <- "`PSTR::WCB_TVTest()` / `PSTR::WCB_HETest()`"
    cli::cli_alert_warning("No evaluation test results found. Run {code} and/or {code2}.")
    return(invisible(self))
  }
  
  cli::cli_h2("Results of the evaluation tests")
  
  # helper: make "stat + p-val" two-line rows per m, without adding extra columns
  mk_stat_pval_table <- function(obj_list, im) {
    out <- NULL
    for (m in 1:im) {
      ttmp <- obj_list[[m]]
      
      stat_row <- c(
        m,
        ttmp$LM1_X, ttmp$LM1_F, ttmp$LM2_X, ttmp$LM2_F
      )
      pval_row <- c(
        "p-val",
        ttmp$PV1_X, ttmp$PV1_F, ttmp$PV2_X, ttmp$PV2_F
      )
      
      out <- rbind(out, stat_row, pval_row)
    }
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    colnames(out) <- c("m", "LM_X", "LM_F", "HAC_X", "HAC_F")
    
    # numeric formatting for all but first column (which is m / p-val)
    for (j in 2:ncol(out)) out[[j]] <- as.numeric(out[[j]])
    out[, 2:ncol(out)] <- signif(out[, 2:ncol(out), drop = FALSE], digits)
    out
  }
  
  mk_boot_table <- function(mat) {
    # mat is self$wcb_tv or self$wcb_ht; columns 2:3 are p-values
    tmp <- mat[, 2:3, drop = FALSE]
    colnames(tmp) <- c("WB_PV", "WCB_PV")
    
    im <- nrow(tmp)
    out <- NULL
    for (m in 1:im) {
      stat_row <- c(m, tmp[m, 1], tmp[m, 2])
      pval_row <- c("p-val", tmp[m, 1], tmp[m, 2])
      out <- rbind(out, stat_row, pval_row)
    }
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    colnames(out) <- c("m", "WB_PV", "WCB_PV")
    out[, 2:3] <- signif(as.data.frame(lapply(out[, 2:3, drop = FALSE], as.numeric)), digits)
    out
  }
  
  # --- Parameter constancy (time-varying) ---
  if (!is.null(self$tv)) {
    cli::cli_h3("Parameter constancy test")
    im <- length(self$tv)
    tab <- mk_stat_pval_table(self$tv, im)
    print(knitr::kable(tab, format = format, row.names = FALSE, ...))
    cat("\n")
  }
  
  if (!is.null(self$wcb_tv)) {
    cli::cli_h3("WB and WCB parameter constancy test")
    tab <- mk_boot_table(self$wcb_tv)
    print(knitr::kable(tab, format = format, row.names = FALSE, ...))
    cat("\n")
  }
  
  # --- No remaining nonlinearity / heterogeneity ---
  if (!is.null(self$ht)) {
    cli::cli_h3("No remaining nonlinearity (heterogeneity) test")
    im <- length(self$ht)
    tab <- mk_stat_pval_table(self$ht, im)
    print(knitr::kable(tab, format = format, row.names = FALSE, ...))
    cat("\n")
  }
  
  if (!is.null(self$wcb_ht)) {
    cli::cli_h3("WB and WCB no remaining nonlinearity (heterogeneity) test")
    tab <- mk_boot_table(self$wcb_ht)
    print(knitr::kable(tab, format = format, row.names = FALSE, ...))
    cat("\n")
  }
  
  invisible(self)
})


#' Plot the transition function of an estimated PSTR model
#'
#' This function plots the estimated transition function
#' \eqn{g(q;\gamma,c)} of a fitted PSTR model.
#'
#' Observed transition values are displayed together with the
#' fitted transition curve. For models with multiple switches,
#' multiple curves are shown.
#'
#' In addition to the exported function
#' \code{plot_transition(obj = ...)}, the same functionality is
#' available as an R6 method via \code{obj$plot_transition(...)}.
#'
#' @param obj An object of class \code{"PSTR"}.
#' @param size Point size.
#' @param color Point colour.
#' @param xlim Optional numeric vector of length 2 specifying x-axis limits.
#' @param ylim Optional numeric vector of length 2 specifying y-axis limits.
#' @param fill Optional colour for highlighting the support of observed q.
#' @param alpha Transparency level for points and shading.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \donttest{
#' pstr <- NewPSTR(Hansen99, dep = "inva", indep = 4:20,
#'                 indep_k = c("vala","debta","cfa","sales"),
#'                 tvars = c("vala"), iT = 14)
#'
#' pstr <- EstPSTR(use = pstr, im = 1, iq = 1,
#'                 useDelta = TRUE, par = c(.63,0), method = "CG")
#'
#' # Exported function
#' plot_transition(pstr)
#'
#' # R6 method
#' pstr$plot_transition()
#' }
#'
#' @export
plot_transition <- function(obj,
                            size = 1.5,
                            color = "blue",
                            xlim = NULL,
                            ylim = NULL,
                            fill = NULL,
                            alpha = NULL) {
  
  if (!inherits(obj, "PSTR"))
    stop(simpleError("The argument 'obj' is not an object of class 'PSTR'"))
  
  obj$plot_transition(size = size,
                      color = color,
                      xlim = xlim,
                      ylim = ylim,
                      fill = fill,
                      alpha = alpha)
}


PSTR$set("public", "plot_transition", function(size = 1.5,
                                               color = "blue",
                                               xlim = NULL,
                                               ylim = NULL,
                                               fill = NULL,
                                               alpha = NULL) {
  
  if (is.null(self$vg))
    stop(simpleError("The PSTR model is not estimated yet."))
  
  iq  <- private$iq
  imm <- private$imm
  mQ  <- private$mQ
  qname <- private$mQ_name
  
  if (is.null(iq))
    stop(simpleError("No transition variable selected."))
  
  qq <- mQ[, iq]
  
  if (is.null(alpha))
    alpha <- 0.2
  
  if (is.null(ylim))
    ylim <- c(0, 1)
  
  # Default xlim
  if (is.null(xlim)) {
    g  <- self$gamma
    cc <- self$c
    
    if (is.null(g) || g <= 0 || any(!is.finite(cc))) {
      xlim <- range(qq, na.rm = TRUE)
    } else {
      z_lo <- -log(1 / 0.002472623 - 1) / g
      z_hi <- -log(1 / 0.9975274 - 1) / g
      cc_rng <- range(cc, na.rm = TRUE)
      xlim <- range(c(qq, cc_rng + z_lo, cc_rng + z_hi), na.rm = TRUE)
    }
  }
  
  vx <- seq(xlim[1], xlim[2], length.out = 1001)
  
  # Build curve
  if (imm == 1L) {
    
    vy <- fTF(vx, self$gamma, self$c)
    
    curve_df <- data.frame(
      x = vx,
      y = as.numeric(vy),
      m = 1L
    )
    
  } else {
    
    mx <- matrix(rep(vx, times = imm),
                 nrow = imm,
                 byrow = TRUE)
    
    vy <- fTF(mx, self$gamma, self$c)
    
    if (is.matrix(vy)) {
      
      curve_df <- data.frame(
        x = rep(vx, times = imm),
        y = as.numeric(t(vy)),
        m = rep(seq_len(imm), each = length(vx))
      )
      
    } else {
      
      curve_df <- data.frame(
        x = vx,
        y = as.numeric(vy),
        m = 1L
      )
      
    }
  }
  
  # Observed points
  vg <- self$vg
  
  if (is.vector(vg)) {
    
    point_df <- data.frame(
      x = qq,
      y = as.numeric(vg)
    )
    
  } else {
    
    point_df <- data.frame(
      x = rep(qq, times = nrow(vg)),
      y = as.numeric(t(vg))
    )
  }
  
  xlab <- if (!is.null(qname)) qname[iq] else paste0("q[", iq, "]")
  
  p <- ggplot2::ggplot() +
    ggplot2::xlim(xlim) +
    ggplot2::ylim(ylim)
  
  if (!is.null(fill)) {
    p <- p + ggplot2::geom_rect(
      ggplot2::aes(
        xmin = min(qq, na.rm = TRUE),
        xmax = max(qq, na.rm = TRUE),
        ymin = 0,
        ymax = 1
      ),
      fill = fill,
      alpha = alpha / 2
    )
  }
  
  p <- p +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = x, y = y, group = m),
      colour = "red"
    ) +
    ggplot2::geom_rug(
      data = point_df,
      ggplot2::aes(x = x),
      sides = "b",
      colour = color
    ) +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = x, y = y),
      colour = color,
      size = size,
      alpha = alpha
    ) +
    ggplot2::labs(
      x = xlab,
      y = "transition function"
    )
  
  p
})


#' Plot the expected response against selected variables
#'
#' This function plots the effect-adjusted expected response for selected
#' nonlinear variables in a PSTR model.
#'
#' If the selected variable differs from the transition variable,
#' a 3-D surface of
#' \deqn{(\beta_{k,0} + \beta_{k,1} g(q;\gamma,c)) z_{k}}
#' is plotted against \eqn{z_k} and the transition variable.
#'
#' If the selected variable coincides with the transition variable,
#' a curve is plotted instead.
#'
#' In addition to the exported function
#' \code{plot_response(obj = ...)}, the same functionality is available
#' as an R6 method via \code{obj$plot_response(...)}.
#'
#' @param obj An object of class \code{"PSTR"}.
#' @param vars Integer vector of column indices from the nonlinear part.
#' @param log_scale Logical scalar or length-2 vector indicating whether
#'   to use log scale for the regressor and/or transition variable.
#' @param length.out Scalar or length-2 numeric vector controlling grid size.
#' @param color Line colour.
#' @param size Line width.
#'
#' @return A named list of \code{ggplot2} (curve) and/or
#'   \code{plotly} (surface) objects.
#'
#' @examples
#' \donttest{
#' pstr <- NewPSTR(Hansen99, dep = "inva", indep = 4:20,
#'                 indep_k = c("vala","debta","cfa","sales"),
#'                 tvars = c("vala","debta","cfa","sales"), iT = 14)
#'
#' pstr <- EstPSTR(use = pstr, im = 1, iq = 1,
#'                 useDelta = TRUE, par = c(.63,0), method = "CG")
#'
#' # Exported interface
#' ret <- plot_response(pstr, vars = 1:4)
#'
#' # R6 method
#' ret2 <- pstr$plot_response(vars = 1:4)
#' }
#'
#' @export
plot_response <- function(obj, vars,
                          log_scale = FALSE,
                          length.out = 20,
                          color = "blue",
                          size = 1.5) {
  
  if (!inherits(obj, "PSTR"))
    stop(simpleError("The argument 'obj' is not an object of class 'PSTR'"))
  
  obj$plot_response(vars = vars,
                    log_scale = log_scale,
                    length.out = length.out,
                    color = color,
                    size = size)
}


PSTR$set("public", "plot_response", function(vars,
                                             log_scale = FALSE,
                                             length.out = 20,
                                             color = "blue",
                                             size = 1.5) {
  
  if (is.null(self$vg))
    stop(simpleError("The PSTR model is not estimated yet."))
  
  if (length(length.out) == 1)
    length.out <- rep(length.out, 2)
  
  if (length(log_scale) == 1)
    log_scale <- rep(log_scale, 2)
  
  iq  <- private$iq
  imm <- private$imm
  mQ  <- private$mQ
  mK  <- private$mK
  
  if (is.null(iq))
    stop(simpleError("No transition variable selected."))
  
  tvar     <- mQ[, iq]
  tvarname <- private$mQ_name[iq]
  mK_name  <- private$mK_name
  
  vy <- seq(min(tvar), max(tvar), length.out = length.out[2])
  mx <- matrix(rep(vy, times = imm), nrow = imm, byrow = TRUE)
  
  vg_grid <- fTF(vx = mx, gamma = self$gamma, vc = self$c)
  vg1 <- if (is.matrix(vg_grid)) vg_grid[1, ] else vg_grid
  
  ret <- list()
  
  for (vter in vars) {
    
    if (vter < 1 || vter > ncol(mK)) next
    
    varname <- mK_name[vter]
    vK <- mK[, vter]
    
    phi0 <- self$beta[paste0(varname, "_0")]
    phi1 <- self$beta[paste0(varname, "_1")]
    
    ftmp <- function(vu) (phi0 + phi1 * vu[2]) * vu[1]
    
    if (varname != tvarname) {
      
      vx <- seq(min(vK), max(vK), length.out = length.out[1])
      mz <- t(matrix(apply(expand.grid(vx, vg1), 1, ftmp),
                     nrow = length(vx)))
      
      scene <- list(
        xaxis = list(title = varname),
        yaxis = list(title = tvarname),
        zaxis = list(title = "response")
      )
      
      if (log_scale[1]) scene$xaxis$type <- "log"
      if (log_scale[2]) scene$yaxis$type <- "log"
      
      p <- plotly::plot_ly(x = vx, y = vy, z = mz) |>
        plotly::add_surface() |>
        plotly::layout(scene = scene)
      
      ret[[varname]] <- p
      
    } else {
      
      vz <- as.numeric(apply(cbind(vy, vg1), 1, ftmp))
      
      p <- ggplot2::ggplot(
        tibble::tibble(vy = vy, vz = vz),
        ggplot2::aes(x = vy, y = vz)
      ) +
        ggplot2::geom_line(colour = color, linewidth = size) +
        ggplot2::labs(x = varname, y = "response")
      
      if (log_scale[2])
        p <- p + ggplot2::scale_x_log10()
      
      ret[[varname]] <- p
    }
  }
  
  ret
})



#' Plot the surface of the target function for nonlinear least squares estimation
#'
#' This function plots a 3-D surface of the nonlinear least squares (NLS) target function
#' used in estimating a PSTR model. It is mainly intended as a diagnostic tool for choosing
#' reasonable initial values for the nonlinear parameters.
#'
#' The target function is evaluated on a two-dimensional grid over two selected parameters,
#' while all other nonlinear parameters are held fixed at values provided by \code{par}.
#' The nonlinear parameter vector is always ordered as
#' \eqn{\delta, c_1, \ldots, c_m}, where \eqn{\gamma = \exp(\delta)} and \eqn{m = im}.
#'
#' In addition to the exported function \code{plot_target(obj = ...)}, the same functionality
#' is available as an R6 method via \code{obj$plot_target(...)}.
#'
#' @param obj An object of class \code{"PSTR"}.
#' @param im Integer. The number of switches \eqn{m} in the transition function used to
#'   construct the target function surface. Default is \code{1}.
#' @param iq Either an integer index (a column number in the transition-variable matrix)
#'   or a character string (a transition-variable name) specifying which transition variable
#'   is used when computing the target function.
#' @param par Numeric vector of length \code{1 + im} giving fixed values of the nonlinear
#'   parameters in the order \code{c(delta, c_1, ..., c_m)}. If \code{NULL}, the function
#'   constructs a default initial vector from quantiles of the chosen transition variable.
#' @param basedon Integer vector of length \code{2}. Positions in \code{par} indicating the
#'   two parameters to vary on the grid. The values at these positions in \code{par} are
#'   ignored and replaced by grid values.
#'   Use \code{1} for \eqn{\delta}, and \code{2:(im+1)} for \eqn{c_1, \ldots, c_m}.
#' @param from Numeric vector of length \code{2}. Lower bounds (start values) for the two
#'   grid parameters specified by \code{basedon}.
#' @param to Numeric vector of length \code{2}. Upper bounds (end values) for the two
#'   grid parameters specified by \code{basedon}.
#' @param length.out Either a scalar or a numeric vector of length \code{2}. If scalar, the
#'   same number of grid points is used for both dimensions. If length 2, the first element
#'   controls the grid resolution for the first parameter in \code{basedon}, and the second
#'   for the second parameter.
#'
#' @return A \code{plotly} object representing a 3-D surface plot of the target function
#'   values evaluated on the specified parameter grid.
#'
#' @seealso \code{\link{NewPSTR}}, \code{\link{LinTest}}, \code{\link{WCB_LinTest}},
#'   \code{\link{EstPSTR}}, \code{\link{EvalTest}}, \code{\link{WCB_TVTest}},
#'   \code{\link{WCB_HETest}}
#'
#' @examples
#' \donttest{
#' pstr <- NewPSTR(Hansen99, dep = "inva", indep = 4:20,
#'                indep_k = c("vala", "debta", "cfa", "sales"),
#'                tvars = c("vala"), iT = 14)
#'
#' # 1) Exported function interface
#' ret <- plot_target(obj = pstr, iq = 1, basedon = c(1, 2),
#'                    from = c(log(1), 6), to = c(log(18), 10),
#'                    length.out = c(40, 40))
#'
#' # 2) R6 method interface
#' ret2 <- pstr$plot_target(iq = 1, basedon = c(1, 2),
#'                          from = c(log(1), 6), to = c(log(18), 10),
#'                          length.out = c(40, 40))
#' }
#'
#' @export
plot_target <- function(obj, im = 1, iq = NULL, par = NULL, basedon = c(1, 2),
                        from, to, length.out = 40) {
  if (!inherits(obj, "PSTR"))
    stop(simpleError("The argument 'obj' is not an object of class 'PSTR'"))
  obj$plot_target(im = im, iq = iq, par = par, basedon = basedon,
                  from = from, to = to, length.out = length.out)
}

# 1) R6 method (inside class)
PSTR$set("public", "plot_target", function(im = 1, iq = NULL, par = NULL,
                                           basedon = c(1, 2), from, to, length.out = 40) {
  if (is.null(iq)) stop(simpleError("Transition variable missing! Please specify iq."))
  if (im < 1) stop(simpleError("The number of switches is invalid."))
  
  # allow iq as name
  if (!is.numeric(iq)) {
    iq <- which(private$mQ_name == iq)
    if (length(iq) != 1) stop(simpleError("Invalid iq name (not found or not unique)."))
  }
  
  iT <- private$iT
  iN <- private$iN
  
  vYb <- private$vYb
  mXb <- private$mXb
  mK  <- private$mK
  ik  <- ncol(mK)
  
  vQ <- private$mQ[, iq]
  mQ <- t(matrix(vQ, iT * iN, im))
  
  ftmp <- function(vx) vx - mean(vx)
  
  pnames <- paste0("c_", basedon - 1)
  pnames[pnames == "c_0"] <- "delta"
  
  if (is.null(par)) {
    tmp <- unname(quantile(vQ, (1:im) / (im + 1)))
    par <- c(log(8 / min(diff(c(0, tmp)))), tmp)
  }
  
  if (length(length.out) == 1) length.out <- rep(length.out, 2)
  
  x <- seq(from = from[1], to = to[1], length.out = length.out[1])
  y <- seq(from = from[2], to = to[2], length.out = length.out[2])
  com <- expand.grid(x, y)
  
  ResiduleSumSquare <- function(vpp) {
    vp <- par
    vp[basedon] <- vpp
    vg <- fTF(vx = mQ, gamma = exp(vp[1]), vc = vp[2:length(vp)])
    
    mXX <- mK * vg
    aXX <- array(c(mXX), dim = c(iT, iN, ik))
    mXXb <- cbind(mXb, matrix(c(apply(aXX, c(2, 3), ftmp)), iT * iN, ik))
    
    tmp <- chol2inv(chol(t(mXXb) %*% mXXb)) %*% t(mXXb) %*% vYb
    vE <- c(vYb - mXXb %*% tmp)
    sum(vE * vE) / iT / iN
  }
  
  val <- apply(com, 1, ResiduleSumSquare)
  val <- t(matrix(val, nrow = length(x)))
  
  ret <- plotly::add_surface(plotly::plot_ly(x = x, y = y, z = val))
  tmpp <- list(
    xaxis = list(title = pnames[1]),
    yaxis = list(title = pnames[2]),
    zaxis = list(title = "target")
  )
  ret <- ret %>% plotly::layout(scene = tmpp)
  
  ret
})



#' Plot coefficients, standard errors, and p-values against the transition variable
#'
#' This function plots three curves against the transition variable:
#' the coefficient function, its standard error, and the corresponding p-value.
#'
#' For each selected variable \eqn{j}, the curves are
#' \deqn{f_1(x) = \beta_{0j} + \beta_{1j} g(x;\gamma,c)}
#' \deqn{f_2(x) = se\{f_1(x)\}}
#' \deqn{f_3(x) = 1 - \Pr\left\{X < \left[f_1(x)/f_2(x)\right]^2\right\}}
#' where \eqn{X} follows a chi-square distribution with one degree of freedom.
#'
#' In addition to the exported function \code{plot_coefficients(obj = ...)},
#' the same functionality is available as an R6 method via
#' \code{obj$plot_coefficients(...)}.
#'
#' @param obj An object of class \code{"PSTR"}.
#' @param vars A vector of column indices or variable names from the nonlinear part.
#' @param length.out Number of grid points over the transition variable.
#' @param color Line colour.
#' @param size Line width.
#'
#' @return A named list of \code{ggplot2} objects.
#'
#' @examples
#' \donttest{
#' pstr <- NewPSTR(Hansen99, dep = "inva", indep = 4:20,
#'                 indep_k = c("vala","debta","cfa","sales"),
#'                 tvars = c("vala","debta","cfa","sales"), iT = 14)
#'
#' pstr <- EstPSTR(use = pstr, im = 1, iq = 1,
#'                 useDelta = TRUE, par = c(.63,0), method = "CG")
#'
#' # Exported function
#' ret <- plot_coefficients(pstr, vars = 1:4)
#'
#' # R6 method
#' ret2 <- pstr$plot_coefficients(vars = 1:4)
#' }
#'
#' @export
plot_coefficients <- function(obj, vars, length.out = 100,
                              color = "blue", size = 1.5) {
  if (!inherits(obj, "PSTR"))
    stop(simpleError("The argument 'obj' is not an object of class 'PSTR'"))
  
  obj$plot_coefficients(vars = vars,
                        length.out = length.out,
                        color = color,
                        size = size)
}

PSTR$set("public", "plot_coefficients", function(vars, length.out = 100, color = "blue", size = 1.5) {
  
  if (is.null(self$vg))
    stop(simpleError("The PSTR model is not estimated yet."))
  if (is.null(self$est) || is.null(self$cov))
    stop(simpleError("Estimates or covariance matrix are not available."))
  
  iq <- private$iq
  if (is.null(iq))
    stop(simpleError("No transition variable selected."))
  
  mQ <- private$mQ
  tvarname <- private$mQ_name[iq]
  tvar <- mQ[, iq]
  
  imm <- private$imm
  mK <- private$mK
  mK_name <- private$mK_name
  
  if (length.out < 2)
    stop(simpleError("length.out must be >= 2."))
  
  pp <- seq(min(tvar), max(tvar), length.out = length.out)
  
  mx <- matrix(rep(pp, times = imm), nrow = imm, byrow = TRUE)
  ratio_grid <- fTF(mx, self$gamma, self$c)
  ratio <- if (is.matrix(ratio_grid)) ratio_grid[1, ] else ratio_grid
  
  tnames <- names(self$est)
  ret <- list()
  
  for (vter in vars) {
    
    if (is.numeric(vter)) {
      if (vter < 1 || vter > ncol(mK)) next
      tchar <- mK_name[vter]
    } else {
      tchar <- as.character(vter)
    }
    
    idx0 <- match(paste0(tchar, "_0"), tnames)
    idx1 <- match(paste0(tchar, "_1"), tnames)
    if (is.na(idx0) || is.na(idx1)) next
    
    w <- rep(0, length(self$est))
    w[idx0] <- 1
    
    bb <- numeric(length(ratio))
    se <- numeric(length(ratio))
    
    for (i in seq_along(ratio)) {
      w[idx1] <- ratio[i]
      bb[i] <- crossprod(w, self$est)
      se[i] <- sqrt(t(w) %*% self$cov %*% w)
    }
    
    pv <- 1 - stats::pchisq((bb / se)^2, df = 1)
    
    dtmp <- data.frame(
      xx = rep(pp, 3),
      yy = c(bb, se, pv),
      gg = factor(
        c(rep("\u03b2", length(bb)),
          rep("s.e.", length(se)),
          rep("p-val", length(pv))),
        levels = c("\u03b2", "s.e.", "p-val")
      )
    )
    
    p <- ggplot2::ggplot(dtmp, ggplot2::aes(x = xx, y = yy)) +
      ggplot2::geom_line(colour = color, linewidth = size) +
      ggplot2::facet_grid(rows = ggplot2::vars(gg), scales = "free") +
      ggplot2::geom_hline(
        data = subset(dtmp, gg == "p-val"),
        ggplot2::aes(yintercept = 0.05),
        colour = "red", linetype = 2
      ) +
      ggplot2::labs(x = tvarname,
                    y = "",
                    title = paste("coefficient", tchar))
    
    ret[[tchar]] <- p
  }
  
  ret
})