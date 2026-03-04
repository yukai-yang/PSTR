test_that("Full PSTR workflow runs without error (light version)", {
  
  skip_if_not_installed("cli")
  skip_if_not_installed("R6")
  skip_if_not_installed("tibble")
  skip_if_not_installed("knitr")
  
  # Only needed if your parallel backend lives in Suggests and is touched indirectly
  # (safe even if not strictly required when parallel = FALSE)
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  
  data("Hansen99", package = "PSTR")
  
  get_private <- function(obj, name) {
    obj$.__enclos_env__$private[[name]]
  }
  
  pstr <- NewPSTR(
    Hansen99,
    dep   = "inva",
    indep = 4:20,
    indep_k = c("vala", "debta", "cfa", "sales"),
    tvars = c("vala"),
    im = 1,
    iT = 14
  )
  expect_s3_class(pstr, "PSTR")
  
  expect_no_error(suppressMessages(pstr$LinTest()))
  expect_true(!is.null(pstr$test))
  
  expect_no_error(suppressMessages(
    WCB_LinTest(pstr, iB = 2, parallel = FALSE, cpus = 1)
  ))
  expect_true(!is.null(pstr$wcb_test))
  
  expect_no_error(capture.output(print(pstr, mode = "tests")))
  
  expect_no_error(suppressMessages(
    EstPSTR(pstr, im = 1, iq = 1, useDelta = TRUE, par = c(-0.462, 0), method = "CG")
  ))
  expect_true(!is.null(pstr$vg))
  expect_true(!is.null(pstr$beta))
  
  expect_no_error(suppressMessages(
    pstr$EstPSTR(im = 1, iq = 1, useDelta = TRUE, par = c(-0.462, 0), method = "CG")
  ))
  
  expect_no_error(capture.output(print(pstr, mode = "estimates")))
  
  vq_internal <- as.numeric(get_private(pstr, "mQ")[, 1])
  expect_true(length(vq_internal) == nrow(get_private(pstr, "mX")))
  
  expect_no_error(suppressMessages(EvalTest(pstr, vq = vq_internal)))
  expect_true(!is.null(pstr$tv) || !is.null(pstr$ht))
  
  expect_no_error(suppressMessages(
    WCB_TVTest(use = pstr, iB = 2, parallel = FALSE, cpus = 1)
  ))
  expect_true(!is.null(pstr$wcb_tv))
  
  expect_no_error(suppressMessages(
    WCB_HETest(use = pstr, vq = vq_internal, iB = 2, parallel = FALSE, cpus = 1)
  ))
  expect_true(!is.null(pstr$wcb_ht))
  
  expect_no_error(capture.output(print(pstr, mode = "evaluation")))
})