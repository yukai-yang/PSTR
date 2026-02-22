test_that("Full PSTR workflow runs without error (light version)", {
  
  skip_if_not_installed("cli")
  skip_if_not_installed("R6")
  skip_if_not_installed("tibble")
  
  data("Hansen99", package = "PSTR")
  
  # helper: read a private field from an R6 instance (testing only)
  get_private <- function(obj, name) {
    obj$.__enclos_env__$private[[name]]
  }
  
  # 1) create
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
  
  # 2) linearity tests
  expect_no_error(suppressMessages(pstr$LinTest()))
  expect_true(!is.null(pstr$test))
  
  expect_no_error(suppressMessages(WCB_LinTest(pstr, iB = 2, parallel = FALSE, cpus = 1)))
  expect_true(!is.null(pstr$wcb_test))
  
  expect_no_error(capture.output(print(pstr, mode = "tests")))
  
  # 3) estimate PSTR (wrapper)
  expect_no_error(suppressMessages(
    EstPSTR(pstr, im = 1, iq = 1, useDelta = TRUE, par = c(-0.462, 0), method = "CG")
  ))
  expect_true(!is.null(pstr$vg))
  expect_true(!is.null(pstr$beta))
  
  # 3b) direct method interface (DO NOT call with default iq=NULL)
  expect_no_error(suppressMessages(
    pstr$EstPSTR(im = 1, iq = 1, useDelta = TRUE, par = c(-0.462, 0), method = "CG")
  ))
  
  expect_no_error(capture.output(print(pstr, mode = "estimates")))
  
  # 4) evaluation tests
  # IMPORTANT: vq must align with NA-removed internal panel rows
  vq_internal <- as.numeric(get_private(pstr, "mQ")[, 1])
  expect_true(length(vq_internal) == nrow(get_private(pstr, "mX")))
  
  expect_no_error(suppressMessages(EvalTest(pstr, vq = vq_internal)))
  expect_true(!is.null(pstr$tv) || !is.null(pstr$ht))
  
  # 5) bootstrap evaluation (light)
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