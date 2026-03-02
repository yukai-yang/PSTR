# Internal helper: check future.apply availability
.pstr_need_future_apply <- function() {
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Parallel option requires package 'future.apply'. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Parallel option requires package 'future'. Please install it.", call. = FALSE)
  }
  invisible(TRUE)
}

# Internal helper: map in serial or parallel
# - X: list or vector of tasks
# - FUN: function applied to each element
# - parallel: TRUE/FALSE
# - workers: integer number of workers (optional)
# - seed: for reproducibility in parallel bootstrap
.pstr_lapply <- function(X, FUN, ..., parallel = FALSE, workers = NULL, seed = TRUE) {
  if (!parallel) {
    return(lapply(X, FUN, ...))
  }
  
  .pstr_need_future_apply()
  
  # Respect existing plan if user has set one, unless they asked us to choose workers
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  if (!is.null(workers)) {
    # Cross-platform safe choice: multisession
    future::plan(future::multisession, workers = workers)
  } else {
    # Do not overwrite user's plan
    # If no plan set, future will use default sequential; enforce multisession as a safe default
    # But only if current plan is sequential
    if (inherits(old_plan, "sequential")) {
      future::plan(future::multisession)
    }
  }
  
  future.apply::future_lapply(X, FUN, ..., future.seed = seed)
}

# Optional: parallel apply that returns simplified vector (if you used sfSapply patterns)
.pstr_vapply <- function(X, FUN, FUN.VALUE, ..., parallel = FALSE, workers = NULL, seed = TRUE) {
  if (!parallel) {
    return(vapply(X, FUN, FUN.VALUE = FUN.VALUE, ...))
  }
  out <- .pstr_lapply(X, FUN, ..., parallel = TRUE, workers = workers, seed = seed)
  vapply(out, identity, FUN.VALUE = FUN.VALUE)
}