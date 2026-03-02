# ------------------------------------------------------------------------------
# Parallel helpers for PSTR (internal)
#
# Goal:
#   Replace snowfall-based parallelism with a future.apply backend while keeping
#   package code clean and consistent.
#
# Design principles:
#   1) Cross-platform safety: use multisession (works on Windows/macOS/Linux).
#   2) Do not force users into parallel: default is serial unless parallel=TRUE.
#   3) Be a good citizen: if we temporarily change the future plan or options,
#      restore them on exit so we do not pollute the user's session.
#   4) Guard rails:
#        - Cap 'workers' to availableCores() to avoid accidental overload.
#        - Use future.globals.maxSize as a "fuse" to fail early when huge objects
#          would otherwise be exported to workers (prevents silent OOM crashes).
#
# How to use in package code:
#   - Prefer .pstr_lapply(seq_len(B), FUN, parallel=parallel, workers=cpus)
#   - Keep 'X' small (indices), and keep large read-only objects as globals
#     (captured by the closure) so each worker receives them at most once.
# ------------------------------------------------------------------------------

# Check that optional parallel dependencies are installed.
# These are in Suggests, so they might not exist on the user's machine.
.pstr_need_future_apply <- function() {
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop(
      "Parallel option requires package 'future.apply'. Please install it.",
      call. = FALSE
    )
  }
  if (!requireNamespace("future", quietly = TRUE)) {
    stop(
      "Parallel option requires package 'future'. Please install it.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Normalise and validate number of workers.
# We cap at future::availableCores() to reduce the chance of accidental overload.
.pstr_normalise_workers <- function(workers) {
  if (is.null(workers)) return(NULL)
  
  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers)) {
    stop("'workers' must be a single non-missing numeric value.", call. = FALSE)
  }
  
  workers <- as.integer(workers)
  if (workers < 1L) stop("'workers' must be >= 1.", call. = FALSE)
  
  maxw <- future::availableCores()
  if (workers > maxw) {
    warning(
      sprintf(
        "'workers'=%d is larger than availableCores()=%d; using %d.",
        workers, maxw, maxw
      ),
      call. = FALSE
    )
    workers <- maxw
  }
  
  workers
}

# Internal helper: run lapply serially or in parallel.
#
# Arguments:
#   X        : list or vector of tasks (best practice: use small indices)
#   FUN      : function to apply to each element of X
#   ...      : forwarded to FUN
#   parallel : TRUE/FALSE, whether to use future.apply
#   workers  : integer, number of workers; NULL means "respect existing plan"
#   seed     : passed to future_lapply via future.seed for reproducible RNG
#
# Behaviour:
#   - If parallel=FALSE -> base::lapply
#   - If parallel=TRUE  -> future.apply::future_lapply with a safe plan
#     (multisession) unless the user has already set a non-sequential plan.
.pstr_lapply <- function(X, FUN, ..., parallel = FALSE, workers = NULL, seed = TRUE) {
  if (!parallel) {
    return(lapply(X, FUN, ...))
  }
  
  .pstr_need_future_apply()
  workers <- .pstr_normalise_workers(workers)
  
  # --------------------------------------------------------------------------
  # "Fuse" for accidentally exporting huge objects (globals) to workers.
  #
  # future inspects the closure and exports referenced objects to each worker.
  # If a user accidentally captures something very large (e.g., the whole data
  # tibble or a big R6 object), this can cause slowdowns or out-of-memory (OOM).
  #
  # Setting a maximum size makes this fail early with a clear error, instead of
  # crashing the session. You can adjust this value if you have a strong reason.
  # --------------------------------------------------------------------------
  old_max <- getOption("future.globals.maxSize")
  if (is.null(old_max) || old_max < 1024^3) {
    old_opt <- options(future.globals.maxSize = 1024^3)  # 1 GB
    on.exit(options(old_opt), add = TRUE)
  }
  
  # --------------------------------------------------------------------------
  # Respect the user's existing future plan where possible.
  #
  # - If 'workers' is provided, we choose multisession(workers=...).
  # - If 'workers' is NULL:
  #     * if the current plan is sequential, we set multisession()
  #     * otherwise we keep the user's plan (e.g., cluster, multisession, etc.)
  #
  # We always restore the original plan when leaving this function.
  # --------------------------------------------------------------------------
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  if (!is.null(workers)) {
    future::plan(future::multisession, workers = workers)
  } else {
    if (inherits(old_plan, "sequential")) {
      future::plan(future::multisession)
    }
  }
  
  future.apply::future_lapply(X, FUN, ..., future.seed = seed)
}

# Convenience wrapper similar to vapply(), with optional parallel evaluation.
#
# Note:
#   vapply enforces type/length consistency. In parallel mode we first compute
#   a list (like lapply) and then coerce/validate via vapply.
.pstr_vapply <- function(X, FUN, FUN.VALUE, ..., parallel = FALSE, workers = NULL, seed = TRUE) {
  if (!parallel) {
    return(vapply(X, FUN, FUN.VALUE = FUN.VALUE, ...))
  }
  
  out <- .pstr_lapply(X, FUN, ..., parallel = TRUE, workers = workers, seed = seed)
  
  # Coerce + validate: this will error if elements do not match FUN.VALUE,
  # which is exactly what vapply is meant to guarantee.
  vapply(out, function(z) z, FUN.VALUE = FUN.VALUE)
}