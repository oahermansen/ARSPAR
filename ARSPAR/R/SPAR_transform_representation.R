#' Attach a compiled transform chain to a SPAR representation
#'
#' Stores a compiled `spar_transform_chain` on a `spar_representation` and,
#' optionally, executes it immediately to populate transformed data.
#'
#' @param x A `spar_representation` object.
#' @param chain A `spar_transform_chain` object.
#' @param transform_id Optional transform identifier. Defaults to
#'   `chain$name`, falling back to an auto-generated id.
#' @param set_active Logical; if `TRUE`, make this transform active.
#' @param run Logical; if `TRUE`, execute the chain on `x$data$X_original`.
#' @param keep_step_cache Logical; if `TRUE`, retain runtime step cache produced
#'   when running the chain.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_set_transform_chain <- function(
    x,
    chain,
    transform_id = NULL,
    set_active = TRUE,
    run = TRUE,
    keep_step_cache = FALSE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  validate_spar_transform_chain(chain)

  chain <- spar_finalize_transform_chain(chain)

  x <- spar_transform_registry_normalize(x)

  if (is.null(transform_id)) {
    transform_id <- chain$name %||% NULL
  }

  if (is.null(transform_id) || !is.character(transform_id) || length(transform_id) != 1L || is.na(transform_id) || !nzchar(transform_id)) {
    transform_id <- sprintf("transform_%d", length(x$transform$chains) + 1L)
  }

  chain$name <- transform_id
  x$transform$chains[[transform_id]] <- chain

  if (!isTRUE(set_active)) {
    if (isTRUE(run)) {
      warning("`run = TRUE` ignored because `set_active = FALSE`.", call. = FALSE)
    }
    return(x)
  }

  x$transform$active <- transform_id

  x$transform$steps <- chain$steps
  x$transform$chain <- chain
  x$transform$cache_requirements <- chain$cache_requirements
  x$transform$analysis <- chain$analysis

  x$transform$forward <- function(X) {
    spar_run_transform_chain(chain, X)$current
  }

  x$transform$inverse <- spar_chain_inverse_or_stub(chain)

  if (isTRUE(run)) {
    x <- spar_apply_transform_chain(x, keep_step_cache = keep_step_cache)
  }

  x
}

#' Build inverse function or a descriptive stub
#'
#' @param chain A `spar_transform_chain` object.
#'
#' @return A function accepting transformed data.
#' @keywords internal
spar_chain_inverse_or_stub <- function(chain) {
  inv <- chain$inverse %||% NULL
  if (is.function(inv)) {
    return(inv)
  }

  id <- chain$name %||% "active"
  function(Xt) {
    stop(sprintf("No inverse transform is available for transform '%s'.", id), call. = FALSE)
  }
}

#' Normalize transform registry structure
#'
#' Ensures transform registry fields exist and migrates legacy single-chain
#' layouts into a chain registry.
#'
#' @param x A `spar_representation` object.
#'
#' @return Updated `spar_representation`.
#' @keywords internal
spar_transform_registry_normalize <- function(x) {
  chains <- x$transform$chains
  if (is.null(chains)) {
    chains <- list()
  }
  if (!is.list(chains)) {
    stop("`transform$chains` must be a list.", call. = FALSE)
  }

  if (length(chains) > 0L) {
    nm <- names(chains)
    if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
      stop("`transform$chains` must be a named list with non-empty ids.", call. = FALSE)
    }
    for (id in nm) {
      validate_spar_transform_chain(chains[[id]])
      chains[[id]]$name <- id
    }
  }

  if (is.null(x$transform$chain) && length(x$transform$steps) > 0L) {
    legacy_id <- x$transform$active %||% "active"
    if (!is.character(legacy_id) || length(legacy_id) != 1L || is.na(legacy_id) || !nzchar(legacy_id)) {
      legacy_id <- "active"
    }
    ch <- new_spar_transform_chain(steps = x$transform$steps, name = legacy_id)
    ch <- spar_finalize_transform_chain(ch)
    chains[[legacy_id]] <- ch
    x$transform$chain <- ch
  }

  active <- x$transform$active %||% NULL
  if (!is.null(active) && (!is.character(active) || length(active) != 1L || is.na(active) || !nzchar(active))) {
    active <- NULL
  }

  if (is.null(active)) {
    if (!is.null(x$transform$chain) && is.list(x$transform$chain)) {
      active <- x$transform$chain$name %||% NULL
    }
    if (is.null(active) && length(chains) == 1L) {
      active <- names(chains)[1L]
    }
    if (is.null(active) && length(x$transform$steps) == 0L) {
      active <- "identity"
    }
  }

  if (!is.null(active) && !(active %in% names(chains)) && length(x$transform$steps) > 0L) {
    ch <- new_spar_transform_chain(steps = x$transform$steps, name = active)
    ch <- spar_finalize_transform_chain(ch)
    chains[[active]] <- ch
  }

  x$transform$chains <- chains
  x$transform$active <- active
  x
}

#' Get active transform identifier
#'
#' @param x A `spar_representation` object.
#'
#' @return Single character transform identifier.
#' @export
spar_active_transform_id <- function(x) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  x <- spar_transform_registry_normalize(x)
  id <- x$transform$active %||% "identity"
  as.character(id)
}

#' Set the active transform by identifier
#'
#' Activates a transform chain stored in `x$transform$chains` and optionally
#' computes its transformed data into `x$transform$transformed_data`.
#'
#' @param x A `spar_representation` object.
#' @param transform_id Transform identifier.
#' @param run Logical; if `TRUE`, recompute transformed data for active chain.
#' @param keep_step_cache Logical; if `TRUE`, retain runtime step cache.
#'
#' @return Updated `spar_representation`.
#' @export
spar_set_active_transform <- function(
    x,
    transform_id,
    run = TRUE,
    keep_step_cache = FALSE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.character(transform_id) || length(transform_id) != 1L || is.na(transform_id) || !nzchar(transform_id)) {
    stop("`transform_id` must be a single non-empty character string.", call. = FALSE)
  }

  x <- spar_transform_registry_normalize(x)

  if (identical(transform_id, "identity")) {
    x$transform$active <- "identity"
    x$transform$steps <- list()
    x$transform$chain <- NULL
    x$transform$forward <- function(X) X
    x$transform$inverse <- function(X) X
    x$transform$cache_requirements <- list(forward = list(), reverse = NULL, combined = list())
    x$transform$analysis <- NULL
    if (isTRUE(run)) {
      Xt <- x$data$X_original
      if (!is.matrix(Xt)) Xt <- as.matrix(Xt)
      x$transform$transformed_data <- Xt
      x$schema$transformed_names <- colnames(Xt)
      x$transform$step_cache <- if (isTRUE(keep_step_cache)) list() else NULL
    }
    return(x)
  }

  chain <- x$transform$chains[[transform_id]]
  if (is.null(chain)) {
    stop(sprintf("Unknown transform id '%s'.", transform_id), call. = FALSE)
  }
  validate_spar_transform_chain(chain)

  x$transform$active <- transform_id
  x$transform$steps <- chain$steps
  x$transform$chain <- chain
  x$transform$cache_requirements <- chain$cache_requirements
  x$transform$analysis <- chain$analysis
  x$transform$forward <- function(X) spar_run_transform_chain(chain, X)$current
  x$transform$inverse <- spar_chain_inverse_or_stub(chain)

  if (isTRUE(run)) {
    x <- spar_apply_transform_chain(x, keep_step_cache = keep_step_cache)
  }

  x
}

#' Get transformed data for a specific transform id
#'
#' Computes transformed data for `transform_id` without changing
#' `x$transform$active`.
#'
#' @param x A `spar_representation` object.
#' @param transform_id Transform identifier.
#' @param format Output format: `"matrix"`, `"data.frame"`, or `"tibble"`.
#'
#' @return Transformed data in requested format.
#' @export
spar_get_transformed_data <- function(
    x,
    transform_id,
    format = c("matrix", "data.frame", "tibble")
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.character(transform_id) || length(transform_id) != 1L || is.na(transform_id) || !nzchar(transform_id)) {
    stop("`transform_id` must be a single non-empty character string.", call. = FALSE)
  }

  format <- match.arg(format)
  x <- spar_transform_registry_normalize(x)

  Xt <- NULL
  if (identical(transform_id, x$transform$active) && !is.null(x$transform$transformed_data)) {
    Xt <- x$transform$transformed_data
  } else if (identical(transform_id, "identity")) {
    Xt <- x$data$X_original
  } else {
    chain <- x$transform$chains[[transform_id]]
    if (is.null(chain)) {
      stop(sprintf("Unknown transform id '%s'.", transform_id), call. = FALSE)
    }
    Xt <- spar_run_transform_chain(chain, x$data$X_original)$current
  }

  if (!is.matrix(Xt)) Xt <- as.matrix(Xt)

  if (identical(format, "matrix")) {
    return(Xt)
  }
  out <- as.data.frame(Xt, stringsAsFactors = FALSE)
  if (identical(format, "tibble")) {
    return(tibble::as_tibble(out))
  }
  out
}

#' Set inverse function for a transform id
#'
#' Registers a user-defined inverse function on a stored transform chain.
#'
#' @param x A `spar_representation` object.
#' @param inverse Inverse function.
#' @param transform_id Optional transform identifier. Defaults to active.
#'
#' @return Updated `spar_representation`.
#' @export
spar_set_transform_inverse <- function(x, inverse, transform_id = NULL) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.function(inverse)) {
    stop("`inverse` must be a function.", call. = FALSE)
  }

  x <- spar_transform_registry_normalize(x)
  if (is.null(transform_id)) {
    transform_id <- x$transform$active %||% "identity"
  }
  if (!is.character(transform_id) || length(transform_id) != 1L || is.na(transform_id) || !nzchar(transform_id)) {
    stop("`transform_id` must be a single non-empty character string.", call. = FALSE)
  }

  if (identical(transform_id, "identity")) {
    x$transform$inverse <- inverse
    return(x)
  }

  chain <- x$transform$chains[[transform_id]]
  if (is.null(chain)) {
    stop(sprintf("Unknown transform id '%s'.", transform_id), call. = FALSE)
  }

  chain$inverse <- inverse
  x$transform$chains[[transform_id]] <- chain

  if (identical(x$transform$active, transform_id)) {
    x$transform$inverse <- inverse
    x$transform$chain <- chain
  }

  x
}

#' Execute the active transform chain on a SPAR representation
#'
#' Runs the transform chain currently attached to `x` and stores the resulting
#' transformed matrix in `x$transform$transformed_data`.
#'
#' @param x A `spar_representation` object.
#' @param keep_step_cache Logical; if `TRUE`, retain runtime step cache.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_apply_transform_chain <- function(x, keep_step_cache = FALSE) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  x <- spar_transform_registry_normalize(x)

  active_id <- x$transform$active %||% "identity"
  if (identical(active_id, "identity")) {
    Xt <- x$data$X_original
    if (!is.matrix(Xt)) {
      Xt <- as.matrix(Xt)
    }

    x$transform$steps <- list()
    x$transform$chain <- NULL
    x$transform$forward <- function(X) X
    x$transform$inverse <- function(X) X
    x$transform$transformed_data <- Xt
    x$schema$transformed_names <- colnames(Xt)
    x$transform$cache_requirements <- list(forward = list(), reverse = NULL, combined = list())
    x$transform$analysis <- NULL
    x$transform$step_cache <- if (isTRUE(keep_step_cache)) list() else NULL
    return(x)
  }

  chain <- x$transform$chains[[active_id]]
  if (!is.null(chain)) {
    x$transform$chain <- chain
  }
  chain <- x$transform$chain
  if (is.null(chain)) {
    if (length(x$transform$steps) == 0L) {
      Xt <- x$data$X_original
      if (!is.matrix(Xt)) {
        Xt <- as.matrix(Xt)
      }

      x$transform$transformed_data <- Xt
      x$schema$transformed_names <- colnames(Xt)
      x$transform$cache_requirements <- list(
        forward = list(),
        reverse = NULL,
        combined = list()
      )
      x$transform$step_cache <- if (isTRUE(keep_step_cache)) list() else NULL
      return(x)
    }

    chain <- new_spar_transform_chain(
      steps = x$transform$steps,
      name = active_id
    )
    chain <- spar_finalize_transform_chain(chain)
    x$transform$chain <- chain
    x$transform$chains[[active_id]] <- chain
  } else {
    validate_spar_transform_chain(chain)
  }

  run <- spar_run_transform_chain(
    chain = chain,
    X_original = x$data$X_original
  )

  Xt <- run$current
  if (!is.matrix(Xt)) {
    Xt <- as.matrix(Xt)
  }

  x$transform$transformed_data <- Xt
  x$schema$transformed_names <- colnames(Xt)
  x$transform$steps <- chain$steps
  x$transform$cache_requirements <- chain$cache_requirements
  x$transform$analysis <- chain$analysis
  x$transform$forward <- function(X) spar_run_transform_chain(chain, X)$current
  x$transform$inverse <- spar_chain_inverse_or_stub(chain)
  x$transform$cache_requirements <- list(
    forward = run$cache_requirements,
    reverse = NULL,
    combined = run$cache_requirements
  )

  if (isTRUE(keep_step_cache)) {
    x$transform$step_cache <- run$step_cache
  } else {
    x$transform$step_cache <- NULL
  }

  x
}

#' Build, attach, and optionally run a transform chain
#'
#' Convenience wrapper that builds a chain from step specifications using the
#' original data in `x`, attaches it to `x`, and optionally runs it.
#'
#' @param x A `spar_representation` object.
#' @param ... Step specifications, typically from [spar_step_mutate()].
#' @param name Optional transform-chain name.
#' @param run Logical; if `TRUE`, execute immediately.
#' @param keep_step_cache Logical; if `TRUE`, retain runtime step cache.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_build_representation_transform <- function(
    x,
    ...,
    name = NULL,
    run = TRUE,
    keep_step_cache = FALSE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  chain <- spar_build_transform_chain(
    X_original = x$data$X_original,
    ...,
    name = name
  )

  spar_set_transform_chain(
    x = x,
    chain = chain,
    run = run,
    keep_step_cache = keep_step_cache
  )
}
