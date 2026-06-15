spar_time_checked <- function(x) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  !is.null(x$index$time_ordered)
}

spar_time_order <- function(x, use_identity_if_missing = TRUE) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  n <- nrow(x$data$X_original)

  if (!isTRUE(x$index$time_present)) {
    if (use_identity_if_missing) {
      return(seq_len(n))
    }
    stop("No time variable is available.", call. = FALSE)
  }

  if (!is.null(x$index$time_order)) {
    return(x$index$time_order)
  }

  order(x$data$time)
}

#' Validate a SPAR representation
#'
#' @param x spar_representation
#' @param warn_time_unordered logical value, if TRUE the validation function will return updated time information.
#' @export
validate_spar_representation <- function(x, warn_time_unordered = FALSE) {
  if (!inherits(x, "spar_representation")) {
    stop("Object is not a `spar_representation`.", call. = FALSE)
  }

  required_top <- c(
    "meta", "schema", "data", "transform", "angular",
    "threshold", "excess", "excursions", "fitted", "index"
  )

  missing_top <- setdiff(required_top, names(x))
  if (length(missing_top) > 0L) {
    stop(
      sprintf("Missing top-level fields: %s", paste(missing_top, collapse = ", ")),
      call. = FALSE
    )
  }

  X <- x$data$X_original
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`data$X_original` must be a numeric matrix.", call. = FALSE)
  }

  n <- nrow(X)
  d <- ncol(X)

  if (!identical(x$schema$original_dim, d)) {
    stop("`schema$original_dim` does not match `ncol(data$X_original)`.", call. = FALSE)
  }

  if (!identical(length(x$schema$original_names), d)) {
    stop("`schema$original_names` has invalid length.", call. = FALSE)
  }

  if (!identical(colnames(X), x$schema$original_names)) {
    stop("Column names of `data$X_original` do not match `schema$original_names`.", call. = FALSE)
  }

  if (!is.null(x$data$time) && length(x$data$time) != n) {
    stop("`data$time` has invalid length.", call. = FALSE)
  }

  if (length(x$data$observation_id) != n) {
    stop("`data$observation_id` has invalid length.", call. = FALSE)
  }

  if (anyDuplicated(x$data$observation_id)) {
    stop("`data$observation_id` must be unique.", call. = FALSE)
  }

  if (!is.null(x$angular$R) && length(x$angular$R) != n) {
    stop("`angular$R` has invalid length.", call. = FALSE)
  }

  if (!is.null(x$angular$phi) && length(x$angular$phi) != n) {
    stop("`angular$phi` has invalid length.", call. = FALSE)
  }

  if (!is.null(x$angular$active) && (!is.character(x$angular$active) || length(x$angular$active) != 1L || is.na(x$angular$active) || !nzchar(x$angular$active))) {
    stop("`angular$active` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  if (!is.null(x$angular$maps)) {
    if (!is.list(x$angular$maps)) {
      stop("`angular$maps` must be a list when present.", call. = FALSE)
    }
    if (length(x$angular$maps) > 0L) {
      nm <- names(x$angular$maps)
      if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
        stop("`angular$maps` must be a named list with non-empty ids.", call. = FALSE)
      }

      for (id in nm) {
        spec <- x$angular$maps[[id]]
        if (!is.list(spec)) {
          stop(sprintf("`angular$maps[['%s']]` must be a list.", id), call. = FALSE)
        }
        if (!is.function(spec$radial_fun) || !is.function(spec$angle_fun)) {
          stop(sprintf("`angular$maps[['%s']]` must contain valid radial/angle functions.", id), call. = FALSE)
        }
        if (!is.null(spec$inverse_fun) && !is.function(spec$inverse_fun)) {
          stop(sprintf("`angular$maps[['%s']]$inverse_fun` must be NULL or a function.", id), call. = FALSE)
        }
        if (!is.character(spec$source) || length(spec$source) != 1L || is.na(spec$source) || !(spec$source %in% c("transformed", "original"))) {
          stop(sprintf("`angular$maps[['%s']]$source` is invalid.", id), call. = FALSE)
        }
        if (!is.character(spec$format) || length(spec$format) != 1L || is.na(spec$format) || !(spec$format %in% c("matrix", "data.frame", "tibble"))) {
          stop(sprintf("`angular$maps[['%s']]$format` is invalid.", id), call. = FALSE)
        }
        as_spar_angle_domain(spec$domain)
      }
    }
  }

  if (!is.null(x$angular$active) && !is.null(x$angular$maps) && length(x$angular$maps) > 0L && !(x$angular$active %in% names(x$angular$maps))) {
    stop("`angular$active` must reference a key in `angular$maps`.", call. = FALSE)
  }

  if (!is.null(x$threshold$per_observation) &&
      length(x$threshold$per_observation) != n) {
    stop("`threshold$per_observation` has invalid length.", call. = FALSE)
  }

  if (!is.null(x$excess$value) && length(x$excess$value) != n) {
    stop("`excess$value` has invalid length.", call. = FALSE)
  }

  if (!is.null(x$excess$is_exceedance) && length(x$excess$is_exceedance) != n) {
    stop("`excess$is_exceedance` has invalid length.", call. = FALSE)
  }

  if (!is.function(x$transform$forward)) {
    stop("`transform$forward` must be a function.", call. = FALSE)
  }

  if (!is.function(x$transform$inverse)) {
    stop("`transform$inverse` must be a function.", call. = FALSE)
  }

  if (!is.null(x$transform$active) && (!is.character(x$transform$active) || length(x$transform$active) != 1L || is.na(x$transform$active) || !nzchar(x$transform$active))) {
    stop("`transform$active` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  if (!is.null(x$transform$chains)) {
    if (!is.list(x$transform$chains)) {
      stop("`transform$chains` must be a list when present.", call. = FALSE)
    }
    if (length(x$transform$chains) > 0L) {
      nm <- names(x$transform$chains)
      if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
        stop("`transform$chains` must be a named list with non-empty ids.", call. = FALSE)
      }
      for (id in nm) {
        validate_spar_transform_chain(x$transform$chains[[id]])
      }
    }
  }

  active_thr <- x$threshold$active
  if (!is.null(active_thr) && !is.character(active_thr) && !is.list(active_thr)) {
    stop("`threshold$active` must be NULL, a character scalar, or a named list map.", call. = FALSE)
  }

  reg <- x$threshold$registry
  if (!is.null(reg) && !is.data.frame(reg)) {
    stop("`threshold$registry` must be NULL or a data.frame.", call. = FALSE)
  }

  if(is.null(x$index$time_ordered)){
    # --- time index state ---
    time <- x$data$time

    if (is.null(time)) {
      x$index$time_present <- FALSE
      x$index$time_ordered <- NA
      x$index$time_order <- NULL
    } else {
      ord <- order(time)
      ordered_flag <- identical(ord, seq_along(time))

      x$index$time_present <- TRUE
      x$index$time_ordered <- ordered_flag
      x$index$time_order <- ord

      if (!ordered_flag && isTRUE(warn_time_unordered)) {
        warning(
          "Observations are not ordered by time. ",
          "Temporal procedures such as declustering may be slower or behave unexpectedly ",
          "unless an ordered index is used.",
          call. = FALSE
        )
      }
    }

  }

  return(x)
}
