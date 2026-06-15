#' Access observation-level data from a spar representation
#'
#' Returns the observation-level data stored in a `spar_representation`
#' object in one of several coordinate spaces.
#'
#' @param x A `spar_representation` object.
#' @param space Which data representation to return: `"original"`,
#'   `"transformed"`, or `"angular"`.
#' @param format Output format: `"matrix"`, `"data.frame"`, or `"tibble"`.
#' @param include_time Logical; if `TRUE`, include time when returning
#'   tabular output and time is available.
#' @param include_observation_id Logical; if `TRUE`, include observation IDs
#'   when returning tabular output.
#' @param order_time Logical; if `TRUE`, return observations ordered by time.
#' @param warn_time_unordered Logical; if `TRUE`, warn when time is present
#'   but observations are not ordered and `order_time = FALSE`.
#'
#' @return A matrix, data.frame, or tibble depending on `format`.
#' @export
spar_data <- function(
    x,
    space = c("original", "transformed", "angular"),
    format = c("matrix", "data.frame", "tibble"),
    include_time = FALSE,
    include_observation_id = FALSE,
    order_time = FALSE,
    warn_time_unordered = FALSE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  space <- match.arg(space)
  format <- match.arg(format)

  # --- extract core data for requested space ---
  core <- switch(
    space,
    original = {
      X <- x$data$X_original
      if (!is.matrix(X)) {
        X <- as.matrix(X)
      }
      X
    },

    transformed = {
      X <- x$transform$transformed_data
      if (is.null(X)) {
        stop(
          "No transformed data available in `x$transform$transformed_data`.",
          call. = FALSE
        )
      }
      if (!is.matrix(X)) {
        X <- as.matrix(X)
      }
      X
    },

    angular = {
      if (is.null(x$angular$R) || is.null(x$angular$phi)) {
        stop(
          "Angular representation has not been computed.",
          call. = FALSE
        )
      }
      X <- cbind(x$angular$R, x$angular$phi)
      X
    }
  )

  if (!is.matrix(core)) {
    core <- as.matrix(core)
  }

  core <- spar_apply_schema_names(x, core, space)

  n <- nrow(core)

  # --- time handling ---
  has_time <- isTRUE(x$index$time_present)
  time_checked <- !is.null(x$index$time_ordered)
  time_ordered <- isTRUE(x$index$time_ordered)

  if (order_time) {
    if (!has_time && format != "matrix") {
      warning(
        "`order_time = TRUE` requested, but no time variable is available. Returning original order.",
        call. = FALSE
      )
      ord <- seq_len(n)
    } else if(has_time) {
      if (is.null(x$index$time_order)) {
        ord <- order(x$data$time)
      } else {
        ord <- x$index$time_order
      }
    }
  } else {
    ord <- seq_len(n)

    if (has_time && warn_time_unordered) {
      if (!time_checked) {
        warning(
          "Time is present but ordering has not been checked. ",
          "Use validation or `order_time = TRUE` for temporal workflows.",
          call. = FALSE
        )
      } else if (!time_ordered) {
        warning(
          "Observations are not ordered by time. ",
          "Temporal procedures such as declustering may require `order_time = TRUE` ",
          "or use of the stored ordered index.",
          call. = FALSE
        )
      }
    }
  }

  core <- core[ord, , drop = FALSE]

  # --- matrix output: only core representation ---
  if (format == "matrix") {
    if (include_time || include_observation_id) {
      warning(
        "`include_time` and `include_observation_id` are ignored when `format = \"matrix\"`.",
        call. = FALSE
      )
    }
    return(core)
  }

  # --- tabular output with optional metadata ---
  out <- as.data.frame(core, stringsAsFactors = FALSE)

  meta <- spar_schema_meta_names(x)
  obs_name <- meta$observation_id
  if (include_observation_id) {
    out[[obs_name]] <- x$data$observation_id[ord]
    out <- out[, c(obs_name, setdiff(names(out), obs_name)), drop = FALSE]
  }

  if (include_time && has_time) {
    time_name <- meta$time
    out[[time_name]] <- x$data$time[ord]

    front <- c()
    if (include_observation_id) {
      front <- c(front, obs_name)
    }
    front <- c(front, time_name)

    other <- setdiff(names(out), front)
    out <- out[, c(front, other), drop = FALSE]
  }

  if (format == "tibble") {
    return(tibble::as_tibble(out))
  }

  out
}

#' Returns the field names for the different spaces of the SPAR model representation.
#'
#' @param x spar_representation
#' @param space Character, either `original`, `transformed` or `angular`
#' @returns The names used when accessed from `spar_data(...)`
#' @export
spar_schema_names <- function(x, space = c("original", "transformed", "angular")) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  space <- match.arg(space)
  sch <- x$schema

  switch(
    space,
    original = sch$original_names,
    transformed = sch$transformed_names %||% NULL,
    angular = c(sch$radial_name, sch$angle_name)
  )
}

#' Returns the time data stored in the SPAR model representation.
#'
#' @param x spar_representation
#' @param field Character, either `time`, `order` or `index`
#' @param format Character, either `vector`, `data.frame`, `tibble`
#' @param include_observation_id Logical, if `TRUE` the field `observation_id` will be added with its name as specified in the schema.
#' @param order_time Logical, if `TRUE` the data will be ordered by time if possible.
#' @returns The names used when accessed from `spar_data(...)`
#' @export
spar_time <- function(
    x,
    field = c("time", "order", "index"),
    format = c("vector", "data.frame", "tibble"),
    include_observation_id = FALSE,
    order_time = FALSE,
    warn_time_unordered = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  field <- match.arg(field)
  format <- match.arg(format)

  if (!isTRUE(x$index$time_present) || is.null(x$data$time)) {
    stop("No time variable is available in this `spar_representation`.", call. = FALSE)
  }

  time <- x$data$time
  n <- length(time)

  time_checked <- !is.null(x$index$time_ordered)
  time_ordered <- isTRUE(x$index$time_ordered)

  ord <- x$index$time_order
  if (is.null(ord)) {
    ord <- order(time)
  }

  if (!order_time && isTRUE(warn_time_unordered)) {
    if (!time_checked) {
      warning(
        "Time is present but ordering has not been checked. ",
        "Use validation or `order_time = TRUE` for temporal workflows.",
        call. = FALSE
      )
    } else if (!time_ordered) {
      warning(
        "Observations are not ordered by time. ",
        "Use `order_time = TRUE` or request `field = \"order\"` for temporal workflows.",
        call. = FALSE
      )
    }
  }

  if (field == "order") {
    return(ord)
  }

  use_idx <- if (order_time) ord else seq_len(n)
  meta_names <- spar_schema_meta_names(x)
  obs_name <- meta_names$observation_id
  time_name <- meta_names$time

  if (field == "time") {
    out_time <- time[use_idx]

    if (format == "vector") {
      return(out_time)
    }

    out <- data.frame(stringsAsFactors = FALSE)

    if (include_observation_id) {
      out[[obs_name]] <- x$data$observation_id[use_idx]
    }

    out[[time_name]] <- out_time

    if (format == "tibble") {
      return(tibble::as_tibble(out))
    }
    return(out)
  }

  # what == "index"
  out <- data.frame(stringsAsFactors = FALSE)

  if (include_observation_id) {
    out[[obs_name]] <- x$data$observation_id[use_idx]
  }

  out$row_index <- use_idx
  out[[time_name]] <- time[use_idx]

  if (format == "tibble") {
    return(tibble::as_tibble(out))
  }

  out
}

spar_component <- function(x, name) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  name <- match.arg(
    name,
    c("data", "time", "observation_id", "angular", "threshold", "excess")
  )

  switch(
    name,

    data = {
      X <- x$data$X_original
      spar_apply_schema_names(x, X, "original")
    },

    time = x$data$time,

    observation_id = x$data$observation_id,

    angular = {
      if (is.null(x$angular$R) || is.null(x$angular$phi)) {
        return(NULL)
      }
      A <- cbind(x$angular$R, x$angular$phi)
      spar_apply_schema_names(x, A, "angular")
    },

    threshold = x$threshold$per_observation,

    excess = {
      if (is.null(x$excess$value) && is.null(x$excess$is_exceedance)) {
        return(NULL)
      }
      list(
        value = x$excess$value,
        is_exceedance = x$excess$is_exceedance
      )
    }
  )
}

spar_subset <- function(
    x,
    i,
    components = c("data", "time", "observation_id", "angular", "threshold", "excess")
) {

  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (missing(i) || is.null(i)) {
    i <- seq_len(nrow(x$data$X_original))
  }

  components <- unique(components)

  out <- lapply(components, function(comp) {
    spar_slice_obs(
      spar_component(x, comp),
      i
    )
  })

  names(out) <- components

  out
}
