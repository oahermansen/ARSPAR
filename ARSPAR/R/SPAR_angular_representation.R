#' Build an angular-map specification
#'
#' @param radial_fun Function mapping observations to radial values.
#' @param angle_fun Function mapping observations to angular values.
#' @param domain Angular domain object.
#' @param source Default source space, one of `"transformed"` or `"original"`.
#' @param format Input format expected by mapping functions.
#' @param inverse_fun Optional inverse map function.
#'
#' @return Angular-map specification list.
#' @keywords internal
spar_new_angular_map_spec <- function(
    radial_fun,
    angle_fun,
    domain,
    source,
    format,
    inverse_fun = NULL
) {
  list(
    radial_fun = radial_fun,
    angle_fun = angle_fun,
    domain = as_spar_angle_domain(domain),
    source = source,
    format = format,
    inverse_fun = inverse_fun
  )
}

#' Validate angular-map specification
#'
#' @param spec Angular-map specification list.
#' @param id Optional map identifier used in error messages.
#'
#' @return `TRUE` invisibly.
#' @keywords internal
spar_validate_angular_map_spec <- function(spec, id = NULL) {
  tag <- if (is.null(id)) "angular map" else sprintf("angular map '%s'", id)
  if (!is.list(spec)) {
    stop(sprintf("%s must be a list.", tag), call. = FALSE)
  }

  if (!is.function(spec$radial_fun)) {
    stop(sprintf("%s requires a valid `radial_fun`.", tag), call. = FALSE)
  }
  if (!is.function(spec$angle_fun)) {
    stop(sprintf("%s requires a valid `angle_fun`.", tag), call. = FALSE)
  }
  if (!is.null(spec$inverse_fun) && !is.function(spec$inverse_fun)) {
    stop(sprintf("%s has invalid `inverse_fun`.", tag), call. = FALSE)
  }

  if (!is.character(spec$source) || length(spec$source) != 1L || is.na(spec$source) || !(spec$source %in% c("transformed", "original"))) {
    stop(sprintf("%s has invalid `source`.", tag), call. = FALSE)
  }

  if (!is.character(spec$format) || length(spec$format) != 1L || is.na(spec$format) || !(spec$format %in% c("matrix", "data.frame", "tibble"))) {
    stop(sprintf("%s has invalid `format`.", tag), call. = FALSE)
  }

  as_spar_angle_domain(spec$domain)
  invisible(TRUE)
}

#' Normalize angular registry structure
#'
#' Ensures angular map registry fields exist and migrates legacy single-map
#' layouts (`gauge`/`angle_map`) into the map registry.
#'
#' @param x A `spar_representation` object.
#'
#' @return Updated `spar_representation`.
#' @keywords internal
spar_angular_registry_normalize <- function(x) {
  maps <- x$angular$maps
  if (is.null(maps)) {
    maps <- list()
  }
  if (!is.list(maps)) {
    stop("`angular$maps` must be a list.", call. = FALSE)
  }

  if (length(maps) > 0L) {
    nm <- names(maps)
    if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
      stop("`angular$maps` must be a named list with non-empty ids.", call. = FALSE)
    }
    for (id in nm) {
      spar_validate_angular_map_spec(maps[[id]], id = id)
    }
  }

  if (length(maps) == 0L && is.function(x$angular$gauge) && is.function(x$angular$angle_map)) {
    legacy_id <- x$angular$active %||% "active"
    if (!is.character(legacy_id) || length(legacy_id) != 1L || is.na(legacy_id) || !nzchar(legacy_id)) {
      legacy_id <- "active"
    }

    source <- x$angular$source %||% "transformed"
    if (!(source %in% c("transformed", "original"))) {
      source <- "transformed"
    }

    format <- x$angular$format %||% "matrix"
    if (!(format %in% c("matrix", "data.frame", "tibble"))) {
      format <- "matrix"
    }

    domain <- x$angular$domain %||% new_spar_angle_domain()
    maps[[legacy_id]] <- spar_new_angular_map_spec(
      radial_fun = x$angular$gauge,
      angle_fun = x$angular$angle_map,
      domain = domain,
      source = source,
      format = format,
      inverse_fun = x$angular$inverse %||% NULL
    )
  }

  active <- x$angular$active %||% NULL
  if (!is.null(active) && (!is.character(active) || length(active) != 1L || is.na(active) || !nzchar(active))) {
    active <- NULL
  }

  if (is.null(active) && length(maps) == 1L) {
    active <- names(maps)[1L]
  }

  if (!is.null(active) && !(active %in% names(maps))) {
    active <- NULL
  }

  x$angular$maps <- maps
  x$angular$active <- active
  x
}

#' Synchronize active angular convenience fields
#'
#' @param x A `spar_representation` object.
#'
#' @return Updated `spar_representation`.
#' @keywords internal
spar_sync_active_angular_fields <- function(x) {
  x <- spar_angular_registry_normalize(x)
  id <- x$angular$active
  if (is.null(id)) {
    x$angular$gauge <- NULL
    x$angular$angle_map <- NULL
    x$angular$domain <- NULL
    x$angular$source <- NULL
    x$angular$format <- NULL
    x$angular$inverse <- NULL
    return(x)
  }

  spec <- x$angular$maps[[id]]
  x$angular$gauge <- spec$radial_fun
  x$angular$angle_map <- spec$angle_fun
  x$angular$domain <- spec$domain
  x$angular$source <- spec$source
  x$angular$format <- spec$format
  x$angular$inverse <- spec$inverse_fun
  x
}

#' Compute angular map data without mutating active state
#'
#' @param x A `spar_representation` object.
#' @param spec Angular-map specification.
#' @param source Optional source override.
#' @param transform_id Optional transform id override for transformed source.
#'
#' @return List with `X`, `R`, `phi`, `source`, and `transform_id`.
#' @keywords internal
spar_compute_angular_map_data <- function(x, spec, source = NULL, transform_id = NULL) {
  spar_validate_angular_map_spec(spec)

  source_use <- source %||% spec$source
  source_use <- match.arg(source_use, c("transformed", "original"))

  transform_use <- NULL
  if (identical(source_use, "original")) {
    X <- spar_data(x, space = "original", format = spec$format)
  } else {
    transform_use <- transform_id
    if (is.null(transform_use)) {
      if (exists("spar_active_transform_id", mode = "function")) {
        transform_use <- spar_active_transform_id(x)
      } else {
        transform_use <- x$transform$active %||% "identity"
      }
    }
    X <- spar_get_transformed_data(x, transform_id = transform_use, format = spec$format)
  }

  n <- nrow(X)
  R <- spec$radial_fun(X)
  phi <- spec$angle_fun(X)

  if (is.matrix(R) || is.data.frame(R)) {
    if (ncol(R) != 1L) {
      stop("`radial_fun` must return a vector or one-column object.", call. = FALSE)
    }
    R <- R[, 1L, drop = TRUE]
  }
  if (is.matrix(phi) || is.data.frame(phi)) {
    if (ncol(phi) != 1L) {
      stop("`angle_fun` must return a vector or one-column object.", call. = FALSE)
    }
    phi <- phi[, 1L, drop = TRUE]
  }

  R <- as.numeric(R)
  phi <- as.numeric(phi)
  if (length(R) != n) {
    stop(sprintf("`radial_fun` returned length %d, expected %d.", length(R), n), call. = FALSE)
  }
  if (length(phi) != n) {
    stop(sprintf("`angle_fun` returned length %d, expected %d.", length(phi), n), call. = FALSE)
  }

  dom <- as_spar_angle_domain(spec$domain)
  phi <- spar_normalize_angle(phi, dom)

  list(
    X = X,
    R = R,
    phi = phi,
    source = source_use,
    transform_id = transform_use,
    domain = dom
  )
}

#' Get active angular map identifier
#'
#' @param x A `spar_representation` object.
#'
#' @return Active angular map identifier or `NULL` if none is active.
#' @export
spar_active_angular_id <- function(x) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  x <- spar_angular_registry_normalize(x)
  x$angular$active
}

#' Attach angular mapping functions to a SPAR representation
#'
#' Stores radial and angle mapping functions on a `spar_representation` and,
#' optionally, executes them immediately.
#'
#' @param x A `spar_representation` object.
#' @param radial_fun Function mapping an observation matrix to a numeric radial
#'   vector of length `nrow(X)`.
#' @param angle_fun Function mapping an observation matrix to a numeric angular
#'   vector of length `nrow(X)`.
#' @param domain Angular domain object, either `spar_angle_domain` or
#'   legacy `angle_domain`.
#' @param name Optional map identifier (legacy alias for `angular_id`).
#' @param angular_id Optional map identifier.
#' @param source Default source space used when applying mappings, one of
#'   `"transformed"` or `"original"`.
#' @param format Function parameter expected format, one of `"matrix"`,
#'   `"data.frame"` or `"tibble"`. `"matrix"` is default.
#' @param inverse_fun Optional inverse map function.
#' @param set_active Logical; if `TRUE`, make this map active.
#' @param run Logical; if `TRUE`, run immediately when active.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_set_angular_map <- function(
    x,
    radial_fun,
    angle_fun,
    domain = new_spar_angle_domain(),
    name = NULL,
    angular_id = NULL,
    source = c("transformed", "original"),
    format = c("matrix", "data.frame", "tibble"),
    inverse_fun = NULL,
    set_active = TRUE,
    run = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.function(radial_fun)) {
    stop("`radial_fun` must be a function.", call. = FALSE)
  }
  if (!is.function(angle_fun)) {
    stop("`angle_fun` must be a function.", call. = FALSE)
  }
  if (!is.null(inverse_fun) && !is.function(inverse_fun)) {
    stop("`inverse_fun` must be NULL or a function.", call. = FALSE)
  }

  source <- match.arg(source)
  format <- match.arg(format)
  domain <- as_spar_angle_domain(domain)

  x <- spar_angular_registry_normalize(x)

  id <- angular_id %||% name %||% NULL
  if (is.null(id)) {
    id <- sprintf("angular_%d", length(x$angular$maps) + 1L)
  }
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("`angular_id`/`name` must be a single non-empty character string.", call. = FALSE)
  }

  spec <- spar_new_angular_map_spec(
    radial_fun = radial_fun,
    angle_fun = angle_fun,
    domain = domain,
    source = source,
    format = format,
    inverse_fun = inverse_fun
  )
  spar_validate_angular_map_spec(spec, id = id)

  x$angular$maps[[id]] <- spec

  if (!isTRUE(set_active)) {
    if (isTRUE(run)) {
      warning("`run = TRUE` ignored because `set_active = FALSE`.", call. = FALSE)
    }
    return(x)
  }

  spar_set_active_angular_map(x, angular_id = id, run = run)
}

#' Set active angular map
#'
#' Activates a stored angular map and optionally computes active `R`/`phi`.
#'
#' @param x A `spar_representation` object.
#' @param angular_id Angular map identifier.
#' @param run Logical; if `TRUE`, compute active angular representation.
#' @param source Optional source override for immediate run.
#' @param transform_id Optional transform id override used when `source`
#'   resolves to `"transformed"`.
#'
#' @return Updated `spar_representation`.
#' @export
spar_set_active_angular_map <- function(
    x,
    angular_id,
    run = TRUE,
    source = NULL,
    transform_id = NULL
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.character(angular_id) || length(angular_id) != 1L || is.na(angular_id) || !nzchar(angular_id)) {
    stop("`angular_id` must be a single non-empty character string.", call. = FALSE)
  }

  x <- spar_angular_registry_normalize(x)
  spec <- x$angular$maps[[angular_id]]
  if (is.null(spec)) {
    stop(sprintf("Unknown angular map id '%s'.", angular_id), call. = FALSE)
  }

  x$angular$active <- angular_id
  x <- spar_sync_active_angular_fields(x)

  if (isTRUE(run)) {
    x <- spar_apply_angular_map(x, source = source, transform_id = transform_id)
  }

  x
}

#' Set inverse function for an angular map
#'
#' @param x A `spar_representation` object.
#' @param inverse_fun Inverse map function.
#' @param angular_id Optional angular map identifier. Defaults to active map.
#'
#' @return Updated `spar_representation`.
#' @export
spar_set_angular_inverse <- function(x, inverse_fun, angular_id = NULL) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.function(inverse_fun)) {
    stop("`inverse_fun` must be a function.", call. = FALSE)
  }

  x <- spar_angular_registry_normalize(x)
  if (is.null(angular_id)) {
    angular_id <- x$angular$active
  }
  if (!is.character(angular_id) || length(angular_id) != 1L || is.na(angular_id) || !nzchar(angular_id)) {
    stop("`angular_id` must be a single non-empty character string.", call. = FALSE)
  }

  spec <- x$angular$maps[[angular_id]]
  if (is.null(spec)) {
    stop(sprintf("Unknown angular map id '%s'.", angular_id), call. = FALSE)
  }

  spec$inverse_fun <- inverse_fun
  x$angular$maps[[angular_id]] <- spec

  if (identical(x$angular$active, angular_id)) {
    x <- spar_sync_active_angular_fields(x)
  }

  x
}

#' Get angular data for a specific angular map without changing active state
#'
#' Computes `R` and `phi` for a specified angular map id without mutating
#' `x$angular$active`.
#'
#' @param x A `spar_representation` object.
#' @param angular_id Angular map identifier.
#' @param transform_id Optional transform id used when map source is
#'   `"transformed"`.
#' @param source Optional source override.
#' @param format Output format: `"matrix"`, `"data.frame"`, or `"tibble"`.
#'
#' @return Angular data in requested format with columns `R` and `phi`.
#' @export
spar_get_angular_data <- function(
    x,
    angular_id,
    transform_id = NULL,
    source = NULL,
    format = c("matrix", "data.frame", "tibble")
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.character(angular_id) || length(angular_id) != 1L || is.na(angular_id) || !nzchar(angular_id)) {
    stop("`angular_id` must be a single non-empty character string.", call. = FALSE)
  }

  format <- match.arg(format)
  x <- spar_angular_registry_normalize(x)
  spec <- x$angular$maps[[angular_id]]
  if (is.null(spec)) {
    stop(sprintf("Unknown angular map id '%s'.", angular_id), call. = FALSE)
  }

  out <- spar_compute_angular_map_data(x, spec = spec, source = source, transform_id = transform_id)
  mat <- cbind(R = out$R, phi = out$phi)

  if (identical(format, "matrix")) {
    return(mat)
  }

  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  if (identical(format, "tibble")) {
    return(tibble::as_tibble(df))
  }
  df
}

#' Apply the active angular mapping to a SPAR representation
#'
#' Executes radial and angle mappings and stores angular representation fields on
#' the object.
#'
#' @param x A `spar_representation` object.
#' @param source Source data space, one of `"transformed"` or `"original"`.
#'   If `NULL`, uses the stored default source.
#' @param transform_id Optional transform id used when `source` resolves to
#'   `"transformed"`.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_apply_angular_map <- function(x, source = NULL, transform_id = NULL) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  x <- spar_angular_registry_normalize(x)
  x <- spar_sync_active_angular_fields(x)

  id <- x$angular$active
  if (is.null(id)) {
    stop("No active angular map is set.", call. = FALSE)
  }

  spec <- x$angular$maps[[id]]
  out <- spar_compute_angular_map_data(x, spec = spec, source = source, transform_id = transform_id)

  x$angular$X <- out$X
  x$angular$R <- out$R
  x$angular$phi <- out$phi
  x$angular$domain <- out$domain
  x$angular$source <- out$source
  x$angular$transform_id <- out$transform_id

  x
}

#' Apply inverse angular map
#'
#' Applies the user-defined inverse map for an angular map id.
#'
#' @param x A `spar_representation` object.
#' @param R Optional radial values. Defaults to active/computed values.
#' @param phi Optional angular values. Defaults to active/computed values.
#' @param angular_id Optional angular map identifier. Defaults to active map.
#'
#' @return Reconstructed data as a matrix.
#' @export
spar_apply_angular_inverse <- function(x, R = NULL, phi = NULL, angular_id = NULL) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  x <- spar_angular_registry_normalize(x)
  if (is.null(angular_id)) {
    angular_id <- x$angular$active
  }
  if (!is.character(angular_id) || length(angular_id) != 1L || is.na(angular_id) || !nzchar(angular_id)) {
    stop("`angular_id` must be a single non-empty character string.", call. = FALSE)
  }

  spec <- x$angular$maps[[angular_id]]
  if (is.null(spec)) {
    stop(sprintf("Unknown angular map id '%s'.", angular_id), call. = FALSE)
  }
  inv <- spec$inverse_fun
  if (!is.function(inv)) {
    stop(sprintf("No inverse map is registered for angular map '%s'.", angular_id), call. = FALSE)
  }

  if (is.null(R) || is.null(phi)) {
    if (identical(angular_id, x$angular$active) && !is.null(x$angular$R) && !is.null(x$angular$phi)) {
      R <- x$angular$R
      phi <- x$angular$phi
    } else {
      ang <- spar_get_angular_data(x, angular_id = angular_id, format = "matrix")
      R <- ang[, 1L]
      phi <- ang[, 2L]
    }
  }

  R <- as.numeric(R)
  phi <- as.numeric(phi)
  if (length(R) != length(phi)) {
    stop("`R` and `phi` must have equal lengths.", call. = FALSE)
  }

  fml <- names(formals(inv))
  if (!is.null(fml) && ("R" %in% fml || "phi" %in% fml || "x" %in% fml || "angular_id" %in% fml)) {
    args <- list()
    if ("R" %in% fml) args$R <- R
    if ("phi" %in% fml) args$phi <- phi
    if ("x" %in% fml) args$x <- x
    if ("angular_id" %in% fml) args$angular_id <- angular_id
    out <- do.call(inv, args)
  } else {
    out <- inv(cbind(R = R, phi = phi))
  }

  out <- as.matrix(out)
  if (!is.numeric(out) || nrow(out) != length(R)) {
    stop("Inverse angular map must return a numeric matrix with one row per observation.", call. = FALSE)
  }

  out
}

#' Build and apply an angular representation on a SPAR object
#'
#' Convenience wrapper that attaches an angular map and applies it.
#'
#' @param x A `spar_representation` object.
#' @param radial_fun Function returning radial values from matrix input.
#' @param angle_fun Function returning angular values from matrix input.
#' @param domain Angular domain object.
#' @param source Source data space, one of `"transformed"` or `"original"`.
#' @param inverse_fun Optional inverse map function.
#' @param name Optional active mapping label.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_build_angular_representation <- function(
    x,
    radial_fun,
    angle_fun,
    domain = new_spar_angle_domain(),
    source = c("transformed", "original"),
    inverse_fun = NULL,
    name = NULL
) {
  source <- match.arg(source)

  spar_set_angular_map(
    x = x,
    radial_fun = radial_fun,
    angle_fun = angle_fun,
    domain = domain,
    name = name,
    source = source,
    inverse_fun = inverse_fun,
    run = TRUE
  )
}
