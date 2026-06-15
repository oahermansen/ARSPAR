#' Resolve a threshold vector for a SPAR representation
#'
#' @param x A `spar_representation` object.
#' @param u Threshold specification. Either a numeric scalar/vector or a
#'   function.
#'
#' @return Numeric threshold vector of length `nrow(x$data$X_original)`.
#'
#' @keywords internal
spar_resolve_threshold_vector <- function(x, u) {
  n <- nrow(x$data$X_original)

  if (is.function(u)) {
    fn_formals <- names(formals(u))
    args <- list()

    if (!is.null(fn_formals)) {
      if ("x" %in% fn_formals) args$x <- x
      if ("R" %in% fn_formals) args$R <- x$angular$R
      if ("phi" %in% fn_formals) args$phi <- x$angular$phi
      if ("n" %in% fn_formals) args$n <- n
    }

    if (length(args) == 0L) {
      u_val <- u(x)
    } else {
      u_val <- do.call(u, args)
    }
  } else {
    u_val <- u
  }

  if (is.matrix(u_val) || is.data.frame(u_val)) {
    if (ncol(u_val) != 1L) {
      stop("Threshold specification must resolve to a scalar, vector, or one-column object.", call. = FALSE)
    }
    u_val <- u_val[, 1L, drop = TRUE]
  }

  u_val <- as.numeric(u_val)

  if (length(u_val) == 1L) {
    u_val <- rep(u_val, n)
  }

  if (length(u_val) != n) {
    stop(sprintf("Threshold length %d does not match number of observations %d.", length(u_val), n), call. = FALSE)
  }

  if (any(!is.finite(u_val))) {
    stop("Threshold vector contains non-finite values.", call. = FALSE)
  }

  u_val
}

#' Update excess fields from radial values and thresholds
#'
#' @param x A `spar_representation` object.
#'
#' @return Updated `spar_representation`.
#'
#' @keywords internal
spar_update_excess_from_threshold <- function(x) {
  if (is.null(x$angular$R)) {
    stop("Cannot compute excesses because `angular$R` is not available.", call. = FALSE)
  }

  u <- x$threshold$per_observation
  if (is.null(u)) {
    stop("Cannot compute excesses because `threshold$per_observation` is not set.", call. = FALSE)
  }
  excess <- pmax(x$angular$R - u, 0)
  x$excess$value <- excess
  x$excess$is_exceedance <- excess > 0
  x$excess$level <- NULL
  x$excess$threshold_name <- spar_threshold_active_get(x)

  x
}

#' Set per-observation thresholds on a SPAR representation
#'
#' @param x A `spar_representation` object.
#' @param u Threshold specification as numeric scalar/vector or function.
#' @param name Name under which the threshold function is stored.
#' @param set_active Logical; whether to set `name` as active threshold.
#' @param compute_excess Logical; if `TRUE`, compute excesses from
#'   `x$angular$R - threshold`.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_set_threshold <- function(
    x,
    u,
    name = "manual",
    set_active = TRUE,
    compute_excess = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }

  u_vec <- spar_resolve_threshold_vector(x, u)
  x$threshold$per_observation <- u_vec

  x$threshold$functions[[name]] <- u

  if (isTRUE(set_active)) {
    x <- spar_threshold_active_set(x, threshold_id = name)
    x <- spar_threshold_registry_upsert(
      x,
      threshold_id = name,
      estimator = list(kind = "manual", storage = "manual")
    )
  }

  if (isTRUE(compute_excess)) {
    x <- spar_update_excess_from_threshold(x)
  }

  x
}

#' Register a threshold estimator
#'
#' @param x A `spar_representation` object.
#' @param name Estimator name.
#' @param estimator Estimator object.
#' @param predict_fun Function of `(estimator, x)` returning threshold vector.
#' @param set_active Logical; whether to mark this estimator active.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_register_threshold_estimator <- function(
    x,
    name,
    estimator,
    predict_fun,
    set_active = FALSE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.function(predict_fun)) {
    stop("`predict_fun` must be a function.", call. = FALSE)
  }

  x$threshold$estimators[[name]] <- estimator
  x$threshold$functions[[name]] <- function(obj) predict_fun(estimator, obj)

  if (isTRUE(set_active)) {
    x <- spar_threshold_active_set(x, threshold_id = name)
    x <- spar_threshold_registry_upsert(x, threshold_id = name, estimator = estimator)
  }

  x
}

spar_extract_threshold_prediction <- function(pred) {
  if (is.data.frame(pred) || is.matrix(pred)) {
    pred_names <- colnames(pred)
    if (!is.null(pred_names) && "location" %in% pred_names) {
      return(as.numeric(pred[, "location", drop = TRUE]))
    }
    return(as.numeric(pred[, 1L, drop = TRUE]))
  }

  as.numeric(pred)
}

## Threshold registry helpers -------------------------------------------------

#' Create an empty threshold registry table
#'
#' @return Empty data frame with canonical threshold-registry columns.
#' @keywords internal
spar_threshold_registry_empty <- function() {
  data.frame(
    representation_id = character(0),
    transform_id = character(0),
    angular_id = character(0),
    threshold_id = character(0),
    kind = character(0),
    storage = character(0),
    tau = numeric(0),
    fit_path = character(0),
    updated_at = as.POSIXct(character(0), tz = "UTC"),
    stringsAsFactors = FALSE
  )
}

#' Normalize threshold registry structure
#'
#' Ensures `x$threshold$registry` exists and has the expected schema.
#'
#' @param x A `spar_representation` object.
#'
#' @return Updated `spar_representation`.
#' @keywords internal
spar_threshold_registry_normalize <- function(x) {
  reg <- x$threshold$registry
  if (is.null(reg)) {
    x$threshold$registry <- spar_threshold_registry_empty()
    return(x)
  }

  if (!is.data.frame(reg)) {
    stop("`threshold$registry` must be a data.frame.", call. = FALSE)
  }

  ref <- spar_threshold_registry_empty()
  for (nm in names(ref)) {
    if (!(nm %in% names(reg))) {
      reg[[nm]] <- ref[[nm]]
    }
  }

  reg <- reg[, names(ref), drop = FALSE]
  reg$representation_id <- as.character(reg$representation_id)
  reg$transform_id <- as.character(reg$transform_id)
  reg$angular_id <- as.character(reg$angular_id)
  reg$threshold_id <- as.character(reg$threshold_id)
  reg$kind <- as.character(reg$kind)
  reg$storage <- as.character(reg$storage)
  reg$tau <- as.numeric(reg$tau)
  reg$fit_path <- as.character(reg$fit_path)
  reg$updated_at <- as.POSIXct(reg$updated_at, tz = "UTC")

  x$threshold$registry <- reg
  x
}

#' Get current active transform identifier
#'
#' Returns a stable identifier for the currently active transform chain.
#'
#' @param x A `spar_representation` object.
#'
#' @return A single character identifier.
#' @keywords internal
spar_active_transform_identifier <- function(x) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (exists("spar_transform_registry_normalize", mode = "function")) {
    x <- spar_transform_registry_normalize(x)
  }

  active <- x$transform$active %||% NULL
  if (is.character(active) && length(active) == 1L && !is.na(active) && nzchar(active)) {
    return(active)
  }

  chain <- x$transform$chain
  if (is.list(chain)) {
    nm <- chain$name %||% NULL
    if (is.character(nm) && length(nm) == 1L && !is.na(nm) && nzchar(nm)) {
      return(nm)
    }
  }

  if (length(x$transform$steps) == 0L) {
    return("identity")
  }

  "active"
}

#' Get current active representation identifier
#'
#' Returns the composite identifier used by threshold active-map keys and
#' threshold registry rows: `transform=<id>|angular=<id>`.
#'
#' @param x A `spar_representation` object.
#'
#' @return A single character identifier.
#' @export
spar_active_representation_identifier <- function(x) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  transform_id <- spar_active_transform_identifier(x)
  angular_id <- x$angular$active %||% "active"
  if (!is.character(angular_id) || length(angular_id) != 1L || is.na(angular_id) || !nzchar(angular_id)) {
    angular_id <- "active"
  }

  sprintf("transform=%s|angular=%s", transform_id, angular_id)
}

#' Build a composite representation identifier
#'
#' @param transform_id Transform identifier.
#' @param angular_id Angular-map identifier.
#'
#' @return A single composite representation identifier string.
#' @export
spar_representation_id <- function(transform_id, angular_id) {
  if (!is.character(transform_id) || length(transform_id) != 1L || is.na(transform_id) || !nzchar(transform_id)) {
    stop("`transform_id` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(angular_id) || length(angular_id) != 1L || is.na(angular_id) || !nzchar(angular_id)) {
    stop("`angular_id` must be a single non-empty character string.", call. = FALSE)
  }

  sprintf("transform=%s|angular=%s", transform_id, angular_id)
}

spar_parse_representation_identifier <- function(representation_id) {
  if (!is.character(representation_id) || length(representation_id) != 1L || is.na(representation_id) || !nzchar(representation_id)) {
    stop("`representation_id` must be a single non-empty character string.", call. = FALSE)
  }

  m <- regexec("^transform=(.+)\\|angular=(.+)$", representation_id)
  parts <- regmatches(representation_id, m)[[1L]]
  if (length(parts) != 3L) {
    stop("`representation_id` must have format 'transform=<id>|angular=<id>'.", call. = FALSE)
  }

  list(
    transform_id = parts[2L],
    angular_id = parts[3L],
    representation_id = representation_id
  )
}

spar_predict_threshold_by_name_phi <- function(x, name, phi, domain = NULL) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }

  phi <- as.numeric(phi)
  if (length(phi) == 0L || any(!is.finite(phi))) {
    stop("`phi` must be a non-empty numeric vector of finite values.", call. = FALSE)
  }

  dom <- domain
  if (!is.null(dom)) {
    dom <- as_spar_angle_domain(dom)
    phi <- spar_normalize_angle(phi, dom)
  }

  est <- x$threshold$estimators[[name]]
  u <- NULL

  if (is.list(est) && identical(est$kind, "evgam_ald") && !is.null(est$fit)) {
    pred <- stats::predict(est$fit, newdata = data.frame(phi = phi), type = "response")
    u <- spar_extract_threshold_prediction(pred)
  }

  if (is.null(u) && is.list(est) && (identical(est$kind, "evgam_ald_compact") || identical(est$kind, "evgam_ald_hybrid") || identical(est$kind, "transported_threshold"))) {
    u <- spar_predict_compact_threshold(estimator = est, phi = phi, domain = dom %||% est$domain)
  }

  if (is.null(u)) {
    fn <- x$threshold$functions[[name]]
    if (!is.function(fn)) {
      stop(sprintf("Threshold '%s' does not expose angle-based prediction.", name), call. = FALSE)
    }

    fml <- names(formals(fn))
    if (is.null(fml) || !("phi" %in% fml)) {
      stop(sprintf("Threshold '%s' does not expose a `phi` argument.", name), call. = FALSE)
    }

    args <- list(phi = phi)
    if ("x" %in% fml) args$x <- x
    if ("n" %in% fml) args$n <- length(phi)
    u <- do.call(fn, args)
  }

  u <- as.numeric(u)
  if (length(u) == 1L) u <- rep(u, length(phi))
  if (length(u) != length(phi)) {
    stop("Predicted threshold length does not match `phi` length.", call. = FALSE)
  }
  if (any(!is.finite(u))) {
    stop("Predicted threshold contains non-finite values.", call. = FALSE)
  }

  u
}

#' Normalize threshold active-map structure
#'
#' Migrates legacy scalar active-threshold state to the representation-keyed
#' map format used by `threshold$active`.
#'
#' @param x A `spar_representation` object.
#'
#' @return Updated `spar_representation`.
#' @keywords internal
spar_threshold_active_normalize <- function(x) {
  active <- x$threshold$active

  if (is.null(active)) {
    x$threshold$active <- list()
    return(x)
  }

  if (is.character(active) && length(active) == 1L) {
    rep_id <- spar_active_representation_identifier(x)
    x$threshold$active <- stats::setNames(list(if (is.na(active) || !nzchar(active)) NULL else active), rep_id)
    return(x)
  }

  if (is.list(active)) {
    nm <- names(active)
    if (is.null(nm)) {
      nm <- rep("", length(active))
    }
    idx_empty <- which(is.na(nm) | !nzchar(nm))
    if (length(idx_empty) > 0L) {
      rep_id <- spar_active_representation_identifier(x)
      nm[idx_empty] <- paste0(rep_id, "#", seq_along(idx_empty))
      names(active) <- nm
    }
    x$threshold$active <- active
    return(x)
  }

  stop("`threshold$active` must be NULL, a character scalar, or a named list.", call. = FALSE)
}

#' Get active threshold identifier for a representation
#'
#' @param x A `spar_representation` object.
#' @param representation_id Optional composite representation identifier. If
#'   `NULL`, uses [spar_active_representation_identifier()].
#'
#' @return Active threshold identifier for the representation, or `NULL`.
#' @keywords internal
spar_threshold_active_get <- function(x, representation_id = NULL) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  x <- spar_threshold_active_normalize(x)
  if (is.null(representation_id)) {
    representation_id <- spar_active_representation_identifier(x)
  }

  if (!is.character(representation_id) || length(representation_id) != 1L || is.na(representation_id) || !nzchar(representation_id)) {
    stop("`representation_id` must be a single non-empty character string.", call. = FALSE)
  }

  id <- x$threshold$active[[representation_id]]
  if (is.null(id)) return(NULL)
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) return(NULL)
  id
}

#' Set active threshold identifier for a representation
#'
#' @param x A `spar_representation` object.
#' @param threshold_id Threshold identifier, or `NULL` to clear.
#' @param representation_id Optional composite representation identifier. If
#'   `NULL`, uses [spar_active_representation_identifier()].
#'
#' @return Updated `spar_representation`.
#' @keywords internal
spar_threshold_active_set <- function(x, threshold_id = NULL, representation_id = NULL) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  x <- spar_threshold_active_normalize(x)
  if (is.null(representation_id)) {
    representation_id <- spar_active_representation_identifier(x)
  }

  if (!is.character(representation_id) || length(representation_id) != 1L || is.na(representation_id) || !nzchar(representation_id)) {
    stop("`representation_id` must be a single non-empty character string.", call. = FALSE)
  }

  if (is.null(threshold_id)) {
    x$threshold$active[[representation_id]] <- NULL
    return(x)
  }

  if (!is.character(threshold_id) || length(threshold_id) != 1L || is.na(threshold_id) || !nzchar(threshold_id)) {
    stop("`threshold_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  x$threshold$active[[representation_id]] <- threshold_id
  x
}

#' Upsert threshold registry row for active representation
#'
#' Writes a registry row keyed by (`representation_id`, `threshold_id`). If a
#' matching row already exists, it is updated; otherwise a new row is appended.
#'
#' @param x A `spar_representation` object.
#' @param threshold_id Threshold identifier.
#' @param estimator Optional estimator object used to populate metadata.
#' @param representation_id Optional representation id override.
#' @param transform_id Optional transform id override.
#' @param angular_id Optional angular id override.
#'
#' @return Updated `spar_representation`.
#' @keywords internal
spar_threshold_registry_upsert <- function(
    x,
    threshold_id,
    estimator = NULL,
    representation_id = NULL,
    transform_id = NULL,
    angular_id = NULL
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.character(threshold_id) || length(threshold_id) != 1L || is.na(threshold_id) || !nzchar(threshold_id)) {
    stop("`threshold_id` must be a single non-empty character string.", call. = FALSE)
  }

  x <- spar_threshold_registry_normalize(x)
  rep_id <- representation_id %||% spar_active_representation_identifier(x)
  if (!is.character(rep_id) || length(rep_id) != 1L || is.na(rep_id) || !nzchar(rep_id)) {
    stop("`representation_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  parsed <- spar_parse_representation_identifier(rep_id)
  transform_id_use <- transform_id %||% parsed$transform_id
  angular_id_use <- angular_id %||% parsed$angular_id
  if (!is.character(transform_id_use) || length(transform_id_use) != 1L || is.na(transform_id_use) || !nzchar(transform_id_use)) {
    stop("`transform_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(angular_id_use) || length(angular_id_use) != 1L || is.na(angular_id_use) || !nzchar(angular_id_use)) {
    stop("`angular_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  kind <- if (is.list(estimator) && !is.null(estimator$kind)) as.character(estimator$kind) else "manual"
  storage <- if (is.list(estimator) && !is.null(estimator$storage)) as.character(estimator$storage) else "manual"
  tau <- if (is.list(estimator) && !is.null(estimator$tau) && length(estimator$tau) == 1L && is.finite(estimator$tau)) as.numeric(estimator$tau) else NA_real_
  fit_path <- if (is.list(estimator) && !is.null(estimator$fit_path)) as.character(estimator$fit_path) else NA_character_
  now <- Sys.time()

  row <- data.frame(
    representation_id = rep_id,
    transform_id = transform_id_use,
    angular_id = angular_id_use,
    threshold_id = threshold_id,
    kind = kind,
    storage = storage,
    tau = tau,
    fit_path = fit_path,
    updated_at = now,
    stringsAsFactors = FALSE
  )

  reg <- x$threshold$registry
  hit <- which(reg$representation_id == rep_id & reg$threshold_id == threshold_id)
  if (length(hit) > 0L) {
    reg[hit[1], names(row)] <- row[1, names(row), drop = FALSE]
    if (length(hit) > 1L) {
      reg <- reg[-hit[-1L], , drop = FALSE]
    }
  } else {
    reg <- rbind(reg, row)
  }

  x$threshold$registry <- reg
  x
}

#' Return threshold registry rows
#'
#' @param x A `spar_representation` object.
#' @param representation_id Optional composite representation identifier. If
#'   `NULL`, returns all registry rows.
#'
#' @return Data frame of registry rows.
#' @export
spar_threshold_registry_get <- function(x, representation_id = NULL) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  x <- spar_threshold_registry_normalize(x)
  reg <- x$threshold$registry
  if (is.null(representation_id)) {
    return(reg)
  }

  if (!is.character(representation_id) || length(representation_id) != 1L || is.na(representation_id) || !nzchar(representation_id)) {
    stop("`representation_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  reg[reg$representation_id == representation_id, , drop = FALSE]
}

#' Transport a threshold fit across representations
#'
#' Builds a threshold estimator for a target representation by evaluating a
#' source threshold over the shared observation index and remapping to target
#' angular coordinates.
#'
#' If multiple threshold values are observed for the same target-angle bin, the
#' transported threshold uses the highest value and emits a warning.
#'
#' @param x A `spar_representation` object.
#' @param threshold_id Source threshold identifier.
#' @param from_representation_id Source representation id. Defaults to active.
#' @param to_representation_id Target representation id. Defaults to active.
#' @param new_threshold_id Identifier for the transported threshold. If `NULL`,
#'   an id is auto-generated.
#' @param phi_grid Optional target angular grid.
#' @param phi_grid_n Number of target grid points when `phi_grid` is `NULL`.
#' @param phi_tol Target-angle bin width used during transport aggregation.
#' @param u_tol Tolerance for deciding whether thresholds differ within an angle
#'   bin.
#' @param warn_multi Logical; if `TRUE`, warn when multiple threshold values are
#'   observed within any target-angle bin.
#' @param set_active Logical; if `TRUE`, set transported threshold active for
#'   `to_representation_id`.
#' @param apply Logical; if `TRUE`, apply the transported threshold immediately
#'   when target representation equals currently active representation.
#' @param compute_excess Passed to [spar_apply_threshold()] when `apply = TRUE`.
#'
#' @return Updated `spar_representation`.
#' @export
spar_transport_threshold <- function(
    x,
    threshold_id,
    from_representation_id = NULL,
    to_representation_id = NULL,
    new_threshold_id = NULL,
    phi_grid = NULL,
    phi_grid_n = 720L,
    phi_tol = 1e-3,
    u_tol = 1e-8,
    warn_multi = TRUE,
    set_active = FALSE,
    apply = FALSE,
    compute_excess = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.character(threshold_id) || length(threshold_id) != 1L || is.na(threshold_id) || !nzchar(threshold_id)) {
    stop("`threshold_id` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.numeric(phi_tol) || length(phi_tol) != 1L || is.na(phi_tol) || phi_tol <= 0) {
    stop("`phi_tol` must be a single numeric value > 0.", call. = FALSE)
  }
  if (!is.numeric(u_tol) || length(u_tol) != 1L || is.na(u_tol) || u_tol < 0) {
    stop("`u_tol` must be a single numeric value >= 0.", call. = FALSE)
  }

  if (exists("spar_transform_registry_normalize", mode = "function")) {
    x <- spar_transform_registry_normalize(x)
  }
  if (exists("spar_angular_registry_normalize", mode = "function")) {
    x <- spar_angular_registry_normalize(x)
  }

  from_representation_id <- from_representation_id %||% spar_active_representation_identifier(x)
  to_representation_id <- to_representation_id %||% spar_active_representation_identifier(x)
  from <- spar_parse_representation_identifier(from_representation_id)
  to <- spar_parse_representation_identifier(to_representation_id)

  if (!(threshold_id %in% names(x$threshold$estimators)) && !(threshold_id %in% names(x$threshold$functions))) {
    stop(sprintf("Threshold '%s' is not registered.", threshold_id), call. = FALSE)
  }

  if (is.null(new_threshold_id)) {
    new_threshold_id <- sprintf("%s__to__%s__%s", threshold_id, to$transform_id, to$angular_id)
  }
  if (!is.character(new_threshold_id) || length(new_threshold_id) != 1L || is.na(new_threshold_id) || !nzchar(new_threshold_id)) {
    stop("`new_threshold_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  source_spec <- x$angular$maps[[from$angular_id]]
  if (is.null(source_spec)) {
    stop(sprintf("Source angular map '%s' is not registered.", from$angular_id), call. = FALSE)
  }
  target_spec <- x$angular$maps[[to$angular_id]]
  if (is.null(target_spec)) {
    stop(sprintf("Target angular map '%s' is not registered.", to$angular_id), call. = FALSE)
  }

  source_data <- spar_compute_angular_map_data(
    x = x,
    spec = source_spec,
    source = source_spec$source,
    transform_id = from$transform_id
  )
  target_data <- spar_compute_angular_map_data(
    x = x,
    spec = target_spec,
    source = target_spec$source,
    transform_id = to$transform_id
  )

  if (length(source_data$phi) != length(target_data$phi)) {
    stop("Source and target representations do not share observation length.", call. = FALSE)
  }

  u_source <- spar_predict_threshold_by_name_phi(
    x = x,
    name = threshold_id,
    phi = source_data$phi,
    domain = source_spec$domain
  )

  domain_to <- as_spar_angle_domain(target_spec$domain)
  phi_target <- spar_normalize_angle(as.numeric(target_data$phi), domain_to)

  if (identical(domain_to$type, "cyclical")) {
    period <- domain_to$upper - domain_to$lower
    phi_base <- ((phi_target - domain_to$lower) %% period) + domain_to$lower
    bin_id <- as.integer(floor((phi_base - domain_to$lower) / phi_tol))
  } else {
    phi_base <- phi_target
    bin_id <- as.integer(floor((phi_base - domain_to$lower) / phi_tol))
  }

  split_idx <- split(seq_along(phi_target), bin_id)
  agg_rows <- lapply(names(split_idx), function(bn) {
    idx <- split_idx[[bn]]
    u_vals <- u_source[idx]
    phi_vals <- phi_target[idx]

    if (u_tol > 0) {
      u_disc <- round(u_vals / u_tol) * u_tol
    } else {
      u_disc <- u_vals
    }

    data.frame(
      bin_id = as.integer(bn),
      phi_center = stats::median(phi_vals),
      u_max = max(u_vals, na.rm = TRUE),
      n_obs = length(idx),
      n_unique_u = length(unique(u_disc)),
      u_spread = diff(range(u_vals, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  agg <- do.call(rbind, agg_rows)
  agg <- agg[order(agg$phi_center), , drop = FALSE]

  amb <- agg[agg$n_unique_u > 1L, , drop = FALSE]
  if (isTRUE(warn_multi) && nrow(amb) > 0L) {
    warning(
      sprintf(
        paste0(
          "Threshold transport detected %d target-angle bins with multiple threshold values. ",
          "Using highest threshold per angle bin (max spread %.6g)."
        ),
        nrow(amb),
        max(amb$u_spread, na.rm = TRUE)
      ),
      call. = FALSE
    )
  }

  phi_grid_use <- spar_prepare_phi_grid(
    domain = domain_to,
    phi_observed = phi_target,
    phi_grid = phi_grid,
    phi_grid_n = phi_grid_n
  )

  if (nrow(agg) == 1L) {
    u_grid <- rep(agg$u_max[[1L]], length(phi_grid_use))
  } else {
    u_grid <- spar_predict_threshold_from_grid(
      phi = phi_grid_use,
      phi_grid = agg$phi_center,
      u_grid = agg$u_max,
      domain = domain_to
    )
  }

  source_est <- x$threshold$estimators[[threshold_id]]
  tau <- if (is.list(source_est) && !is.null(source_est$tau) && length(source_est$tau) == 1L && is.finite(source_est$tau)) as.numeric(source_est$tau) else NA_real_

  estimator <- list(
    kind = "transported_threshold",
    storage = "compact",
    tau = tau,
    domain = domain_to,
    grid = list(phi = phi_grid_use, u = as.numeric(u_grid)),
    smooth = list(method = "none", model = NULL),
    source_threshold_id = threshold_id,
    from_representation_id = from_representation_id,
    to_representation_id = to_representation_id,
    aggregation = list(
      method = "max",
      phi_tol = phi_tol,
      u_tol = u_tol,
      n_bins = nrow(agg),
      n_ambiguous_bins = nrow(amb),
      max_ambiguous_spread = if (nrow(amb) > 0L) max(amb$u_spread, na.rm = TRUE) else 0,
      detail = agg
    )
  )

  threshold_fun <- function(x, phi = x$angular$phi) {
    ref_domain <- x$angular$domain %||% estimator$domain
    spar_predict_compact_threshold(
      estimator = estimator,
      phi = as.numeric(phi),
      domain = ref_domain
    )
  }

  x$threshold$estimators[[new_threshold_id]] <- estimator
  x$threshold$functions[[new_threshold_id]] <- threshold_fun
  x <- spar_threshold_registry_upsert(
    x,
    threshold_id = new_threshold_id,
    estimator = estimator,
    representation_id = to_representation_id,
    transform_id = to$transform_id,
    angular_id = to$angular_id
  )

  if (isTRUE(set_active)) {
    x <- spar_threshold_active_set(
      x,
      threshold_id = new_threshold_id,
      representation_id = to_representation_id
    )
  }

  if (isTRUE(apply)) {
    active_rep <- spar_active_representation_identifier(x)
    if (identical(active_rep, to_representation_id)) {
      x <- spar_apply_threshold(
        x,
        name = new_threshold_id,
        set_active = isTRUE(set_active),
        compute_excess = compute_excess
      )
    } else {
      warning(
        "`apply = TRUE` requested, but target representation is not active. Threshold registered but not applied.",
        call. = FALSE
      )
    }
  }

  x
}

#' Check whether an estimator stores a full fit
#'
#' @param est Estimator object.
#'
#' @return Logical scalar.
#' @keywords internal
spar_is_full_threshold_estimator <- function(est) {
  if (!is.list(est)) {
    return(FALSE)
  }

  storage <- est$storage %||% if (identical(est$kind, "evgam_ald")) "full" else NULL
  identical(storage, "full")
}

#' Build a default angular grid for compact threshold storage
#'
#' @param domain Optional angular domain.
#' @param phi_observed Observed angular values.
#' @param n Number of grid points.
#'
#' @return Numeric vector of grid angles.
#' @keywords internal
spar_default_phi_grid <- function(domain, phi_observed, n = 720L) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 8) {
    stop("`phi_grid_n` must be a single numeric value >= 8.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.null(domain)) {
    domain <- as_spar_angle_domain(domain)
    return(as.numeric(spar_angle_grid(n = n, domain = domain)))
  }

  rng <- range(phi_observed[is.finite(phi_observed)], na.rm = TRUE)
  if (length(rng) != 2L || any(!is.finite(rng))) {
    stop("Cannot construct default phi grid from observed angles.", call. = FALSE)
  }
  as.numeric(seq(rng[1], rng[2], length.out = n))
}

#' Validate and normalize user-supplied phi grid
#'
#' @param domain Optional angular domain.
#' @param phi_observed Observed angular values.
#' @param phi_grid Optional user-supplied grid.
#' @param phi_grid_n Grid size when `phi_grid` is `NULL`.
#'
#' @return Sorted unique numeric phi grid.
#' @keywords internal
spar_prepare_phi_grid <- function(domain, phi_observed, phi_grid = NULL, phi_grid_n = 720L) {
  if (is.null(phi_grid)) {
    phi_grid <- spar_default_phi_grid(domain = domain, phi_observed = phi_observed, n = phi_grid_n)
  } else {
    if (is.list(phi_grid)) {
      phi_grid <- unlist(phi_grid, use.names = FALSE)
    }
    phi_grid <- as.numeric(phi_grid)
    if (length(phi_grid) < 8L || any(!is.finite(phi_grid))) {
      stop("`phi_grid` must contain at least 8 finite values.", call. = FALSE)
    }
  }

  if (!is.null(domain)) {
    domain <- as_spar_angle_domain(domain)
    phi_grid <- spar_normalize_angle(phi_grid, domain)
  }

  phi_grid <- sort(unique(phi_grid))
  if (length(phi_grid) < 8L) {
    stop("Normalized `phi_grid` has fewer than 8 unique values.", call. = FALSE)
  }
  phi_grid
}

#' Predict threshold values from a stored grid
#'
#' @param phi Angles at which to predict.
#' @param phi_grid Stored grid angles.
#' @param u_grid Stored threshold values on `phi_grid`.
#' @param domain Optional angular domain.
#'
#' @return Numeric threshold vector.
#' @keywords internal
spar_predict_threshold_from_grid <- function(phi, phi_grid, u_grid, domain = NULL) {
  phi <- as.numeric(phi)
  phi_grid <- as.numeric(phi_grid)
  u_grid <- as.numeric(u_grid)

  if (length(phi_grid) != length(u_grid)) {
    stop("`phi_grid` and `u_grid` must have equal length.", call. = FALSE)
  }
  if (length(phi_grid) < 2L || any(!is.finite(phi_grid)) || any(!is.finite(u_grid))) {
    stop("Grid inputs must contain at least two finite points.", call. = FALSE)
  }

  ord <- order(phi_grid)
  phi_grid <- phi_grid[ord]
  u_grid <- u_grid[ord]
  keep <- !duplicated(phi_grid)
  phi_grid <- phi_grid[keep]
  u_grid <- u_grid[keep]

  if (!is.null(domain)) {
    domain <- as_spar_angle_domain(domain)
    phi <- spar_normalize_angle(phi, domain)

    if (identical(domain$type, "cyclical")) {
      period <- domain$upper - domain$lower
      x <- c(phi_grid - period, phi_grid, phi_grid + period)
      y <- c(u_grid, u_grid, u_grid)
      return(as.numeric(stats::approx(x = x, y = y, xout = phi, rule = 2, ties = "ordered")$y))
    }
  }

  as.numeric(stats::approx(x = phi_grid, y = u_grid, xout = phi, rule = 2, ties = "ordered")$y)
}

#' Summarize approximation error between full and compact thresholds
#'
#' @param u_full Full-fit threshold values.
#' @param u_approx Compact/surrogate threshold values.
#'
#' @return One-row data frame of error metrics.
#' @keywords internal
spar_threshold_approx_summary <- function(u_full, u_approx) {
  err <- as.numeric(u_approx) - as.numeric(u_full)
  abs_err <- abs(err)
  data.frame(
    mae = mean(abs_err, na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    max_abs = max(abs_err, na.rm = TRUE),
    q50_abs = as.numeric(stats::quantile(abs_err, probs = 0.50, na.rm = TRUE, names = FALSE)),
    q90_abs = as.numeric(stats::quantile(abs_err, probs = 0.90, na.rm = TRUE, names = FALSE)),
    q99_abs = as.numeric(stats::quantile(abs_err, probs = 0.99, na.rm = TRUE, names = FALSE)),
    stringsAsFactors = FALSE
  )
}

#' Predict using compact or hybrid threshold estimator
#'
#' @param estimator Compact or hybrid estimator object.
#' @param phi Angles at which to predict.
#' @param domain Optional angular domain.
#'
#' @return Numeric threshold vector.
#' @keywords internal
spar_predict_compact_threshold <- function(estimator, phi, domain = NULL) {
  grid <- estimator$grid
  if (!is.list(grid) || is.null(grid$phi) || is.null(grid$u)) {
    stop("Compact estimator does not contain grid data.", call. = FALSE)
  }

  ref_domain <- domain %||% estimator$domain
  smooth_model <- estimator$smooth$model
  if (!is.null(smooth_model) && inherits(smooth_model, "smooth.spline")) {
    phi_eval <- as.numeric(phi)
    if (!is.null(ref_domain)) {
      ref_domain <- as_spar_angle_domain(ref_domain)
      phi_eval <- spar_normalize_angle(phi_eval, ref_domain)
    }
    return(as.numeric(stats::predict(smooth_model, x = phi_eval)$y))
  }

  spar_predict_threshold_from_grid(
    phi = as.numeric(phi),
    phi_grid = as.numeric(grid$phi),
    u_grid = as.numeric(grid$u),
    domain = ref_domain
  )
}

#' Configure threshold fit size warnings
#'
#' Controls whether warnings about repeatedly storing full EVGAM threshold fits
#' are suppressed.
#'
#' @param x A `spar_representation` object.
#' @param suppress Logical; if `TRUE`, suppress size warnings.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_set_size_warning_suppression <- function(x, suppress = TRUE) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.logical(suppress) || length(suppress) != 1L || is.na(suppress)) {
    stop("`suppress` must be a single non-missing logical value.", call. = FALSE)
  }

  x$meta$suppress_size_warnings <- isTRUE(suppress)
  x
}

#' Fit ALD threshold model with evgam
#'
#' Fits an angle-conditional ALD model for radial thresholds and registers it as
#' a threshold estimator.
#'
#' @param x A `spar_representation` object with angular fields available.
#' @param name Estimator name.
#' @param formula Model formula specification passed to `evgam::evgam()`. If
#'   `NULL`, a default is constructed as
#'   `list(R ~ s(phi, bs = <auto>, k = k), ~ 1)` where `<auto>` is cyclic
#'   (`"cc"`) for cyclical angle domains and thin-plate (`"tp"`) otherwise.
#' @param tau Target quantile level.
#' @param k Basis dimension for the default location smooth when `formula` is
#'   `NULL`.
#' @param trace Trace level forwarded to `evgam::evgam()`. If `NULL`, defaults
#'   to `1` when `verbose = TRUE` and `0` otherwise.
#' @param verbose Logical; whether to print start/finish timing messages.
#' @param ald.args Optional list of ALD family arguments passed to
#'   `evgam::evgam()`. If `tau` is not supplied in `ald.args`, it is set from
#'   `tau`.
#' @param ... Additional arguments forwarded to `evgam::evgam()`.
#' @param set_active Logical; whether to set this estimator active.
#' @param storage Storage mode for estimator internals. `"full"` keeps the
#'   full `evgam` fit in memory (default), `"compact"` stores only an angular
#'   grid representation, and `"hybrid"` stores compact in memory plus the full
#'   fit on disk.
#' @param hybrid_path Optional directory path used when `storage = "hybrid"`.
#'   Defaults to `.spar_cache/thresholds/evgam` under `getwd()`.
#' @param phi_grid Optional numeric vector (or list-like collection) of angular
#'   grid points for compact/hybrid storage.
#' @param phi_grid_n Number of grid points used when `phi_grid` is `NULL` for
#'   compact/hybrid storage.
#' @param compact_smoother Optional compact prediction smoother. One of
#'   `"none"` (piecewise linear interpolation on stored grid) or `"spline"`.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_fit_threshold_ald_evgam <- function(
    x,
    name = "ald_radial",
    formula = NULL,
    tau = 0.95,
    k = 20,
    trace = NULL,
    verbose = interactive(),
    ald.args = list(),
    ...,
    set_active = TRUE,
    storage = c("full", "compact", "hybrid"),
    hybrid_path = NULL,
    phi_grid = NULL,
    phi_grid_n = 720L,
    compact_smoother = c("none", "spline")
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (is.null(x$angular$R) || is.null(x$angular$phi)) {
    stop("Angular representation is required: `angular$R` and `angular$phi` must be available.", call. = FALSE)
  }

  if (!requireNamespace("evgam", quietly = TRUE)) {
    stop("Package 'evgam' is required for ALD threshold fitting.", call. = FALSE)
  }

  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("`tau` must be a single numeric value in (0, 1).", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 3) {
    stop("`k` must be a single numeric value >= 3.", call. = FALSE)
  }

  k <- as.integer(k)

  if (!is.null(trace)) {
    if (!is.numeric(trace) || length(trace) != 1L || is.na(trace)) {
      stop("`trace` must be NULL or a single numeric value.", call. = FALSE)
    }
    trace <- as.numeric(trace)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single non-missing logical value.", call. = FALSE)
  }

  storage <- match.arg(storage)
  compact_smoother <- match.arg(compact_smoother)

  if (!is.numeric(phi_grid_n) || length(phi_grid_n) != 1L || is.na(phi_grid_n) || phi_grid_n < 8) {
    stop("`phi_grid_n` must be a single numeric value >= 8.", call. = FALSE)
  }
  if (!is.null(hybrid_path) && (!is.character(hybrid_path) || length(hybrid_path) != 1L || is.na(hybrid_path) || !nzchar(hybrid_path))) {
    stop("`hybrid_path` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  if (identical(storage, "full")) {
    active_name <- spar_threshold_active_get(x)
    active_est <- if (!is.null(active_name)) x$threshold$estimators[[active_name]] else NULL
    suppress_warn <- isTRUE(x$meta$suppress_size_warnings)

    if (!suppress_warn && spar_is_full_threshold_estimator(active_est)) {
      warning(
        sprintf(
          paste0(
            "Active threshold '%s' already stores a full EVGAM fit. ",
            "Creating another full fit can consume substantial memory. ",
            "Use spar_set_size_warning_suppression(x, TRUE) to suppress this warning."
          ),
          active_name
        ),
        call. = FALSE
      )
    }
  }

  fit_args <- list(...)
  train_df <- data.frame(R = x$angular$R, phi = x$angular$phi)
  domain <- x$angular$domain
  if (is.null(domain)) {
    domain <- new_spar_angle_domain()
  } else {
    domain <- as_spar_angle_domain(domain)
  }

  if (is.null(formula)) {
    bs_type <- if (identical(domain$type, "cyclical")) "cc" else "tp"
    loc_txt <- sprintf("R ~ s(phi, bs = '%s', k = %d)", bs_type, k)
    formula <- list(
      stats::as.formula(loc_txt),
      stats::as.formula("~ 1")
    )
  }

  if (!(is.list(formula) || inherits(formula, "formula"))) {
    stop("`formula` must be NULL, a formula, or list of formulas.", call. = FALSE)
  }

  if (!is.list(ald.args)) {
    stop("`ald.args` must be a list.", call. = FALSE)
  }

  if ("ald.args" %in% names(fit_args)) {
    stop("Supply `ald.args` via the dedicated argument, not through `...`.", call. = FALSE)
  }

  if ("trace" %in% names(fit_args)) {
    if (!is.null(trace)) {
      stop("Supply `trace` either via the dedicated argument or `...`, not both.", call. = FALSE)
    }
    trace <- as.numeric(fit_args$trace)
    fit_args$trace <- NULL
  }

  if (is.null(trace)) {
    trace <- if (isTRUE(verbose)) 1 else 0
  }

  if (identical(domain$type, "cyclical") && !("knots" %in% names(fit_args))) {
    fit_args$knots <- list(phi = c(domain$lower, domain$upper))
  }

  ald.args <- utils::modifyList(list(tau = tau), ald.args)

  base_call <- c(
    list(
      formula = formula,
      data = train_df,
      family = "ald",
      trace = trace,
      ald.args = ald.args
    ),
    fit_args
  )

  t0 <- proc.time()[["elapsed"]]
  if (isTRUE(verbose)) {
    message(sprintf("Starting evgam ALD fit '%s' (n = %d).", name, nrow(train_df)))
  }

  fit <- tryCatch(
    do.call(evgam::evgam, base_call),
    error = function(e) {
      if (isTRUE(verbose)) {
        dt <- proc.time()[["elapsed"]] - t0
        message(sprintf("evgam ALD fit '%s' failed after %.2f seconds.", name, dt))
      }
      stop(e)
    }
  )

  if (isTRUE(verbose)) {
    dt <- proc.time()[["elapsed"]] - t0
    message(sprintf("Finished evgam ALD fit '%s' in %.2f seconds.", name, dt))
  }

  estimator <- list(
    kind = "evgam_ald",
    storage = storage,
    tau = tau,
    formula = formula,
    domain = domain,
    fit = NULL,
    fit_path = NULL,
    approximation = NULL,
    grid = NULL,
    smooth = list(method = compact_smoother, model = NULL)
  )

  if (identical(storage, "full")) {
    estimator$fit <- fit
    threshold_fun <- function(x, phi = x$angular$phi) {
      pred <- stats::predict(fit, newdata = data.frame(phi = as.numeric(phi)), type = "response")
      spar_extract_threshold_prediction(pred)
    }
  } else {
    phi_grid_use <- spar_prepare_phi_grid(
      domain = domain,
      phi_observed = train_df$phi,
      phi_grid = phi_grid,
      phi_grid_n = phi_grid_n
    )

    u_grid_full <- spar_extract_threshold_prediction(
      stats::predict(fit, newdata = data.frame(phi = phi_grid_use), type = "response")
    )
    estimator$grid <- list(phi = phi_grid_use, u = u_grid_full)

    if (identical(compact_smoother, "spline")) {
      estimator$smooth$model <- stats::smooth.spline(
        x = phi_grid_use,
        y = u_grid_full,
        all.knots = TRUE
      )
    }

    u_full_train <- spar_extract_threshold_prediction(
      stats::predict(fit, newdata = data.frame(phi = train_df$phi), type = "response")
    )
    u_approx_train <- spar_predict_compact_threshold(
      estimator = estimator,
      phi = train_df$phi,
      domain = domain
    )
    approx_summary <- spar_threshold_approx_summary(
      u_full = u_full_train,
      u_approx = u_approx_train
    )
    estimator$approximation <- list(
      metric = "u_compact - u_full",
      n_train = nrow(train_df),
      phi_grid_n = length(phi_grid_use),
      error_summary = approx_summary
    )

    message(sprintf(
      paste0(
        "Compact approximation summary for '%s' (storage=%s): ",
        "MAE=%.6g, RMSE=%.6g, max|err|=%.6g"
      ),
      name,
      storage,
      approx_summary$mae,
      approx_summary$rmse,
      approx_summary$max_abs
    ))

    if (identical(storage, "hybrid")) {
      target_dir <- hybrid_path %||% file.path(getwd(), ".spar_cache", "thresholds", "evgam")
      dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      pid <- Sys.getpid()
      rnd <- sprintf("%06d", sample.int(1e6, 1L) - 1L)
      fit_file <- file.path(target_dir, sprintf("%s_%s_pid%s_%s.rds", name, stamp, pid, rnd))
      saveRDS(fit, fit_file)
      estimator$fit_path <- normalizePath(fit_file, winslash = "/", mustWork = FALSE)
      estimator$kind <- "evgam_ald_hybrid"
    } else {
      estimator$kind <- "evgam_ald_compact"
    }

    threshold_fun <- function(x, phi = x$angular$phi) {
      ref_domain <- x$angular$domain %||% estimator$domain
      spar_predict_compact_threshold(
        estimator = estimator,
        phi = as.numeric(phi),
        domain = ref_domain
      )
    }
  }

  x$threshold$estimators[[name]] <- estimator
  x$threshold$functions[[name]] <- threshold_fun

  if (isTRUE(set_active)) {
    x <- spar_threshold_active_set(x, threshold_id = name)
    x <- spar_threshold_registry_upsert(x, threshold_id = name, estimator = estimator)
  }

  x
}

#' Fit and optionally apply a threshold estimator
#'
#' High-level wrapper for fitting threshold estimators on a
#' `spar_representation`.
#'
#' @param x A `spar_representation` object.
#' @param method Threshold fitting method. Currently supports `"evgam_ald"`.
#' @param name Estimator name.
#' @param apply Logical; if `TRUE`, applies the fitted threshold immediately.
#' @param set_active Logical; whether to set the estimator as active.
#' @param compute_excess Logical; if `TRUE` and `apply = TRUE`, update excess
#'   fields.
#' @param ... Additional arguments forwarded to the method-specific fit
#'   function.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_fit_threshold <- function(
    x,
    method = c("evgam_ald"),
    name = "ald_radial",
    apply = TRUE,
    set_active = TRUE,
    compute_excess = TRUE,
    ...
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  method <- match.arg(method)

  x <- switch(
    method,
    evgam_ald = spar_fit_threshold_ald_evgam(
      x = x,
      name = name,
      set_active = set_active,
      ...
    )
  )

  if (isTRUE(apply)) {
    x <- spar_apply_threshold(
      x,
      name = name,
      set_active = set_active,
      compute_excess = compute_excess
    )
  }

  x
}

#' Predict threshold values at supplied or grid angles
#'
#' @param x A `spar_representation` object.
#' @param phi Optional numeric vector of angles at which to predict threshold.
#'   If `NULL`, a grid is generated.
#' @param n Number of grid points when `phi` is `NULL`.
#' @param name Threshold estimator/function name. Defaults to active threshold.
#'
#' @return A data frame with columns `phi` and `threshold`.
#'
#' @export
spar_predict_threshold <- function(
    x,
    phi = NULL,
    n = 200L,
    name = NULL
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (is.null(name)) {
    name <- spar_threshold_active_get(x)
  }

  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("No threshold name provided and no active threshold is set.", call. = FALSE)
  }

  if (is.null(phi)) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 2) {
      stop("`n` must be a single numeric value >= 2.", call. = FALSE)
    }

    n <- as.integer(n)

    domain <- x$angular$domain
    if (!is.null(domain)) {
      domain <- as_spar_angle_domain(domain)
      phi <- spar_angle_grid(n = n, domain = domain)
    } else if (!is.null(x$angular$phi)) {
      rng <- range(x$angular$phi[is.finite(x$angular$phi)])
      phi <- seq(rng[1], rng[2], length.out = n)
    } else {
      stop("Cannot build angle grid because neither `angular$domain` nor `angular$phi` is available.", call. = FALSE)
    }
  } else {
    phi <- as.numeric(phi)
    if (length(phi) == 0L || any(!is.finite(phi))) {
      stop("`phi` must be a non-empty numeric vector of finite values.", call. = FALSE)
    }
  }

  domain <- x$angular$domain
  if (!is.null(domain)) {
    domain <- as_spar_angle_domain(domain)
    phi <- spar_normalize_angle(phi, domain)
  }

  est <- x$threshold$estimators[[name]]
  u <- NULL

  if (is.list(est) && identical(est$kind, "evgam_ald") && !is.null(est$fit)) {
    pred <- stats::predict(est$fit, newdata = data.frame(phi = phi), type = "response")
    u <- spar_extract_threshold_prediction(pred)
  }

  if (is.null(u) && is.list(est) && (identical(est$kind, "evgam_ald_compact") || identical(est$kind, "evgam_ald_hybrid") || identical(est$kind, "transported_threshold"))) {
    u <- spar_predict_compact_threshold(
      estimator = est,
      phi = phi,
      domain = domain %||% est$domain
    )
  }

  if (is.null(u)) {
    fn <- x$threshold$functions[[name]]
    if (!is.function(fn)) {
      stop(sprintf("Threshold function '%s' is not registered.", name), call. = FALSE)
    }

    fml <- names(formals(fn))
    if (!is.null(fml) && "phi" %in% fml) {
      args <- list(phi = phi)
      if ("x" %in% fml) args$x <- x
      if ("n" %in% fml) args$n <- length(phi)
      u <- do.call(fn, args)
    } else {
      stop(
        sprintf(
          "Threshold '%s' does not expose angle-based prediction. Use an evgam ALD estimator or a function with a `phi` argument.",
          name
        ),
        call. = FALSE
      )
    }
  }

  u <- as.numeric(u)
  if (length(u) == 1L) {
    u <- rep(u, length(phi))
  }

  if (length(u) != length(phi)) {
    stop("Predicted threshold length does not match `phi` length.", call. = FALSE)
  }

  data.frame(phi = phi, threshold = u)
}

#' Apply a threshold function or estimator
#'
#' @param x A `spar_representation` object.
#' @param name Name of registered threshold function/estimator. If `NULL`, uses
#'   active threshold.
#' @param set_active Logical; if `TRUE`, set `name` active.
#' @param compute_excess Logical; if `TRUE`, update excess fields.
#'
#' @return Updated `spar_representation`.
#'
#' @export
spar_apply_threshold <- function(
    x,
    name = NULL,
    set_active = TRUE,
    compute_excess = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (is.null(name)) {
    name <- spar_threshold_active_get(x)
  }

  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("No threshold name provided and no active threshold is set.", call. = FALSE)
  }

  fn <- x$threshold$functions[[name]]
  if (!is.function(fn)) {
    stop(sprintf("Threshold function '%s' is not registered.", name), call. = FALSE)
  }

  u_vec <- spar_resolve_threshold_vector(x, fn)
  x$threshold$per_observation <- u_vec

  if (isTRUE(set_active)) {
    x <- spar_threshold_active_set(x, threshold_id = name)
    x <- spar_threshold_registry_upsert(x, threshold_id = name, estimator = x$threshold$estimators[[name]])
  }

  if (isTRUE(compute_excess)) {
    x <- spar_update_excess_from_threshold(x)
  }

  x
}
