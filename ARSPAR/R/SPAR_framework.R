#' Null-coalescing helper
#'
#' Returns `y` if `x` is `NULL`, otherwise returns `x`.
#'
#' @param x An object.
#' @param y A fallback object.
#'
#' @return `x` if non-`NULL`, otherwise `y`.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Create a new SPAR representation object
#'
#' Low-level constructor for a `spar_representation` object.
#' This function creates the internal object structure and normalizes the
#' original observational data into the canonical internal form.
#'
#' In general, users should prefer [spar_representation()], which wraps this
#' constructor and optionally validates the resulting object.
#'
#' @param X_original Original observational data. Must be coercible to a numeric
#'   matrix with one row per observation.
#' @param time Optional observation-aligned time vector. If provided, must have
#'   length equal to `nrow(X_original)`.
#' @param observation_id Optional observation identifiers. If provided, must
#'   have length equal to `nrow(X_original)` and contain unique values.
#' @param original_names Optional names for the original-space variables.
#'   Defaults to `colnames(X_original)` when available, otherwise generic names
#'   are created.
#' @param description Optional free-text description of the representation.
#' @param id Optional identifier for the representation.
#'
#' @return An object of class `"spar_representation"`.
#'
#' @seealso [spar_representation()], [validate_spar_representation()],
#'   [new_spar_angle_domain()]
#'
#' @keywords internal
new_spar_representation <- function(
    X_original,
    time = NULL,
    observation_id = NULL,
    original_names = colnames(X_original),
    description = NULL,
    id = NULL
) {
  if (is.data.frame(X_original) || tibble::is_tibble(X_original)) {
    X_original <- as.matrix(X_original)
  } else if (is.vector(X_original) && is.numeric(X_original)) {
    X_original <- matrix(X_original, ncol = 1L)
  } else {
    X_original <- as.matrix(X_original)
  }
  
  if (!is.numeric(X_original)) {
    stop("`X_original` must be coercible to a numeric matrix.", call. = FALSE)
  }
  
  n <- nrow(X_original)
  d <- ncol(X_original)
  
  if (is.null(n) || is.null(d) || n == 0L || d == 0L) {
    stop("`X_original` must have positive dimensions.", call. = FALSE)
  }
  
  if (is.null(original_names)) {
    original_names <- paste0("X", seq_len(d))
  }
  
  if (length(original_names) != d) {
    stop(
      "`original_names` must have length equal to `ncol(X_original)`.",
      call. = FALSE
    )
  }
  
  colnames(X_original) <- original_names
  
  if (!is.null(time) && length(time) != n) {
    stop(
      "`time` must have length equal to `nrow(X_original)`.",
      call. = FALSE
    )
  }
  
  if (is.null(observation_id)) {
    observation_id <- seq_len(n)
  } else {
    if (length(observation_id) != n) {
      stop(
        "`observation_id` must have length equal to `nrow(X_original)`.",
        call. = FALSE
      )
    }
    
    if (anyDuplicated(observation_id)) {
      stop("`observation_id` must contain unique values.", call. = FALSE)
    }
  }
  
  identity_map <- function(X) X
  
  structure(
    list(
      meta = list(
        version = "0.1.0",
        created_at = Sys.time(),
        id = id,
        description = description
      ),
      
      schema = list(
        original_names = original_names,
        transformed_names = original_names,
        original_dim = d,
        time_name = if (!is.null(time)) "time" else NULL,
        observation_id_name = "observation_id",
        radial_name = "R",
        angle_name = "phi"
      ),
      
      data = list(
        X_original = X_original,
        time = time,
        observation_id = observation_id
      ),
      
      index = list(
        time_present = !is.null(time),
        time_ordered = if (is.null(time)) NA else NULL,
        time_order = NULL
      ),
      
      workspace = list(
        active = NULL,
        layouts = list(),
        pointers = NULL,
        state = list(
          active_transform = NULL,
          active_threshold = NULL,
          active_declustering = NULL
        )
      ),
      
        transform = list(
          active = "identity",
          chains = list(),
          steps = list(),
          forward = identity_map,
          inverse = identity_map,
          transformed_data = X_original,
        chain = NULL,
        cache_requirements = NULL,
        analysis = NULL,
        step_cache = NULL
      ),
      
        angular = list(
          X = NULL,
          R = NULL,
          phi = NULL,
          maps = list(),
          gauge = NULL,
          angle_map = NULL,
          inverse = NULL,
          domain = NULL,
          active = NULL,
          source = NULL,
          format = NULL,
          transform_id = NULL
        ),
      
      threshold = list(
        per_observation = NULL,
        by_level = list(),
        estimators = list(),
        functions = list(),
        active = list(),
        registry = data.frame(
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
      ),
      
      excess = list(
        value = NULL,
        is_exceedance = NULL,
        level = NULL,
        threshold_name = NULL
      ),
      
      excursions = list(
        pointwise = NULL,
        clusters = NULL,
        paths = NULL,
        envelopes = NULL
      ),
      
      fitted = list(
        models = list(),
        diagnostics = list()
      )
    ),
    class = "spar_representation"
  )
}


#' Construct a SPAR representation
#'
#' Creates a `spar_representation` object from original observational data and
#' optional observation metadata.
#'
#' A `spar_representation` is the central object used to store:
#' \itemize{
#'   \item original-space data,
#'   \item naming/schema information,
#'   \item observation indexing and time-order metadata,
#'   \item transformation definitions and caches,
#'   \item angular/radial representation data and angular-domain semantics,
#'   \item threshold and excess information,
#'   \item excursion-derived objects,
#'   \item and placeholders for future optimized workspace layouts.
#' }
#'
#' The constructor normalizes `X_original` to the canonical internal
#' representation as a numeric matrix with one row per observation.
#'
#' @param X_original Original observational data. Must be coercible to a numeric
#'   matrix with one row per observation.
#' @param time Optional observation-aligned time vector. If provided, must have
#'   length equal to `nrow(X_original)`.
#' @param observation_id Optional observation identifiers. If provided, must
#'   have length equal to `nrow(X_original)` and contain unique values. If not
#'   supplied, sequential identifiers are generated.
#' @param original_names Optional names for the original-space variables.
#'   Defaults to `colnames(X_original)` when available; otherwise generic names
#'   are created.
#' @param description Optional free-text description of the representation.
#' @param id Optional identifier for the representation.
#' @param validate Logical; if `TRUE`, validate the constructed object with
#'   [validate_spar_representation()].
#' @param warn_time_unordered Logical; if `TRUE`, validation may warn when a
#'   time variable is present but observations are not ordered by time.
#'
#' @details
#' The resulting object contains canonical source observations in `data`,
#' user-facing naming information in `schema`, and placeholder structures for
#' derived layers such as transformations, angular/radial coordinates,
#' thresholds, excesses, excursions, and future active workspace layouts.
#'
#' The constructor itself does not compute transformed coordinates, angular
#' coordinates, thresholds, excesses, or excursion objects. Those are intended
#' to be populated by later mutator or builder functions.
#'
#' @return An object of class `"spar_representation"`.
#'
#' @examples
#' X <- cbind(
#'   Hs = c(1.2, 2.4, 1.8),
#'   Tm = c(5.0, 6.2, 5.7)
#' )
#'
#' obj <- spar_representation(X)
#'
#' tt <- as.POSIXct(
#'   c("2026-01-01 00:00:00", "2026-01-01 01:00:00", "2026-01-01 02:00:00"),
#'   tz = "UTC"
#' )
#'
#' obj_time <- spar_representation(
#'   X_original = X,
#'   time = tt,
#'   description = "Example wave observations"
#' )
#'
#' @seealso [new_spar_representation()], [validate_spar_representation()],
#'   [new_spar_angle_domain()]
#'
#' @export
spar_representation <- function(
    X_original,
    time = NULL,
    observation_id = NULL,
    original_names = colnames(X_original),
    description = NULL,
    id = NULL,
    validate = TRUE,
    warn_time_unordered = TRUE
) {
  obj <- new_spar_representation(
    X_original = X_original,
    time = time,
    observation_id = observation_id,
    original_names = original_names,
    description = description,
    id = id
  )
  
  if (isTRUE(validate)) {
    obj <- validate_spar_representation(
      obj,
      warn_time_unordered = warn_time_unordered
    )
  }
  
  obj
}

#' Print method for spar_representation
#' @export
print.spar_representation <- function(x, ..., max_items = 5L) {
  stopifnot(inherits(x, "spar_representation"))

  n <- nrow(x$data$X_original)
  d <- ncol(x$data$X_original)

  cat("<spar_representation>\n")

  # --- meta ---
  version <- x$meta$version %||% NA_character_
  desc <- x$meta$description %||% NULL
  obj_id <- x$meta$id %||% NULL

  cat("  Version: ", version, "\n", sep = "")
  if (!is.null(obj_id)) {
    cat("  ID:      ", obj_id, "\n", sep = "")
  }
  if (!is.null(desc)) {
    cat("  Desc:    ", desc, "\n", sep = "")
  }

  # --- schema / dimensions ---
  cat("  Data:    ", n, " observations x ", d, " variables\n", sep = "")

  vars <- x$schema$original_names
  if (!is.null(vars)) {
    shown <- head(vars, max_items)
    more <- length(vars) - length(shown)
    var_text <- paste(shown, collapse = ", ")
    if (more > 0L) {
      var_text <- paste0(var_text, ", ... (+", more, " more)")
    }
    cat("  Vars:    ", var_text, "\n", sep = "")
  }

  # --- time info ---
  time <- x$data$time
  if (is.null(time)) {
    cat("  Time:    <none>\n")
  } else {
    time_msg <- if (isTRUE(x$index$time_ordered)) "ordered" else "not ordered"

    t_min <- tryCatch(min(time, na.rm = TRUE), error = function(e) NA)
    t_max <- tryCatch(max(time, na.rm = TRUE), error = function(e) NA)

    cat("  Time:    ", time_msg, "\n", sep = "")
    cat("  Range:   ", format(t_min), " -> ", format(t_max), "\n", sep = "")
  }

  # --- transform info ---
  n_steps <- length(x$transform$steps)
  active_transform <- x$transform$active %||% NULL
  n_transforms <- length(x$transform$chains %||% list())
  if (!is.null(active_transform)) {
    cat("  Active transform: ", active_transform, "\n", sep = "")
  }
  cat("  Transform registry size: ", n_transforms, "\n", sep = "")
  cat("  Transform steps: ", n_steps, "\n", sep = "")

  if (n_steps > 0L) {
    step_names <- vapply(
      x$transform$steps,
      function(s) {
        if (!is.null(s$name)) s$name else "<unnamed>"
      },
      character(1)
    )
    shown <- head(step_names, max_items)
    more <- length(step_names) - length(shown)
    step_text <- paste(shown, collapse = ", ")
    if (more > 0L) {
      step_text <- paste0(step_text, ", ... (+", more, " more)")
    }
    cat("  Steps:   ", step_text, "\n", sep = "")
  }

  # --- angular state ---
  has_angular <- !is.null(x$angular$R) && !is.null(x$angular$phi)
  n_angular_maps <- length(x$angular$maps %||% list())
  cat("  Angular: ", if (has_angular) "available" else "not computed", "\n", sep = "")
  cat("  Angular registry size: ", n_angular_maps, "\n", sep = "")

  if (!is.null(x$angular$active)) {
    ang_name <- x$angular$active %||% "<unnamed>"
    cat("  Active angular representation: ", ang_name, "\n", sep = "")
  }

  # --- threshold state ---
  n_estimators <- length(x$threshold$estimators)
  n_functions <- length(x$threshold$functions)
  has_threshold_obs <- !is.null(x$threshold$per_observation)

  cat("  Threshold estimators: ", n_estimators, "\n", sep = "")
  cat("  Threshold functions:  ", n_functions, "\n", sep = "")
  cat("  Per-observation threshold: ",
      if (has_threshold_obs) "available" else "not set",
      "\n", sep = "")

  rep_id <- tryCatch(spar_active_representation_identifier(x), error = function(e) NULL)
  active_threshold <- tryCatch(spar_threshold_active_get(x, representation_id = rep_id), error = function(e) NULL)
  if (!is.null(rep_id)) {
    cat("  Active representation id: ", rep_id, "\n", sep = "")
  }
  if (!is.null(active_threshold)) {
    cat("  Active threshold: ", active_threshold, "\n", sep = "")
  }

  # --- excess state ---
  has_excess <- !is.null(x$excess$value)
  cat("  Excesses: ", if (has_excess) "available" else "not computed", "\n", sep = "")

  if (has_excess && !is.null(x$excess$is_exceedance)) {
    n_exc <- sum(x$excess$is_exceedance, na.rm = TRUE)
    cat("  Exceedances: ", n_exc, "\n", sep = "")
  }

  # --- excursions state ---
  exc <- x$excursions
  cat("  Excursions:\n")
  cat("    pointwise: ", if (is.null(exc$pointwise)) "none" else "available", "\n", sep = "")
  cat("    clusters:  ", if (is.null(exc$clusters)) "none" else length(exc$clusters), "\n", sep = "")
  cat("    paths:     ", if (is.null(exc$paths)) "none" else length(exc$paths), "\n", sep = "")
  cat("    envelopes: ", if (is.null(exc$envelopes)) "none" else length(exc$envelopes), "\n", sep = "")

  # --- fitted objects ---
  n_models <- length(x$fitted$models)
  n_diags <- length(x$fitted$diagnostics)
  cat("  Fitted:\n")
  cat("    models:      ", n_models, "\n", sep = "")
  cat("    diagnostics: ", n_diags, "\n", sep = "")

  invisible(x)
}
