#' Extract exceedance index from a SPAR representation
#'
#' @param x A `spar_representation` object.
#' @param t Optional observation index vector. Defaults to sequential index.
#' @param R Optional radial vector. Defaults to `x$angular$R`.
#' @param phi Optional angle vector. Defaults to `x$angular$phi`.
#' @param u Optional threshold vector. Defaults to
#'   `x$threshold$per_observation`.
#' @param excess Optional excess vector. Defaults to `x$excess$value`, otherwise
#'   computed as `pmax(R - u, 0)`.
#'
#' @return An `exceedance_index` object.
#'
#' @export
spar_extract_exceedance_index <- function(
    x,
    t = NULL,
    R = NULL,
    phi = NULL,
    u = NULL,
    excess = NULL
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  n <- nrow(x$data$X_original)

  if (is.null(R)) R <- x$angular$R
  if (is.null(phi)) phi <- x$angular$phi
  if (is.null(u)) u <- x$threshold$per_observation

  if (is.null(R) || is.null(phi)) {
    stop("`angular$R` and `angular$phi` are required (or provide `R` and `phi`).", call. = FALSE)
  }

  if (is.null(u)) {
    stop("Threshold vector `u` is required (or set `x$threshold$per_observation`).", call. = FALSE)
  }

  R <- as.numeric(R)
  phi <- as.numeric(phi)
  u <- as.numeric(u)

  if (length(u) == 1L) {
    u <- rep(u, n)
  }

  if (length(R) != n || length(phi) != n || length(u) != n) {
    stop("`R`, `phi`, and `u` must all have length equal to number of observations.", call. = FALSE)
  }

  if (is.null(t)) {
    t <- seq_len(n)
  }
  t <- as.numeric(t)
  if (length(t) != n) {
    stop("`t` must have length equal to number of observations.", call. = FALSE)
  }

  if (is.null(excess)) {
    excess <- x$excess$value
    if (is.null(excess)) {
      excess <- pmax(R - u, 0)
    }
  }
  excess <- as.numeric(excess)
  if (length(excess) != n) {
    stop("`excess` must have length equal to number of observations.", call. = FALSE)
  }

  keep <- excess > 0

  new_exceedance_index(
    idx = t[keep],
    t = t[keep],
    R = R[keep],
    phi = phi[keep],
    u = u[keep],
    excess = excess[keep],
    metadata = list(source = "spar_representation")
  )
}

#' Build excursion path group from a SPAR representation
#'
#' @param x A `spar_representation` object.
#' @param clustered_spans Optional `clustered_excursion_spans`. If `NULL`, spans
#'   are derived from exceedances using `gap_rule`.
#' @param gap_rule Non-negative declustering gap used when `clustered_spans` is
#'   `NULL`.
#' @param space Data space used for path coordinates, one of `"transformed"` or
#'   `"original"`.
#' @param keep_prev_next Logical; whether to retain neighboring points in path
#'   metadata.
#' @param store Logical; if `TRUE`, returns an updated `spar_representation`
#'   with `excursions` fields populated; otherwise returns an
#'   `excursion_path_group`.
#'
#' @return An `excursion_path_group` or updated `spar_representation`.
#'
#' @export
spar_build_excursion_path_group <- function(
    x,
    clustered_spans = NULL,
    gap_rule = 0,
    space = c("transformed", "original"),
    keep_prev_next = TRUE,
    store = FALSE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  space <- match.arg(space)

  if (is.null(x$angular$R) || is.null(x$angular$phi)) {
    stop("Angular representation is required to build excursion paths.", call. = FALSE)
  }

  X <- spar_data(x, space = space, format = "matrix")
  if (is.null(X)) {
    stop(sprintf("No data available in '%s' space.", space), call. = FALSE)
  }
  X <- as.matrix(X)

  n <- nrow(X)
  u <- x$threshold$per_observation
  if (is.null(u)) {
    stop("Per-observation threshold is required (`x$threshold$per_observation`).", call. = FALSE)
  }
  u <- as.numeric(u)
  if (length(u) == 1L) u <- rep(u, n)
  if (length(u) != n) {
    stop("Threshold vector length does not match number of observations.", call. = FALSE)
  }

  R <- as.numeric(x$angular$R)
  phi <- as.numeric(x$angular$phi)

  if (length(R) != n || length(phi) != n) {
    stop("Angular vectors have invalid length relative to selected data space.", call. = FALSE)
  }

  excess <- x$excess$value
  if (is.null(excess)) excess <- pmax(R - u, 0)
  excess <- as.numeric(excess)

  t <- seq_len(n)

  exc_index <- spar_extract_exceedance_index(
    x,
    t = t,
    R = R,
    phi = phi,
    u = u,
    excess = excess
  )

  spans <- extract_exceedance_spans(exc_index, angle_domain = x$angular$domain)
  if (is.null(clustered_spans)) {
    clustered_spans <- extract_clustered_excursion_spans(
      spans,
      gap_rule = gap_rule,
      angle_domain = x$angular$domain
    )
  }

  path_df <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
  path_df$R <- R
  path_df$phi <- phi
  path_df$u <- u
  path_df$t <- t
  path_df$excess <- excess

  radial_fun <- x$angular$gauge %||% L1_radial
  angle_fun <- x$angular$angle_map %||% L1_angle

  schema <- list(
    X = colnames(X),
    R = "R",
    phi = "phi",
    u = "u",
    t = "t",
    excess = "excess",
    radial_fun = "radial_fun",
    angle_fun = "angle_fun"
  )

  group <- build_excursion_path_group(
    data = path_df,
    clustered_spans = clustered_spans,
    mapping = NULL,
    schema = schema,
    radial_fun = radial_fun,
    angle_fun = angle_fun,
    keep_prev_next = keep_prev_next,
    metadata = list(
        source = "spar_representation",
        data_space = space,
        threshold_name = spar_threshold_active_get(x)
      )
  )

  if (!isTRUE(store)) {
    return(group)
  }

  x$excursions$pointwise <- exc_index
  x$excursions$clusters <- clustered_spans
  x$excursions$paths <- group
  x
}

#' Decluster exceedances on a SPAR representation
#'
#' @param x A `spar_representation` object.
#' @param gap_rule Non-negative declustering gap.
#' @param store Logical; if `TRUE`, store exceedance and cluster objects in
#'   `x$excursions` and return updated `x`. Otherwise return a list.
#'
#' @return Updated `spar_representation` or list with `pointwise`, `spans`, and
#'   `clusters`.
#'
#' @export
spar_decluster_excursions <- function(x, gap_rule = 0, store = TRUE) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (!is.numeric(gap_rule) || length(gap_rule) != 1L || is.na(gap_rule) || gap_rule < 0) {
    stop("`gap_rule` must be a single non-negative numeric value.", call. = FALSE)
  }

  pointwise <- spar_extract_exceedance_index(x)
  spans <- extract_exceedance_spans(pointwise, angle_domain = x$angular$domain)
  clusters <- extract_clustered_excursion_spans(
    spans,
    gap_rule = gap_rule,
    angle_domain = x$angular$domain
  )

  if (!isTRUE(store)) {
    return(list(pointwise = pointwise, spans = spans, clusters = clusters))
  }

  x$excursions$pointwise <- pointwise
  x$excursions$clusters <- clusters
  x
}

#' Build excursion envelope sample from a SPAR representation
#'
#' @param x A `spar_representation` object.
#' @param path_group Optional `excursion_path_group`. If `NULL`, uses stored
#'   paths or builds paths from current representation state.
#' @param gap_rule Declustering gap used if paths must be built.
#' @param phi_grid Optional angle grid.
#' @param n_phi Grid size when `phi_grid` is `NULL`.
#' @param lambda Mixing weight for envelope interpolation.
#' @param n_sub Number of segment interpolation points.
#' @param method Envelope construction method: `"optimized"` (default) or
#'   `"brute_force"`.
#' @param space Data space used if paths must be built.
#' @param store Logical; if `TRUE`, stores envelope objects in
#'   `x$excursions$envelopes` and returns updated `x`.
#'
#' @return `envelope_sample` or updated `spar_representation`.
#'
#' @export
spar_build_excursion_envelope_sample <- function(
    x,
    path_group = NULL,
    gap_rule = 0,
    phi_grid = NULL,
    n_phi = 200L,
    lambda = 0,
    n_sub = 100L,
    method = c("optimized", "brute_force"),
    space = c("transformed", "original"),
    store = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  space <- match.arg(space)

  if (is.null(path_group)) {
    path_group <- x$excursions$paths
  }

  if (is.null(path_group)) {
    path_group <- spar_build_excursion_path_group(
      x,
      clustered_spans = NULL,
      gap_rule = gap_rule,
      space = space,
      keep_prev_next = TRUE,
      store = FALSE
    )
  }

  if (!inherits(path_group, "excursion_path_group")) {
    stop("`path_group` must be an `excursion_path_group`.", call. = FALSE)
  }

  if (is.null(phi_grid)) {
    if (!is.numeric(n_phi) || length(n_phi) != 1L || is.na(n_phi) || n_phi < 2) {
      stop("`n_phi` must be a single numeric value >= 2.", call. = FALSE)
    }

    n_phi <- as.integer(n_phi)
    dom <- x$angular$domain
    if (!is.null(dom)) {
      dom <- as_spar_angle_domain(dom)
      phi_grid <- spar_angle_grid(n = n_phi, domain = dom)
    } else {
      all_phi <- unlist(lapply(path_group$paths, function(p) p$phi), use.names = FALSE)
      all_phi <- all_phi[is.finite(all_phi)]
      if (length(all_phi) == 0L) {
        stop("Cannot infer `phi_grid`; no finite path angles available.", call. = FALSE)
      }
      phi_grid <- seq(min(all_phi), max(all_phi), length.out = n_phi)
    }
  } else {
    phi_grid <- as.numeric(phi_grid)
    if (length(phi_grid) < 2L || any(!is.finite(phi_grid))) {
      stop("`phi_grid` must be a numeric vector with at least two finite values.", call. = FALSE)
    }
  }

  method <- match.arg(method)

  env_group <- lapply(path_group$paths, function(path) {
    build_cluster_envelopes(
      path = path,
      phi_grid = phi_grid,
      lambda = lambda,
      angle_domain = x$angular$domain,
      X_prev = path$metadata$X_prev,
      X_next = path$metadata$X_next,
      n_sub = n_sub,
      method = method
    )
  })
  names(env_group) <- names(path_group$paths)

  if (length(env_group) == 0L) {
    df_all <- data.frame(cluster_id = integer(0), phi = numeric(0), excess = numeric(0), weight = numeric(0))
  } else {
    cluster_ids <- seq_along(env_group)
    df_list <- lapply(cluster_ids, function(i) {
      envelope_to_dataset(env_group[[i]]$mixed, cluster_id = i, cluster_weight = 1)
    })
    df_all <- do.call(rbind, df_list)
  }

  angle_domain <- x$angular$domain
  if (is.null(angle_domain)) {
    angle_domain <- new_spar_angle_domain(type = "interval", lower = min(phi_grid), upper = max(phi_grid))
  }

  sample <- new_envelope_sample(
    df = df_all,
    angle_domain = angle_domain,
      metadata = list(
        source = "spar_representation",
        threshold_name = spar_threshold_active_get(x),
        lambda = lambda,
        n_sub = n_sub,
        method = method
      )
  )

  if (!isTRUE(store)) {
    return(sample)
  }

  x$excursions$envelopes <- list(
    sample = sample,
    group = env_group,
    phi_grid = phi_grid,
    lambda = lambda,
    method = method
  )
  x
}

#' Compute threshold calibration diagnostic from representation excursions
#'
#' @param x A `spar_representation` object.
#' @param sample Optional `envelope_sample`. If `NULL`, uses stored sample or
#'   builds one.
#' @param target_rate Optional target exceedance rate. If `NULL` and the active
#'   threshold estimator is evgam ALD, uses `1 - tau`; otherwise defaults to
#'   `0.2`.
#' @param mode Extraction mode for angle distribution.
#' @param bandwidth Optional bandwidth for `mode = "window"`.
#' @param phi_grid Optional angle grid.
#' @param progress Logical; if `TRUE`, show a progress bar while calibration is
#'   evaluated over angles.
#' @param store Logical; if `TRUE`, stores result in `x$fitted$diagnostics` and
#'   returns updated `x`.
#'
#' @return `threshold_calibration` or updated `spar_representation`.
#'
#' @export
spar_compute_threshold_calibration <- function(
    x,
    sample = NULL,
    target_rate = NULL,
    mode = c("interpolate", "window"),
    bandwidth = NULL,
    phi_grid = NULL,
    progress = FALSE,
    store = FALSE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  mode <- match.arg(mode)

  if (is.null(sample)) {
    env <- x$excursions$envelopes
    if (is.list(env) && inherits(env$sample, "envelope_sample")) {
      sample <- env$sample
    } else {
      sample <- spar_build_excursion_envelope_sample(x, store = FALSE)
    }
  }

  if (!inherits(sample, "envelope_sample")) {
    stop("`sample` must be an `envelope_sample`.", call. = FALSE)
  }

  if (is.null(target_rate)) {
    active <- spar_threshold_active_get(x)
    est <- if (!is.null(active)) x$threshold$estimators[[active]] else NULL
    if (is.list(est) && is.finite(est$tau)) {
      target_rate <- 1 - est$tau
    } else {
      target_rate <- 0.2
    }
  }

  cal <- compute_threshold_calibration(
    sample = sample,
    phi_grid = if (is.null(phi_grid)) {
      dom <- x$angular$domain
      n_phi_default <- sample$metadata$n_phi
      if (!is.finite(n_phi_default) || length(n_phi_default) != 1L) {
        n_phi_default <- 180L
      }
      if (inherits(dom, "spar_angle_domain")) {
        spar_angle_grid(n = as.integer(n_phi_default), domain = dom)
      } else {
        NULL
      }
    } else {
      phi_grid
    },
    mode = mode,
    bandwidth = bandwidth,
    target_rate = target_rate,
    progress = progress
  )

  if (!isTRUE(store)) {
    return(cal)
  }

  x$fitted$diagnostics$threshold_calibration <- cal
  x
}

#' Compute declustered angular density
#'
#' Computes angle-wise storm support and corresponding time-normalized
#' declustered rates. By default, support is computed directly from declustered
#' excursion spans (no envelope construction).
#'
#' @param x A `spar_representation` object.
#' @param sample Optional `envelope_sample`. If `NULL`, uses stored sample or
#'   builds one (used only when `source = "envelope"`).
#' @param time_span Optional time-span denominator. If `NULL`, uses elapsed
#'   `x$data$time` when available, otherwise `nrow(x$data$X_original)`.
#' @param source Diagnostic source: `"cluster_support"` (default, no envelopes)
#'   or `"envelope"` (legacy envelope-based behavior).
#' @param clusters Optional `clustered_excursion_spans` object. Used only when
#'   `source = "cluster_support"`.
#' @param gap_rule Declustering gap used when clusters are not already available
#'   and `source = "cluster_support"`.
#' @param use_stored Logical; if `TRUE`, reuses `x$excursions$clusters` when
#'   available for `source = "cluster_support"`.
#' @param store Logical; if `TRUE`, stores result in
#'   `x$fitted$diagnostics$angular_density` and returns updated `x`.
#'
#' @return `angular_density_diagnostic` or updated `spar_representation`.
#'
#' @export
spar_compute_angular_density <- function(
    x,
    sample = NULL,
    time_span = NULL,
    phi_grid = NULL,
    n_phi = NULL,
    density_method = c("equal_point_mass", "legacy_support_rate"),
    eps = 1e-12,
    store = FALSE,
    source = c("cluster_support", "envelope"),
    clusters = NULL,
    gap_rule = 0,
    use_stored = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  source <- match.arg(source)

  if (is.null(time_span)) {
    if (!is.null(x$data$time)) {
      tvec <- x$data$time
      time_span <- suppressWarnings(as.numeric(max(tvec, na.rm = TRUE) - min(tvec, na.rm = TRUE)))
      if (!is.finite(time_span) || time_span <= 0) {
        time_span <- nrow(x$data$X_original)
      }
    } else {
      time_span <- nrow(x$data$X_original)
    }
  }

  if (!is.numeric(time_span) || length(time_span) != 1L || is.na(time_span) || time_span <= 0) {
    stop("`time_span` must be a single positive numeric value.", call. = FALSE)
  }

  if (source == "cluster_support") {
    if (is.null(x$angular$phi) || is.null(x$excess$value)) {
      stop("`x$angular$phi` and `x$excess$value` are required for `source = 'cluster_support'`.", call. = FALSE)
    }

    if (is.null(clusters) && isTRUE(use_stored)) {
      clusters <- x$excursions$clusters
    }
    if (is.null(clusters) || !inherits(clusters, "clustered_excursion_spans")) {
      dec <- spar_decluster_excursions(x, gap_rule = gap_rule, store = FALSE)
      clusters <- dec$clusters
    }

    spans <- clusters$spans

    dom <- x$angular$domain
    phi_obs <- as.numeric(x$angular$phi)
    excess <- as.numeric(x$excess$value)

    if (is.null(phi_grid)) {
      n_phi_default <- if (is.null(n_phi)) 180L else n_phi
      if (!is.finite(n_phi_default) || length(n_phi_default) != 1L || n_phi_default < 2) {
        n_phi_default <- 180L
      }

      if (inherits(dom, "spar_angle_domain")) {
        phi_grid <- spar_angle_grid(n = as.integer(n_phi_default), domain = dom)
      } else {
        p <- phi_obs[is.finite(phi_obs)]
        if (length(p) == 0L) {
          stop("Cannot infer `phi_grid`; no finite angles available.", call. = FALSE)
        }
        phi_grid <- seq(min(p), max(p), length.out = as.integer(n_phi_default))
      }
    } else {
      phi_grid <- as.numeric(phi_grid)
      if (length(phi_grid) < 2L || any(!is.finite(phi_grid))) {
        stop("`phi_grid` must be a numeric vector with at least two finite values.", call. = FALSE)
      }
    }

    if (!inherits(dom, "spar_angle_domain")) {
      dom <- new_spar_angle_domain(type = "interval", lower = min(phi_grid), upper = max(phi_grid))
    }
    dom <- as_spar_angle_domain(dom)

    nearest_bins <- function(phi_vals, grid, domain) {
      phi_vals <- as.numeric(phi_vals)
      out <- integer(length(phi_vals))
      ok <- is.finite(phi_vals)
      if (!any(ok)) return(out)

      g <- as.numeric(grid)
      if (domain$type == "cyclical") {
        g <- spar_normalize_angle(g, domain)
        p <- spar_normalize_angle(phi_vals[ok], domain)
        for (i in seq_along(p)) {
          d <- abs(p[i] - g)
          d <- pmin(d, domain$width - d)
          out[which(ok)[i]] <- which.min(d)
        }
      } else {
        mids <- 0.5 * (g[-length(g)] + g[-1L])
        brks <- c(-Inf, mids, Inf)
        out[ok] <- findInterval(phi_vals[ok], brks, rightmost.closed = TRUE)
      }

      out
    }

    n_clusters_total <- nrow(spans)
    support_count <- numeric(length(phi_grid))

    if (n_clusters_total > 0L) {
      for (i in seq_len(n_clusters_total)) {
        idx <- seq.int(spans$First[i], spans$Last[i])
        ok <- is.finite(phi_obs[idx]) & is.finite(excess[idx]) & (excess[idx] > eps)
        if (!any(ok)) next

        bins <- unique(nearest_bins(phi_obs[idx][ok], phi_grid, dom))
        bins <- bins[bins >= 1L & bins <= length(phi_grid)]
        if (length(bins) == 0L) next
        support_count[bins] <- support_count[bins] + 1
      }
    }

    conditional_support <- if (n_clusters_total > 0) support_count / n_clusters_total else rep(NA_real_, length(phi_grid))
    declustered_rate <- support_count / time_span
    w <- angular_grid_weights(phi_grid, normalize = TRUE)
    scale <- sum(declustered_rate * w, na.rm = TRUE)
    density <- if (is.finite(scale) && scale > 0) declustered_rate / scale else rep(NA_real_, length(declustered_rate))

    # On cyclical domains, lower/upper endpoints represent the same angle.
    # Keep endpoint values consistent to avoid artificial edge drops in plots.
    if (dom$type == "cyclical" && length(phi_grid) >= 2L) {
      is_wrapped_endpoint <- is.finite(dom$width) &&
        isTRUE(all.equal(phi_grid[length(phi_grid)] - phi_grid[1], dom$width, tolerance = 1e-8))
      if (is_wrapped_endpoint) {
        support_count[length(support_count)] <- support_count[1]
        conditional_support[length(conditional_support)] <- conditional_support[1]
        declustered_rate[length(declustered_rate)] <- declustered_rate[1]
        density[length(density)] <- density[1]
      }
    }

    global_rate <- if (n_clusters_total > 0) n_clusters_total / time_span else 0

    out <- new_angular_density_diagnostic(
      phi_grid = phi_grid,
      support_count = support_count,
      conditional_support = conditional_support,
      declustered_rate = declustered_rate,
      density = density,
      global_rate = global_rate,
      time_span = time_span,
      metadata = list(source = "cluster_support", gap_rule = clusters$gap_rule %||% gap_rule, eps = eps)
    )

    if (!isTRUE(store)) {
      return(out)
    }

    x$fitted$diagnostics$angular_density <- out
    return(x)
  }

  if (is.null(sample)) {
    env <- x$excursions$envelopes
    if (is.list(env) && inherits(env$sample, "envelope_sample")) {
      sample <- env$sample
    } else {
      sample <- spar_build_excursion_envelope_sample(x, store = FALSE)
    }
  }

  if (!inherits(sample, "envelope_sample")) {
    stop("`sample` must be an `envelope_sample`.", call. = FALSE)
  }

  df <- sample$data
  cl_name <- sample$cluster_name

  if (is.null(phi_grid)) {
    dom <- x$angular$domain
    n_phi_default <- if (is.null(n_phi)) sample$metadata$n_phi else n_phi
    if (!is.finite(n_phi_default) || length(n_phi_default) != 1L || n_phi_default < 2) {
      n_phi_default <- 180L
    }

    if (inherits(dom, "spar_angle_domain")) {
      phi_grid <- spar_angle_grid(n = as.integer(n_phi_default), domain = dom)
    } else {
      phi_name <- sample$angle_name
      phi_grid <- sort(unique(df[[phi_name]]))
    }
  } else {
    phi_grid <- as.numeric(phi_grid)
    if (length(phi_grid) < 2L || any(!is.finite(phi_grid))) {
      stop("`phi_grid` must be a numeric vector with at least two finite values.", call. = FALSE)
    }
  }

  density_method <- match.arg(density_method)
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0) {
    stop("`eps` must be a single non-negative numeric value.", call. = FALSE)
  }

  n_clusters_total <- length(unique(df[[cl_name]]))

  support_count <- numeric(length(phi_grid))

  if (identical(density_method, "legacy_support_rate")) {
    support_count <- vapply(phi_grid, function(phi0) {
      d0 <- extract_angle_distribution(
        sample = sample,
        phi0 = phi0,
        mode = "interpolate"
      )
      if (nrow(d0) == 0L) return(0)
      sum(is.finite(d0$excess) & (d0$excess > 0))
    }, numeric(1))

    conditional_support <- if (n_clusters_total > 0) support_count / n_clusters_total else rep(NA_real_, length(phi_grid))
    declustered_rate <- support_count / time_span

    w <- angular_grid_weights(phi_grid, normalize = TRUE)
    mass <- declustered_rate
    scale <- sum(mass * w, na.rm = TRUE)
    density <- if (is.finite(scale) && scale > 0) mass / scale else rep(NA_real_, length(mass))
  } else {
    dom <- x$angular$domain
    if (!inherits(dom, "spar_angle_domain")) {
      dom <- sample$angle_domain
    }
    if (!inherits(dom, "spar_angle_domain")) {
      dom <- new_spar_angle_domain(type = "interval", lower = min(phi_grid), upper = max(phi_grid))
    }
    dom <- as_spar_angle_domain(dom)

    phi_grid_u <- if (dom$type == "cyclical") {
      phi_n <- spar_normalize_angle(phi_grid, dom)
      start <- dom$lower
      pu <- phi_n
      pu[pu < start] <- pu[pu < start] + dom$width
      pu
    } else {
      as.numeric(phi_grid)
    }

    grid_step <- {
      d <- diff(sort(unique(phi_grid_u)))
      d <- d[is.finite(d) & d > 0]
      if (length(d) == 0L) 1e-6 else median(d)
    }

    density_sum <- numeric(length(phi_grid))

    split_df <- split(df, df[[cl_name]])

    for (cid in names(split_df)) {
      dfi <- split_df[[cid]]
      phi <- as.numeric(dfi[[sample$angle_name]])
      y <- as.numeric(dfi[[sample$response_name]])

      ok <- is.finite(phi) & is.finite(y)
      phi <- phi[ok]
      y <- y[ok]
      if (length(phi) == 0L) next

      phi <- spar_normalize_angle(phi, dom)

      if (dom$type == "cyclical") {
        uw <- unwrap_phi_for_surface(phi, angle_domain = dom)
        phi_u <- uw$phi_unwrap
        start <- uw$start
      } else {
        phi_u <- phi
        start <- min(phi_u)
      }

      ord <- order(phi_u)
      phi_u <- phi_u[ord]
      y <- y[ord]

      supp <- y > eps
      K <- sum(supp)
      if (K == 0L) next

      r <- rle(supp)
      ends <- cumsum(r$lengths)
      starts <- ends - r$lengths + 1L

      cluster_density <- numeric(length(phi_grid))
      cluster_support <- rep(FALSE, length(phi_grid))
      m_i <- 1 / K

      for (ri in seq_along(r$values)) {
        if (!isTRUE(r$values[ri])) next
        idx <- starts[ri]:ends[ri]
        p <- phi_u[idx]
        if (length(p) == 0L) next

        if (length(p) == 1L) {
          left <- p - 0.5 * grid_step
          right <- p + 0.5 * grid_step
          delta <- right - left
          if (!is.finite(delta) || delta <= 0) next
          fval <- m_i / delta

          sel <- phi_grid_u >= left & phi_grid_u <= right
          cluster_density[sel] <- fval
          cluster_support[sel] <- TRUE
          next
        }

        mids <- 0.5 * (p[-length(p)] + p[-1L])
        lefts <- c(p[1] - 0.5 * (p[2] - p[1]), mids)
        rights <- c(mids, p[length(p)] + 0.5 * (p[length(p)] - p[length(p) - 1L]))

        for (j in seq_along(p)) {
          delta <- rights[j] - lefts[j]
          if (!is.finite(delta) || delta <= 0) next
          fval <- m_i / delta

          if (j < length(p)) {
            sel <- phi_grid_u >= lefts[j] & phi_grid_u < rights[j]
          } else {
            sel <- phi_grid_u >= lefts[j] & phi_grid_u <= rights[j]
          }
          cluster_density[sel] <- fval
          cluster_support[sel] <- TRUE
        }
      }

      density_sum <- density_sum + cluster_density
      support_count <- support_count + as.numeric(cluster_support)
    }

    conditional_support <- if (n_clusters_total > 0) support_count / n_clusters_total else rep(NA_real_, length(phi_grid))
    declustered_rate <- support_count / time_span

    density_cluster_mean <- if (n_clusters_total > 0) density_sum / n_clusters_total else rep(NA_real_, length(phi_grid))
    w <- angular_grid_weights(phi_grid, normalize = TRUE)
    scale <- sum(density_cluster_mean * w, na.rm = TRUE)
    density <- if (is.finite(scale) && scale > 0) density_cluster_mean / scale else rep(NA_real_, length(phi_grid))
  }

  clusters <- x$excursions$clusters
  global_rate <- if (!is.null(clusters) && inherits(clusters, "clustered_excursion_spans")) {
    nrow(clusters$spans) / time_span
  } else {
    NA_real_
  }

  out <- new_angular_density_diagnostic(
    phi_grid = phi_grid,
    support_count = support_count,
    conditional_support = conditional_support,
    declustered_rate = declustered_rate,
    density = density,
    global_rate = global_rate,
    time_span = time_span,
    metadata = list(source = "spar_representation", density_method = density_method, eps = eps)
  )

  if (!isTRUE(store)) {
    return(out)
  }

  x$fitted$diagnostics$angular_density <- out
  x
}

#' Extract regional declustered maxima per excursion span
#'
#' Returns one maximum exceedance per declustered span/cluster within a target
#' angular interval.
#'
#' @param x A `spar_representation` object.
#' @param region_phi Numeric length-2 angular interval.
#' @param gap_rule Declustering gap used when clusters are not already stored.
#' @param closed Logical; whether angular interval bounds are closed.
#' @param use_stored Logical; if `TRUE`, reuses stored clustered spans when
#'   available.
#'
#' @return A data frame with columns `cluster_id` and `max_excess`.
#'
#' @export
spar_extract_region_declustered_maxima <- function(
    x,
    region_phi,
    gap_rule = 0,
    closed = TRUE,
    use_stored = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (!is.numeric(region_phi) || length(region_phi) != 2L || anyNA(region_phi)) {
    stop("`region_phi` must be a numeric length-2 vector.", call. = FALSE)
  }

  if (is.null(x$angular$phi)) {
    stop("`x$angular$phi` is required.", call. = FALSE)
  }

  if (is.null(x$excess$value)) {
    stop("`x$excess$value` is required.", call. = FALSE)
  }

  clusters <- if (isTRUE(use_stored)) x$excursions$clusters else NULL
  if (is.null(clusters) || !inherits(clusters, "clustered_excursion_spans")) {
    dec <- spar_decluster_excursions(x, gap_rule = gap_rule, store = FALSE)
    clusters <- dec$clusters
  }

  spans <- clusters$spans
  if (nrow(spans) == 0L) {
    return(data.frame(cluster_id = integer(0), max_excess = numeric(0)))
  }

  dom <- x$angular$domain
  if (is.null(dom)) {
    dom <- new_spar_angle_domain(type = "interval", lower = min(x$angular$phi, na.rm = TRUE), upper = max(x$angular$phi, na.rm = TRUE))
  }

  in_region <- spar_angle_in_range(as.numeric(x$angular$phi), region_phi, dom, closed = closed)
  exc <- as.numeric(x$excess$value)

  out <- lapply(seq_len(nrow(spans)), function(i) {
    idx <- seq.int(spans$First[i], spans$Last[i])
    keep <- in_region[idx] & is.finite(exc[idx]) & exc[idx] > 0
    if (!any(keep)) {
      return(NULL)
    }
    data.frame(cluster_id = spans$cluster_id[i], max_excess = max(exc[idx][keep], na.rm = TRUE))
  })

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) {
    return(data.frame(cluster_id = integer(0), max_excess = numeric(0)))
  }

  do.call(rbind, out)
}

#' Extract regional envelope maxima per cluster
#'
#' Returns one maximum envelope excess per cluster within a target angular
#' interval.
#'
#' @param x A `spar_representation` object.
#' @param region_phi Numeric length-2 angular interval.
#' @param sample Optional `envelope_sample`. If `NULL`, uses stored sample or
#'   builds one.
#' @param gap_rule Declustering gap used when building sample.
#' @param n_phi Angle-grid size when building sample.
#' @param lambda Envelope mixing weight when building sample.
#' @param n_sub Segment interpolation density when building sample.
#' @param space Data space used when building sample.
#' @param closed Logical; whether angular interval bounds are closed.
#'
#' @return A data frame with columns `cluster_id` and `max_excess`.
#'
#' @export
spar_extract_region_envelope_maxima <- function(
    x,
    region_phi,
    sample = NULL,
    gap_rule = 0,
    n_phi = 200L,
    lambda = 0,
    n_sub = 100L,
    space = c("transformed", "original"),
    closed = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (!is.numeric(region_phi) || length(region_phi) != 2L || anyNA(region_phi)) {
    stop("`region_phi` must be a numeric length-2 vector.", call. = FALSE)
  }

  space <- match.arg(space)

  if (is.null(sample)) {
    env <- x$excursions$envelopes
    if (is.list(env) && inherits(env$sample, "envelope_sample")) {
      sample <- env$sample
    } else {
      sample <- spar_build_excursion_envelope_sample(
        x,
        gap_rule = gap_rule,
        n_phi = n_phi,
        lambda = lambda,
        n_sub = n_sub,
        space = space,
        store = FALSE
      )
    }
  }

  if (!inherits(sample, "envelope_sample")) {
    stop("`sample` must be an `envelope_sample`.", call. = FALSE)
  }

  df <- sample$data
  if (nrow(df) == 0L) {
    return(data.frame(cluster_id = integer(0), max_excess = numeric(0)))
  }

  angle_name <- sample$angle_name
  response_name <- sample$response_name
  cluster_name <- sample$cluster_name

  dom <- sample$angle_domain %||% x$angular$domain
  if (is.null(dom)) {
    dom <- new_spar_angle_domain(type = "interval", lower = min(df[[angle_name]], na.rm = TRUE), upper = max(df[[angle_name]], na.rm = TRUE))
  }

  keep <- spar_angle_in_range(as.numeric(df[[angle_name]]), region_phi, dom, closed = closed) &
    is.finite(df[[response_name]]) & (df[[response_name]] > 0)

  dsub <- df[keep, c(cluster_name, response_name), drop = FALSE]
  if (nrow(dsub) == 0L) {
    return(data.frame(cluster_id = integer(0), max_excess = numeric(0)))
  }

  vals <- tapply(dsub[[response_name]], dsub[[cluster_name]], max)
  data.frame(cluster_id = as.integer(names(vals)), max_excess = as.numeric(vals))
}

new_excursion_upper_path <- function(cluster_id,
                                     idx,
                                     t,
                                     phi,
                                     excess,
                                     R = NULL,
                                     X = NULL,
                                     metadata = list()) {
  stopifnot(length(idx) == length(phi), length(phi) == length(excess), length(t) == length(phi))
  structure(
    list(
      cluster_id = cluster_id,
      idx = idx,
      t = t,
      phi = phi,
      excess = excess,
      R = R,
      X = X,
      metadata = metadata
    ),
    class = "excursion_upper_path"
  )
}

print.excursion_upper_path <- function(x, ...) {
  cat("excursion_upper_path\n")
  cat("  cluster_id:", x$cluster_id, "\n")
  cat("  points    :", length(x$phi), "\n")
  cat("  max excess:", max(x$excess, na.rm = TRUE), "\n")
  invisible(x)
}

new_excursion_upper_path_group <- function(paths,
                                           spans = NULL,
                                           metadata = list()) {
  stopifnot(is.list(paths))
  structure(
    list(
      paths = paths,
      spans = spans,
      metadata = metadata
    ),
    class = "excursion_upper_path_group"
  )
}

unwrap_phi_for_surface <- function(phi, angle_domain = NULL) {
  phi <- as.numeric(phi)
  if (is.null(angle_domain)) {
    return(list(phi_norm = phi, phi_unwrap = phi, start = NA_real_))
  }

  dom <- as_spar_angle_domain(angle_domain)
  phi_norm <- spar_normalize_angle(phi, dom)

  if (dom$type != "cyclical") {
    return(list(phi_norm = phi_norm, phi_unwrap = phi_norm, start = min(phi_norm, na.rm = TRUE)))
  }

  p <- sort(phi_norm)
  if (length(p) <= 1L) {
    return(list(phi_norm = phi_norm, phi_unwrap = phi_norm, start = p[1] %||% dom$lower))
  }

  gaps <- diff(c(p, p[1] + dom$width))
  i_gap <- which.max(gaps)
  start <- p[(i_gap %% length(p)) + 1L]

  phi_u <- phi_norm
  phi_u[phi_u < start] <- phi_u[phi_u < start] + dom$width

  list(phi_norm = phi_norm, phi_unwrap = phi_u, start = start)
}

transition_angle_position <- function(phi, p1, p2, angle_domain, tol = 1e-10) {
  dom <- as_spar_angle_domain(angle_domain)
  r <- spar_angle_range(c(p1, p2), dom)
  if (anyNA(r) || !isTRUE(spar_angle_in_range(phi, r, dom, closed = FALSE))) {
    return(NA_real_)
  }

  if (dom$type == "cyclical") {
    span <- spar_smallest_angular_span(p1, p2, dom)
    total <- span$delta
    step <- spar_angle_delta(phi, span$from, dom)
  } else {
    total <- p2 - p1
    step <- phi - p1
  }

  if (!is.finite(total) || abs(total) <= tol) {
    return(NA_real_)
  }

  alpha <- step / total
  if (!is.finite(alpha) || alpha <= tol || alpha >= 1 - tol) {
    return(NA_real_)
  }

  alpha
}

transition_dominates_angular <- function(phi_i, excess_i, phi1, phi2, y1, y2,
                                         angle_domain, tol = 1e-10) {
  alpha <- transition_angle_position(phi_i, phi1, phi2, angle_domain, tol = tol)
  if (!is.finite(alpha)) {
    return(FALSE)
  }

  y_line <- y1 + alpha * (y2 - y1)
  is.finite(y_line) && y_line >= excess_i - tol
}

transition_dominates_transformed <- function(X_i, R_i, phi_i, X1, X2, R1, R2,
                                             phi1, phi2, angle_domain, tol = 1e-10) {
  alpha <- transition_angle_position(phi_i, phi1, phi2, angle_domain, tol = tol)
  if (!is.finite(alpha)) {
    return(FALSE)
  }

  if (!is.finite(R_i) || R_i > max(R1, R2, na.rm = TRUE) + tol) {
    return(FALSE)
  }

  X_i <- as.numeric(X_i)
  X1 <- as.numeric(X1)
  X2 <- as.numeric(X2)
  if (length(X_i) != 2L || length(X1) != 2L || length(X2) != 2L ||
      any(!is.finite(c(X_i, X1, X2))) || abs(R_i) <= tol) {
    return(FALSE)
  }

  dir <- X_i / R_i
  v <- X2 - X1
  A <- matrix(c(dir, -v), nrow = 2L)
  det_A <- A[1L, 1L] * A[2L, 2L] - A[1L, 2L] * A[2L, 1L]
  if (!is.finite(det_A) || abs(det_A) <= tol) {
    return(FALSE)
  }

  sol <- solve(A, X1)
  a <- sol[1L]
  s <- sol[2L]

  is.finite(a) && is.finite(s) && s > tol && s < 1 - tol && a >= R_i - tol
}

extract_upper_excursion_path <- function(path,
                                          angle_domain = NULL,
                                          interpolation = c("transformed", "angular"),
                                          tol = 1e-10) {
  stopifnot(inherits(path, "excursion_path"))

  interpolation <- match.arg(interpolation)

  phi <- as.numeric(path$phi)
  excess <- as.numeric(path$excess)
  R <- as.numeric(path$R)
  n <- length(phi)

  if (is.null(path$t)) {
    t <- seq_len(n)
  } else {
    t <- as.numeric(path$t)
  }

  ok <- is.finite(phi) & is.finite(excess)
  if (!any(ok)) {
    return(new_excursion_upper_path(
      cluster_id = path$metadata$cluster_id %||% NA_integer_,
      idx = integer(0),
      t = numeric(0),
      phi = numeric(0),
      excess = numeric(0),
      R = numeric(0),
      X = matrix(numeric(0), nrow = 0)
    ))
  }

  idx0 <- seq_len(n)[ok]
  phi <- phi[ok]
  excess <- excess[ok]
  R <- R[ok]
  t <- t[ok]
  X <- path$X[ok, , drop = FALSE]

  if (is.null(angle_domain)) {
    phi_lower <- min(phi, na.rm = TRUE)
    phi_upper <- max(phi, na.rm = TRUE)
    if (phi_upper <= phi_lower) {
      phi_upper <- phi_lower + 1
    }
    angle_domain <- new_spar_angle_domain(type = "interval", lower = phi_lower, upper = phi_upper)
  }
  dom <- as_spar_angle_domain(angle_domain)
  phi_n <- spar_normalize_angle(phi, dom)

  if (identical(interpolation, "transformed") && ncol(X) != 2L) {
    stop("`interpolation = 'transformed'` currently requires bivariate `path$X`.", call. = FALSE)
  }

  dominated <- rep(FALSE, length(phi_n))
  m <- length(phi_n)
  if (m >= 3L) {
    valid_transition <- is.finite(t[-m]) & is.finite(t[-1L]) & abs(diff(t) - 1) <= tol &
      is.finite(phi_n[-m]) & is.finite(phi_n[-1L]) &
      is.finite(excess[-m]) & is.finite(excess[-1L]) &
      is.finite(R[-m]) & is.finite(R[-1L])

    if (identical(interpolation, "transformed")) {
      D <- X / R
      valid_direction <- is.finite(R) & abs(R) > tol & rowSums(!is.finite(D)) == 0L
      valid_transition <- valid_transition &
        rowSums(!is.finite(X[-m, , drop = FALSE])) == 0L &
        rowSums(!is.finite(X[-1L, , drop = FALSE])) == 0L
    }

    for (j in which(valid_transition)) {
      candidates <- which(!dominated)
      candidates <- candidates[candidates != j & candidates != j + 1L]
      if (length(candidates) == 0L) {
        break
      }

      tr <- spar_angle_range(phi_n[c(j, j + 1L)], dom)
      if (anyNA(tr)) {
        next
      }

      candidates <- candidates[spar_angle_in_range(phi_n[candidates], tr, dom, closed = FALSE)]
      if (length(candidates) == 0L) {
        next
      }

      if (identical(interpolation, "angular")) {
        if (dom$type == "cyclical") {
          span <- spar_smallest_angular_span(phi_n[j], phi_n[j + 1L], dom)
          total <- span$delta
          alpha <- spar_angle_delta(phi_n[candidates], span$from, dom) / total
        } else {
          total <- phi_n[j + 1L] - phi_n[j]
          alpha <- (phi_n[candidates] - phi_n[j]) / total
        }

        ok_alpha <- is.finite(alpha) & alpha > tol & alpha < 1 - tol
        if (!any(ok_alpha)) {
          next
        }

        candidates <- candidates[ok_alpha]
        alpha <- alpha[ok_alpha]
        y_line <- excess[j] + alpha * (excess[j + 1L] - excess[j])
        dominated[candidates[is.finite(y_line) & y_line >= excess[candidates] - tol]] <- TRUE
        next
      }

      candidates <- candidates[valid_direction[candidates] & R[candidates] <= max(R[j], R[j + 1L]) + tol]
      if (length(candidates) == 0L) {
        next
      }

      x1 <- X[j, 1L]; y1 <- X[j, 2L]
      x2 <- X[j + 1L, 1L]; y2 <- X[j + 1L, 2L]
      vx <- x2 - x1; vy <- y2 - y1
      if (!is.finite(vx) || !is.finite(vy) || (abs(vx) <= tol && abs(vy) <= tol)) {
        next
      }

      L1 <- y1 - y2
      L2 <- x2 - x1
      L3 <- x1 * y2 - y1 * x2
      denom <- L1 * D[candidates, 1L] + L2 * D[candidates, 2L]
      ok_denom <- is.finite(denom) & abs(denom) > tol
      if (!any(ok_denom)) {
        next
      }

      candidates <- candidates[ok_denom]
      denom <- denom[ok_denom]
      lambda <- -L3 / denom
      ok_lambda <- is.finite(lambda) & lambda >= R[candidates] - tol
      if (!any(ok_lambda)) {
        next
      }

      candidates <- candidates[ok_lambda]
      lambda <- lambda[ok_lambda]

      if (abs(vx) >= abs(vy)) {
        s <- (lambda * D[candidates, 1L] - x1) / vx
      } else {
        s <- (lambda * D[candidates, 2L] - y1) / vy
      }

      dominated[candidates[is.finite(s) & s > tol & s < 1 - tol]] <- TRUE
    }
  }

  keep <- !dominated
  idx0 <- idx0[keep]
  t <- t[keep]
  phi_n <- phi_n[keep]
  excess <- excess[keep]
  R <- R[keep]
  X <- X[keep, , drop = FALSE]

  new_excursion_upper_path(
    cluster_id = path$metadata$cluster_id %||% NA_integer_,
    idx = idx0,
    t = t,
    phi = phi_n,
    excess = excess,
    R = R,
    X = X,
    metadata = list(source = "path_transition_subset", interpolation = interpolation)
  )
}

#' Build upper-path subsets from excursion paths
#'
#' @param x A `spar_representation` object.
#' @param path_group Optional `excursion_path_group`. If `NULL`, uses stored
#'   paths or builds paths from current state.
#' @param gap_rule Declustering gap used if path group must be built.
#' @param space Data space if path group must be built.
#' @param interpolation Interpolation space used for transition domination.
#' @param store Logical; if `TRUE`, stores upper paths in
#'   `x$excursions$upper_paths` and returns updated `x`.
#'
#' @return `excursion_upper_path_group` or updated `spar_representation`.
#'
#' @export
spar_build_upper_excursion_paths <- function(
    x,
    path_group = NULL,
    gap_rule = 0,
    space = c("transformed", "original"),
    interpolation = c("transformed", "angular"),
    store = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  space <- match.arg(space)
  interpolation <- match.arg(interpolation)

  if (is.null(path_group)) {
    path_group <- x$excursions$paths
  }

  if (is.null(path_group)) {
    path_group <- spar_build_excursion_path_group(
      x,
      gap_rule = gap_rule,
      space = space,
      keep_prev_next = TRUE,
      store = FALSE
    )
  }

  if (!inherits(path_group, "excursion_path_group")) {
    stop("`path_group` must be an `excursion_path_group`.", call. = FALSE)
  }

  dom <- x$angular$domain
  uppers <- lapply(path_group$paths, function(path) {
    extract_upper_excursion_path(path, angle_domain = dom, interpolation = interpolation)
  })
  names(uppers) <- names(path_group$paths)

  out <- new_excursion_upper_path_group(
    paths = uppers,
    spans = path_group$spans,
    metadata = list(source = "spar_representation", data_space = space, interpolation = interpolation)
  )

  if (!isTRUE(store)) {
    return(out)
  }

  x$excursions$upper_paths <- out
  x
}

#' Extract regional maxima from upper excursion paths
#'
#' @param x A `spar_representation` object.
#' @param region_phi Numeric length-2 angular interval.
#' @param upper_paths Optional `excursion_upper_path_group`. If `NULL`, uses
#'   stored upper paths or builds them.
#' @param gap_rule Declustering gap used when building upper paths.
#' @param space Data space used when building upper paths.
#' @param closed Logical; whether angular interval bounds are closed.
#'
#' @return Data frame with `cluster_id` and `max_excess`.
#'
#' @export
spar_extract_region_upper_path_maxima <- function(
    x,
    region_phi,
    upper_paths = NULL,
    gap_rule = 0,
    space = c("transformed", "original"),
    closed = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  if (!is.numeric(region_phi) || length(region_phi) != 2L || anyNA(region_phi)) {
    stop("`region_phi` must be a numeric length-2 vector.", call. = FALSE)
  }

  space <- match.arg(space)

  if (is.null(upper_paths)) {
    upper_paths <- x$excursions$upper_paths
  }

  if (is.null(upper_paths) || !inherits(upper_paths, "excursion_upper_path_group")) {
    upper_paths <- spar_build_upper_excursion_paths(
      x,
      path_group = NULL,
      gap_rule = gap_rule,
      space = space,
      store = FALSE
    )
  }

  dom <- x$angular$domain
  if (is.null(dom)) {
    dom <- new_spar_angle_domain(type = "interval", lower = min(x$angular$phi, na.rm = TRUE), upper = max(x$angular$phi, na.rm = TRUE))
  }

  out <- lapply(upper_paths$paths, function(up) {
    if (length(up$phi) == 0L) return(NULL)
    keep <- spar_angle_in_range(up$phi, region_phi, dom, closed = closed) & is.finite(up$excess) & up$excess > 0
    if (!any(keep)) return(NULL)
    data.frame(cluster_id = up$cluster_id, max_excess = max(up$excess[keep], na.rm = TRUE))
  })

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) {
    return(data.frame(cluster_id = integer(0), max_excess = numeric(0)))
  }

  do.call(rbind, out)
}

#' Extract canonical GPD comparison sets from a representation
#'
#' Returns the three standard regional samples used in diagnostics:
#' all exceedance observations, declustered span maxima, and upper-path maxima.
#'
#' @param x A `spar_representation` object.
#' @param region_phi Numeric length-2 angular interval.
#' @param gap_rule Declustering gap rule used for declustered and upper-path sets.
#' @param closed Logical; whether angular interval bounds are closed.
#'
#' @return Named list with entries:
#'   `all_obs`, `declustered_span_max`, `upper_path_peak_max`.
#'
#' @export
spar_extract_gpd_comparison_sets <- function(
    x,
    region_phi,
    gap_rule = 6,
    closed = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  if (!is.numeric(region_phi) || length(region_phi) != 2L || anyNA(region_phi)) {
    stop("`region_phi` must be a numeric length-2 vector.", call. = FALSE)
  }
  if (is.null(x$angular$phi) || is.null(x$angular$domain)) {
    stop("`x$angular$phi` and `x$angular$domain` are required.", call. = FALSE)
  }
  if (is.null(x$excess$value)) {
    stop("`x$excess$value` is required.", call. = FALSE)
  }

  in_region <- spar_angle_in_range(x$angular$phi, region_phi, x$angular$domain, closed = closed)
  all_obs <- as.numeric(x$excess$value[in_region & x$excess$value > 0])

  decl_df <- spar_extract_region_declustered_maxima(
    x,
    region_phi = region_phi,
    gap_rule = gap_rule,
    closed = closed,
    use_stored = TRUE
  )

  upper_df <- spar_extract_region_upper_path_maxima(
    x,
    region_phi = region_phi,
    gap_rule = gap_rule,
    closed = closed
  )

  list(
    all_obs = as.numeric(all_obs),
    declustered_span_max = as.numeric(decl_df$max_excess),
    upper_path_peak_max = as.numeric(upper_df$max_excess)
  )
}
