# ============================================================
# Excursion-envelope implementation skeleton for R
#
# Goal:
#   1. Store a cluster path in original space
#   2. Build an envelope from interpolation in original space
#   3. Build an envelope from interpolation in angular-radial space
#   4. Mix the two envelopes with lambda in [0,1]
#        lambda = 0  -> original-space interpolation only  (default)
#        lambda = 1  -> angular-radial interpolation only
#   5. Export a weighted grid dataset suitable for evgam
#
# Notes:
#   - This is a skeleton, not a final package design.
#   - It assumes a 2D original space for plotting convenience,
#     but the original-space interpolation itself works for d >= 2
#     as long as radial_fun() and angle_fun() are defined.
#   - Angle is treated as non-circular here. If your pseudo-angle
#     can wrap, add an unwrap step before building envelopes.

# ============================================================
# 0. Utilities

`%||%` <- function(x, y) if (is.null(x)) y else x

clamp01 <- function(x) pmin(1, pmax(0, x))


# ============================================================
# 6. Envelope from sampled segment curves
#
# Given a collection of sampled curves (phi, Y), construct an
# upper envelope H on a fixed angular grid.
#
# For each subsegment between consecutive sampled points, we linearly
# interpolate within that small piece and update H(grid) wherever
# the segment spans the grid cell.

empty_envelope <- function(phi_grid) {
  structure(
    list(
      phi = phi_grid,
      H = rep(0, length(phi_grid))
    ),
    class = "excursion_envelope"
  )
}

print.excursion_envelope <- function(x, ...) {
  cat("excursion_envelope\n")
  cat("  grid size :", length(x$phi), "\n")
  cat("  max H     :", max(x$H, na.rm = TRUE), "\n")
  invisible(x)
}

split_segment_on_domain <- function(p1, p2, y1, y2, angle_domain = NULL) {
  if (is.null(angle_domain)) {
    return(list(list(p1 = p1, p2 = p2, y1 = y1, y2 = y2)))
  }

  angle_domain <- as_spar_angle_domain(angle_domain)

  if (angle_domain$type != "cyclical") {
    return(list(list(p1 = p1, p2 = p2, y1 = y1, y2 = y2)))
  }

  span <- spar_smallest_angular_span(p1, p2, angle_domain)
  seg <- span$segments

  if (nrow(seg) == 1L) {
    return(list(list(p1 = seg$start[1], p2 = seg$end[1], y1 = y1, y2 = y2)))
  }

  d_total <- span$delta
  if (!is.finite(d_total) || abs(d_total) < .Machine$double.eps^0.5) {
    return(list(list(p1 = seg$start[1], p2 = seg$start[1], y1 = y1, y2 = y2)))
  }

  t_break <- abs(seg$end[1] - seg$start[1]) / abs(d_total)
  y_break <- y1 + t_break * (y2 - y1)

  list(
    list(p1 = seg$start[1], p2 = seg$end[1], y1 = y1, y2 = y_break),
    list(p1 = seg$start[2], p2 = seg$end[2], y1 = y_break, y2 = y2)
  )
}

segment_index_range <- function(phi_grid, a, b) {
  lo <- min(a, b)
  hi <- max(a, b)

  if (hi < phi_grid[1] || lo > phi_grid[length(phi_grid)]) {
    return(integer(0))
  }

  i_lo <- findInterval(lo, phi_grid)
  if (i_lo < 1L) i_lo <- 1L
  if (phi_grid[i_lo] < lo && i_lo < length(phi_grid)) i_lo <- i_lo + 1L

  i_hi <- findInterval(hi, phi_grid)
  if (i_hi < 1L) return(integer(0))
  if (i_hi > length(phi_grid)) i_hi <- length(phi_grid)

  if (i_lo > i_hi) return(integer(0))
  i_lo:i_hi
}

update_envelope_from_polyline <- function(H, phi_grid, phi_vec, Y_vec, angle_domain = NULL) {
  stopifnot(length(phi_vec) == length(Y_vec), length(H) == length(phi_grid))

  ok <- is.finite(phi_vec) & is.finite(Y_vec)
  phi_vec <- phi_vec[ok]
  Y_vec <- Y_vec[ok]

  if (length(phi_vec) < 2L) return(H)

  for (i in seq_len(length(phi_vec) - 1L)) {
    H <- update_envelope_from_segment(
      H = H,
      phi_grid = phi_grid,
      p1 = phi_vec[i],
      p2 = phi_vec[i + 1L],
      y1 = Y_vec[i],
      y2 = Y_vec[i + 1L],
      angle_domain = angle_domain
    )
  }

  H
}

update_envelope_from_segment <- function(H, phi_grid, p1, p2, y1, y2, angle_domain = NULL) {
  if (!is.finite(p1) || !is.finite(p2) || !is.finite(y1) || !is.finite(y2)) {
    return(H)
  }

  segs <- split_segment_on_domain(p1, p2, y1, y2, angle_domain = angle_domain)

  for (seg in segs) {
    idx <- segment_index_range(phi_grid, seg$p1, seg$p2)
    if (length(idx) == 0L) next

    if (abs(seg$p2 - seg$p1) < .Machine$double.eps^0.5) {
      y_line <- pmax(seg$y1, seg$y2)
      H[idx] <- pmax(H[idx], y_line)
    } else {
      t <- (phi_grid[idx] - seg$p1) / (seg$p2 - seg$p1)
      y_line <- seg$y1 + t * (seg$y2 - seg$y1)
      H[idx] <- pmax(H[idx], y_line)
    }
  }

  H
}

build_segment_original_space_vectorized <- function(Xa, Xb, u, radial_fun, angle_fun,
                                                    n_sub = 100L) {
  s <- seq(0, 1, length.out = n_sub)
  Xa <- as.numeric(Xa)
  Xb <- as.numeric(Xb)

  Xs <- tcrossprod(1 - s, Xa) + tcrossprod(s, Xb)

  u_s <- (1 - s) * u[1] + s * u[2]
  tr <- transform_points(Xs, radial_fun, angle_fun, u_s)

  list(phi = tr$phi, Y = tr$Y)
}

# ============================================================
# 7. Build original-space envelope
#
# Default closure:
#   if no neighboring sub-threshold points are supplied,
#   use vertical closure at the start/end path angles.

envelope_original_space_bruteforce <- function(path,
                                                phi_grid,
                                                X_prev = NULL,
                                                X_next = NULL,
                                                angle_domain = NULL,
                                                n_sub = 100L) {
  stopifnot(inherits(path, "excursion_path"))

  X <- path$X
  u <- path$u
  phi <- path$phi
  R <- path$R
  excess <- path$excess
  radial_fun <- path$radial_fun
  angle_fun <- path$angle_fun


  H <- rep(0, length(phi_grid))

  # infer start/end closure angles if possible
  clos <- infer_closure_original(
    X_prev = X_prev,
    X1 = X[1, ],
    X_last = X[nrow(X), ],
    X_next = X_next,
    u = u,
    radial_fun = radial_fun,
    angle_fun = angle_fun
  )

  phi0   <- clos$phi0   %||% phi[1]
  phiend <- clos$phi_end %||% phi[length(phi)]

  # start closure triangle edge in angle-excess plane
  H <- update_envelope_from_polyline(
    H, phi_grid,
    phi_vec = c(phi0, phi[1]),
    Y_vec   = c(0,    excess[1]),
    angle_domain = angle_domain
  )

  # interior segments via original-space interpolation
  if (nrow(X) >= 2L) {
    for (i in seq_len(nrow(X) - 1L)) {
      seg <- build_segment_original_space(
        Xa = X[i, ],
        Xb = X[i + 1L, ],
        u = u[i:(i+1)],
        radial_fun = radial_fun,
        angle_fun = angle_fun,
        n_sub = n_sub
      )
      H <- update_envelope_from_polyline(H, phi_grid, seg$phi, seg$Y, angle_domain = angle_domain)
    }
  }

  # end closure triangle edge
  H <- update_envelope_from_polyline(
    H, phi_grid,
    phi_vec = c(phi[length(phi)], phiend),
    Y_vec   = c(excess[length(excess)],     0),
    angle_domain = angle_domain
  )

  structure(
    list(phi = phi_grid, H = H, mode = "original_space"),
    class = "excursion_envelope"
  )
}

envelope_original_space_optimized <- function(path,
                                               phi_grid,
                                               X_prev = NULL,
                                               X_next = NULL,
                                               angle_domain = NULL,
                                               n_sub = 100L) {
  stopifnot(inherits(path, "excursion_path"))

  X <- path$X
  u <- path$u
  phi <- path$phi
  excess <- path$excess
  radial_fun <- path$radial_fun
  angle_fun <- path$angle_fun

  H <- rep(0, length(phi_grid))

  clos <- infer_closure_original(
    X_prev = X_prev,
    X1 = X[1, ],
    X_last = X[nrow(X), ],
    X_next = X_next,
    u = u,
    radial_fun = radial_fun,
    angle_fun = angle_fun
  )

  phi0 <- clos$phi0 %||% phi[1]
  phiend <- clos$phi_end %||% phi[length(phi)]

  H <- update_envelope_from_segment(
    H = H,
    phi_grid = phi_grid,
    p1 = phi0,
    p2 = phi[1],
    y1 = 0,
    y2 = excess[1],
    angle_domain = angle_domain
  )

  if (nrow(X) >= 2L) {
    for (i in seq_len(nrow(X) - 1L)) {
      seg <- build_segment_original_space_vectorized(
        Xa = X[i, ],
        Xb = X[i + 1L, ],
        u = u[i:(i + 1L)],
        radial_fun = radial_fun,
        angle_fun = angle_fun,
        n_sub = n_sub
      )
      H <- update_envelope_from_polyline(H, phi_grid, seg$phi, seg$Y, angle_domain = angle_domain)
    }
  }

  H <- update_envelope_from_segment(
    H = H,
    phi_grid = phi_grid,
    p1 = phi[length(phi)],
    p2 = phiend,
    y1 = excess[length(excess)],
    y2 = 0,
    angle_domain = angle_domain
  )

  structure(
    list(phi = phi_grid, H = H, mode = "original_space"),
    class = "excursion_envelope"
  )
}

envelope_original_space <- envelope_original_space_optimized

# ============================================================
# 8. Build angular-radial envelope
#
# Closure uses straight lines in angle-excess plane.
# If phi0 / phi_end not supplied, default to vertical closure.

envelope_angular_radial_bruteforce <- function(path,
                                               phi_grid,
                                               phi0 = NULL,
                                               phi_end = NULL,
                                               angle_domain = NULL,
                                               n_sub = 100L) {
  stopifnot(inherits(path, "excursion_path"))

  phi <- path$phi
  excess   <- path$excess

  H <- rep(0, length(phi_grid))

  phi0   <- phi0   %||% phi[1]
  phi_end <- phi_end %||% phi[length(phi)]

  # start closure
  H <- update_envelope_from_polyline(
    H, phi_grid,
    phi_vec = c(phi0, phi[1]),
    Y_vec   = c(0,    excess[1]),
    angle_domain = angle_domain
  )

  # interior linear segments in (phi, Y)
  if (length(phi) >= 2L) {
    for (i in seq_len(length(phi) - 1L)) {
      seg <- build_segment_angular_radial(
        phi_a = phi[i],     Y_a = excess[i],
        phi_b = phi[i + 1], Y_b = excess[i + 1],
        n_sub = n_sub
      )
      H <- update_envelope_from_polyline(H, phi_grid, seg$phi, seg$Y, angle_domain = angle_domain)
    }
  }

  # end closure
  H <- update_envelope_from_polyline(
    H, phi_grid,
    phi_vec = c(phi[length(phi)], phi_end),
    Y_vec   = c(excess[length(excess)],     0),
    angle_domain = angle_domain
  )

  structure(
    list(phi = phi_grid, H = H, mode = "angular_radial"),
    class = "excursion_envelope"
  )
}

envelope_angular_radial_optimized <- function(path,
                                              phi_grid,
                                              phi0 = NULL,
                                              phi_end = NULL,
                                              angle_domain = NULL,
                                              n_sub = 100L) {
  stopifnot(inherits(path, "excursion_path"))

  phi <- path$phi
  excess <- path$excess

  H <- rep(0, length(phi_grid))

  phi0 <- phi0 %||% phi[1]
  phi_end <- phi_end %||% phi[length(phi)]

  H <- update_envelope_from_segment(
    H = H,
    phi_grid = phi_grid,
    p1 = phi0,
    p2 = phi[1],
    y1 = 0,
    y2 = excess[1],
    angle_domain = angle_domain
  )

  if (length(phi) >= 2L) {
    for (i in seq_len(length(phi) - 1L)) {
      H <- update_envelope_from_segment(
        H = H,
        phi_grid = phi_grid,
        p1 = phi[i],
        p2 = phi[i + 1L],
        y1 = excess[i],
        y2 = excess[i + 1L],
        angle_domain = angle_domain
      )
    }
  }

  H <- update_envelope_from_segment(
    H = H,
    phi_grid = phi_grid,
    p1 = phi[length(phi)],
    p2 = phi_end,
    y1 = excess[length(excess)],
    y2 = 0,
    angle_domain = angle_domain
  )

  structure(
    list(phi = phi_grid, H = H, mode = "angular_radial"),
    class = "excursion_envelope"
  )
}

envelope_angular_radial <- envelope_angular_radial_optimized

# ============================================================
# 9. Mixed envelope
# lambda = 0 -> original space only
# lambda = 1 -> angular-radial only

mix_envelopes <- function(env_original, env_ar, lambda = 0) {
  stopifnot(inherits(env_original, "excursion_envelope"))
  stopifnot(inherits(env_ar, "excursion_envelope"))
  stopifnot(length(env_original$phi) == length(env_ar$phi))
  stopifnot(max(abs(env_original$phi - env_ar$phi)) < 1e-12)

  lambda <- clamp01(lambda)

  structure(
    list(
      phi = env_original$phi,
      H = (1 - lambda) * env_original$H + lambda * env_ar$H,
      mode = "mixed",
      lambda = lambda
    ),
    class = "excursion_envelope"
  )
}

# ============================================================
# 10. Plot methods

plot.excursion_envelope <- function(x, ..., main = NULL, xlab = "phi", ylab = "H(phi)") {
  main <- main %||% paste("Envelope:", x$mode)
  plot(x$phi, x$H, type = "l", lwd = 2, main = main, xlab = xlab, ylab = ylab, ...)
}

plot_excursion_path_2d <- function(path, main = "Path in original space", ...) {
  stopifnot(inherits(path, "excursion_path"))
  stopifnot(ncol(path$X) == 2)

  plot(path$X[, 1], path$X[, 2], type = "b", pch = 19,
       xlab = "x1", ylab = "x2", main = main, ...)
}

plot_path_in_angle_excess <- function(path, main = "Path in angle-excess space", ...) {
  stopifnot(inherits(path, "excursion_path"))
  plot(path$phi, path$excess, type = "b", pch = 19,
       xlab = "phi", ylab = "excess", main = main, ...)
}

# ============================================================
# 11. Quadrature weights on a fixed angular grid
#
# These are the weights you would pass to evgam once cluster
# envelopes are stacked into a dataset.
#
# If the grid is uniform, these are all equal after normalization.

angular_grid_weights <- function(phi_grid, normalize = TRUE) {
  m <- length(phi_grid)
  stopifnot(m >= 1)

  if (m == 1L) return(1)

  w <- numeric(m)
  w[1] <- (phi_grid[2] - phi_grid[1]) / 2
  w[m] <- (phi_grid[m] - phi_grid[m - 1]) / 2

  if (m >= 3L) {
    for (j in 2:(m - 1L)) {
      w[j] <- (phi_grid[j + 1] - phi_grid[j - 1]) / 2
    }
  }

  if (normalize) w <- w / sum(w)
  w
}

# ============================================================
# 12. Export one envelope to an evgam-style dataset
#
# Each cluster contributes one normalized angular quadrature.
# Weight sum per cluster = cluster_weight (default 1).

envelope_to_dataset <- function(env,
                                cluster_id = 1L,
                                cluster_weight = 1) {
  stopifnot(inherits(env, "excursion_envelope"))

  w <- angular_grid_weights(env$phi, normalize = TRUE) * cluster_weight

  data.frame(
    cluster_id = cluster_id,
    phi = env$phi,
    excess = env$H,
    weight = w
  )
}

# ============================================================
# 13. Cluster-level convenience wrapper
#
# Builds:
#   - original-space envelope
#   - angular-radial envelope
#   - mixed envelope
# Returns all three for inspection.

build_envelope_group <- function(path_group,
                                 phi_grid = NULL,
                                 lambda = 0,
                                 angle_domain = NULL,
                                 n_phi = 180L,
                                 n_sub = 100L,
                                 method = c("optimized", "brute_force")) {
  stopifnot(inherits(path_group, "excursion_path_group"))

  library(Fragman)

  envs <- lapply_pb(path_group$paths, function(path) {
    build_cluster_envelopes(
      path = path,
      phi_grid = phi_grid,
      lambda = lambda,
      angle_domain = angle_domain,
      X_prev = path$metadata$X_prev,
      X_next = path$metadata$X_next,
      n_phi = n_phi,
      n_sub = n_sub,
      method = method
    )
  })

  names(envs) <- names(path_group$paths)
  envs
}

build_cluster_envelopes <- function(path,
                                    phi_grid = NULL,
                                    lambda = 0,
                                    X_prev = NULL,
                                    X_next = NULL,
                                    phi0 = NULL,
                                    phi_end = NULL,
                                    angle_domain = NULL,
                                    n_phi = 180L,
                                    n_sub = 100L,
                                    method = c("optimized", "brute_force")) {
  method <- match.arg(method)
  lambda <- clamp01(lambda)

  if (is.null(phi_grid)) {
    phi_grid <- infer_envelope_phi_grid(
      path = path,
      angle_domain = angle_domain,
      X_prev = X_prev,
      X_next = X_next,
      n_phi = n_phi
    )
  }

  build_x <- if (identical(method, "optimized")) envelope_original_space_optimized else envelope_original_space_bruteforce
  build_ar <- if (identical(method, "optimized")) envelope_angular_radial_optimized else envelope_angular_radial_bruteforce

  need_original <- lambda < 1
  need_angular <- lambda > 0

  env_x <- NULL
  env_ar <- NULL

  if (isTRUE(need_original)) {
    env_x <- build_x(
      path = path,
      phi_grid = phi_grid,
      X_prev = X_prev,
      X_next = X_next,
      angle_domain = angle_domain,
      n_sub = n_sub
    )
  }

  if (isTRUE(need_angular)) {
    env_ar <- build_ar(
      path = path,
      phi_grid = phi_grid,
      phi0 = phi0,
      phi_end = phi_end,
      angle_domain = angle_domain,
      n_sub = n_sub
    )
  }

  env_mix <- if (!is.null(env_x) && !is.null(env_ar)) {
    mix_envelopes(env_x, env_ar, lambda = lambda)
  } else if (!is.null(env_x)) {
    env_x
  } else {
    env_ar
  }

  list(
    original = env_x,
    angular_radial = env_ar,
    mixed = env_mix,
    method = method,
    lambda = lambda,
    phi_grid = phi_grid
  )
}

infer_envelope_phi_grid <- function(path,
                                    angle_domain = NULL,
                                    X_prev = NULL,
                                    X_next = NULL,
                                    n_phi = 180L) {
  stopifnot(inherits(path, "excursion_path"))
  n_phi <- as.integer(n_phi)
  if (!is.finite(n_phi) || n_phi < 3L) n_phi <- 180L

  phi_base <- as.numeric(path$phi)
  phi_base <- phi_base[is.finite(phi_base)]

  if (length(phi_base) == 0L) {
    if (!is.null(angle_domain)) {
      return(spar_angle_grid(n = n_phi, domain = as_spar_angle_domain(angle_domain)))
    }
    stop("Cannot infer envelope phi-grid: path has no finite angles.", call. = FALSE)
  }

  if (is.null(angle_domain)) {
    lo <- min(phi_base)
    hi <- max(phi_base)
    if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
      eps <- 1e-6
      lo <- lo - eps
      hi <- hi + eps
    }
    return(seq(lo, hi, length.out = n_phi))
  }

  dom <- as_spar_angle_domain(angle_domain)
  phi_base <- spar_normalize_angle(phi_base, dom)

  clos <- infer_closure_original(
    X_prev = X_prev,
    X1 = path$X[1, ],
    X_last = path$X[nrow(path$X), ],
    X_next = X_next,
    u = path$u,
    radial_fun = path$radial_fun,
    angle_fun = path$angle_fun
  )
  phi_closure <- c(clos$phi0, clos$phi_end)
  phi_closure <- phi_closure[is.finite(phi_closure)]

  phi_all <- c(phi_base, phi_closure)
  phi_all <- spar_normalize_angle(phi_all, dom)

  if (dom$type != "cyclical") {
    lo <- min(phi_all)
    hi <- max(phi_all)
    if (!is.finite(lo) || !is.finite(hi) || lo == hi) {
      eps <- max(1e-6, 1e-6 * abs(lo))
      lo <- lo - eps
      hi <- hi + eps
    }
    return(seq(lo, hi, length.out = n_phi))
  }

  width <- dom$upper - dom$lower
  x <- sort(unique(as.numeric(phi_all)))
  if (length(x) == 1L) {
    eps <- width / max(n_phi, 1000L)
    lo <- x[1] - eps
    hi <- x[1] + eps
    return(spar_normalize_angle(seq(lo, hi, length.out = n_phi), dom))
  }

  gaps <- c(diff(x), (x[1] + width) - x[length(x)])
  ig <- which.max(gaps)
  start <- if (ig < length(x)) x[ig + 1L] else x[1]
  span <- width - gaps[ig]

  if (!is.finite(span) || span <= 0) {
    return(spar_angle_grid(n = n_phi, domain = dom))
  }

  phi_unwrapped <- seq(start, start + span, length.out = n_phi)
  spar_normalize_angle(phi_unwrapped, dom)
}

test_excursion_envelope <- function() {
  # 14. Example data

  set.seed(1)

  # Example 2D path, assumed already to be a connected excursion path.
  # In practice this comes from your cluster/path extraction logic.
  X_path <- rbind(
    c(3.2, 1.1),
    c(4.5, 1.7),
    c(5.0, 2.6),
    c(4.2, 2.4),
    c(3.8, 1.8)
  )

  u <- 3.5

  path <- new_excursion_path(
    X = X_path,
    u = u
  )

  print(path)
  phi_grid <- seq(min(path$phi, na.rm = TRUE) - 0.05,
                  max(path$phi, na.rm = TRUE) + 0.05,
                  length.out = 150)

  envs <- build_cluster_envelopes(
    path = path,
    phi_grid = phi_grid,
    lambda = 0,   # default: original-space interpolation only
    n_sub = 200L
  )

  # 15. Inspect

  par(mfrow = c(1, 3))
  plot_excursion_path_2d(path, main = "Original space")
  plot_path_in_angle_excess(path, main = "Angle-excess points")
  plot(envs$mixed, main = "Mixed envelope (lambda = 0)")

  par(mfrow = c(1, 1))
  plot(envs$original, main = "Original-space envelope")
  lines(envs$angular_radial$phi, envs$angular_radial$H, lty = 2)
  legend("topright",
         legend = c("original-space", "angular-radial"),
         lty = c(1, 2), bty = "n")

  # ============================================================
  # 16. Build evgam-style dataset for one cluster

  cluster_df <- envelope_to_dataset(envs$mixed, cluster_id = 1L, cluster_weight = 1)
  head(cluster_df)
  sum(cluster_df$weight)

  # ============================================================
  # 17. Stack multiple clusters
  #
  # For illustration, create a second path and stack.
  # In practice you'd loop over your extracted clusters.

  X_path2 <- rbind(
    c(2.9, 1.4),
    c(3.6, 2.1),
    c(4.0, 2.8),
    c(3.4, 2.4)
  )

  path2 <- new_excursion_path(X_path2, u = u)

  envs2 <- build_cluster_envelopes(
    path = path2,
    phi_grid = phi_grid,
    lambda = 0,
    n_sub = 200L
  )

  df_all <- rbind(
    envelope_to_dataset(envs$mixed,  cluster_id = 1L, cluster_weight = 1),
    envelope_to_dataset(envs2$mixed, cluster_id = 2L, cluster_weight = 1)
  )

  aggregate(weight ~ cluster_id, data = df_all, sum)

  # ============================================================
  # 27. Example usage with the earlier envelope skeleton

  # Suppose df_all came from stacking cluster envelopes:
  #   df_all <- rbind(
  #     envelope_to_dataset(envs$mixed,  cluster_id = 1L, cluster_weight = 1),
  #     envelope_to_dataset(envs2$mixed, cluster_id = 2L, cluster_weight = 1)
  #   )

  sample_obj <- new_envelope_sample(df_all)
  print(sample_obj)

  # A single-angle diagnostic:
  diag1 <- diagnose_angle_distribution(
    sample = sample_obj,
    phi0 = median(df_all$phi),
    mode = "interpolate"   # recommended: one interpolated value per cluster
  )

  print(diag1)
  plot(diag1)

  # Window-based alternative:
  diag2 <- diagnose_angle_distribution(
    sample = sample_obj,
    phi0 = median(df_all$phi),
    bandwidth = 0.02,
    mode = "window"
  )

  print(diag2)
  plot(diag2)

  # Several angles at once:
  multi <- diagnose_multiple_angles(
    sample = sample_obj,
    phi_vec = seq(min(df_all$phi), max(df_all$phi), length.out = 5),
    mode = "interpolate"
  )

  print(multi$summary)
}



# 18. Skeleton evgam call
# This requires the evgam package installed.
# The exact family/formula may need adjustment depending on how
# you want to interpret the envelope values.
#
# A natural first pass is to treat H(phi) as an excess-like response.
#
# library(evgam)
#
# fit <- evgam(
#   excess ~ s(phi, k = 20),
#   data = df_all,
#   family = "gpd",
#   weights = weight
# )
#
# summary(fit)
# plot(fit)
#
# 19. Suggested extension points
# - angle unwrapping for circular pseudo-angle coordinates
# - non-linear interpolation in original space
# - exact threshold-crossing closure using adjacent sub-threshold points
# - envelope normalization:
#       H_norm = H / max(H)
# - richer metadata on each cluster
# - bootstrap by cluster_id

# ============================================================
# 20. Stacked envelope data structure
#
# Stores a collection of cluster-level envelope datasets.
# Each row typically corresponds to one angular grid point from
# one cluster envelope.

new_envelope_sample <- function(df,
                                angle_name = "phi",
                                angle_domain = new_spar_angle_domain(
                                  type = "cyclical",
                                  lower = -2,
                                  upper = 2,
                                  label = "default pseudo-angle"
                                ),
                                response_name = "excess",
                                cluster_name = "cluster_id",
                                weight_name = "weight",
                                metadata = list()) {
  stopifnot(is.data.frame(df))
  required <- c(angle_name, response_name, cluster_name, weight_name)
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  structure(
    list(
      data = df,
      angle_name = angle_name,
      angle_domain = angle_domain,
      response_name = response_name,
      cluster_name = cluster_name,
      weight_name = weight_name,
      metadata = metadata
    ),
    class = "envelope_sample"
  )
}

print.envelope_sample <- function(x, ...) {
  df <- x$data
  angle_name <- x$angle_name
  response_name <- x$response_name
  cluster_name <- x$cluster_name
  weight_name <- x$weight_name

  cat("envelope_sample\n")
  cat("  rows            :", nrow(df), "\n")
  cat("  clusters        :", length(unique(df[[cluster_name]])), "\n")
  cat("  angle range     :", paste0(
    signif(min(df[[angle_name]], na.rm = TRUE), 4), " to ",
    signif(max(df[[angle_name]], na.rm = TRUE), 4)
  ), "\n")
  cat("  response range  :", paste0(
    signif(min(df[[response_name]], na.rm = TRUE), 4), " to ",
    signif(max(df[[response_name]], na.rm = TRUE), 4)
  ), "\n")
  cat("  total weight    :", signif(sum(df[[weight_name]], na.rm = TRUE), 6), "\n")
  invisible(x)
}



# ============================================================
# 35. GP reference builders

gp_survival <- function(y, sigma, xi) {
  if (!is.finite(y) || !is.finite(sigma) || !is.finite(xi) || sigma <= 0 || y < 0) {
    return(NA_real_)
  }

  if (abs(xi) < 1e-10) {
    return(exp(-y / sigma))
  }

  val <- 1 + xi * y / sigma
  if (val <= 0) return(0)
  val^(-1 / xi)
}

gp_density <- function(y, sigma, xi) {
  if (!is.finite(y) || !is.finite(sigma) || !is.finite(xi) || sigma <= 0 || y < 0) {
    return(NA_real_)
  }

  if (abs(xi) < 1e-10) {
    return((1 / sigma) * exp(-y / sigma))
  }

  val <- 1 + xi * y / sigma
  if (val <= 0) return(0)
  (1 / sigma) * val^(-1 / xi - 1)
}

make_gp_reference <- function(scale_fun, shape_fun) {
  list(
    survival = function(phi, y) gp_survival(y, sigma = scale_fun(phi), xi = shape_fun(phi)),
    density  = function(phi, y) gp_density(y, sigma = scale_fun(phi), xi = shape_fun(phi))
  )
}



# ============================================================
# 42. GP simulation helpers
#
# Draw excess Y ~ GP(sigma, xi)

rgp <- function(n, sigma, xi) {
  stopifnot(length(sigma) == 1L, length(xi) == 1L, sigma > 0, n >= 0)

  u <- runif(n)

  if (abs(xi) < 1e-10) {
    return(-sigma * log(1 - u))
  }

  sigma / xi * ((1 - u)^(-xi) - 1)
}

# Vectorized version for angle-varying parameters
rgp_vec <- function(sigma, xi) {
  stopifnot(length(sigma) == length(xi))
  n <- length(sigma)

  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- rgp(1, sigma = sigma[i], xi = xi[i])
  }
  out
}

# ============================================================
# 43. Simulate one synthetic event-level sample from a GP reference
#
# This produces one value per (cluster, angle-grid point), which is
# often the cleanest comparison object for your diagnostic tooling.
#
# Inputs:
#   phi_grid      : angular grid
#   n_clusters    : number of synthetic events
#   scale_fun     : sigma(phi)
#   shape_fun     : xi(phi)
#   weight_mode   : "equal" or "quadrature"
#
# Output:
#   data.frame with columns cluster_id, phi, excess, weight

simulate_gp_reference_sample <- function(phi_grid,
                                         n_clusters,
                                         scale_fun,
                                         shape_fun,
                                         weight_mode = c("quadrature", "equal")) {
  weight_mode <- match.arg(weight_mode)

  stopifnot(is.numeric(phi_grid), length(phi_grid) >= 1L, n_clusters >= 1L)

  if (weight_mode == "quadrature") {
    w_grid <- angular_grid_weights(phi_grid, normalize = TRUE)
  } else {
    w_grid <- rep(1 / length(phi_grid), length(phi_grid))
  }

  rows <- vector("list", n_clusters)

  for (cid in seq_len(n_clusters)) {
    sigma <- vapply(phi_grid, scale_fun, numeric(1))
    xi    <- vapply(phi_grid, shape_fun, numeric(1))

    y <- rgp_vec(sigma = sigma, xi = xi)

    rows[[cid]] <- data.frame(
      cluster_id = cid,
      phi = phi_grid,
      excess = y,
      weight = w_grid
    )
  }

  do.call(rbind, rows)
}

# ============================================================
# 44. Simulate from a GP reference object created earlier
#
# Works with:
#   ref <- make_gp_reference(scale_fun = ..., shape_fun = ...)
#
# But because make_gp_reference only stored survival/density, we also
# pass scale_fun and shape_fun explicitly here.

simulate_gp_reference_sample_from_funs <- function(phi_grid,
                                                   n_clusters,
                                                   scale_fun,
                                                   shape_fun,
                                                   metadata = list()) {
  df <- simulate_gp_reference_sample(
    phi_grid = phi_grid,
    n_clusters = n_clusters,
    scale_fun = scale_fun,
    shape_fun = shape_fun,
    weight_mode = "quadrature"
  )

  new_envelope_sample(
    df = df,
    metadata = c(
      list(source = "simulated_gp_reference",
           n_clusters = n_clusters),
      metadata
    )
  )
}
