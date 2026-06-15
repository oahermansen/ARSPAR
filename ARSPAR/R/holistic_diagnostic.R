# ============================================================
# 28. Holistic calibration map object

new_holistic_diagnostic <- function(phi_grid,
                                    y_grid,
                                    value_matrix,
                                    mode = c("log_survival_ratio",
                                             "survival_difference",
                                             "log_density_ratio"),
                                    metadata = list()) {
  mode <- match.arg(mode)

  stopifnot(is.numeric(phi_grid), length(phi_grid) >= 1)
  stopifnot(is.numeric(y_grid), length(y_grid) >= 1)
  stopifnot(is.matrix(value_matrix))
  stopifnot(nrow(value_matrix) == length(phi_grid))
  stopifnot(ncol(value_matrix) == length(y_grid))

  structure(
    list(
      phi_grid = phi_grid,
      y_grid = y_grid,
      Z = value_matrix,
      mode = mode,
      metadata = metadata
    ),
    class = "holistic_diagnostic"
  )
}

print.holistic_diagnostic <- function(x, ...) {
  z <- x$Z
  cat("holistic_diagnostic\n")
  cat("  mode         :", x$mode, "\n")
  cat("  phi grid     :", length(x$phi_grid), "\n")
  cat("  excess grid  :", length(x$y_grid), "\n")
  cat("  min / max    :", signif(min(z, na.rm = TRUE), 5), "/",
      signif(max(z, na.rm = TRUE), 5), "\n")
  invisible(x)
}
# ============================================================
# 29. Empirical survival helper at one angle
#
# Uses the extracted event-level distribution at phi0.

empirical_survival_at_angle <- function(sample,
                                        phi0,
                                        y_grid,
                                        mode = c("interpolate", "window"),
                                        bandwidth = NULL) {
  mode <- match.arg(mode)

  df <- extract_angle_distribution(
    sample = sample,
    phi0 = phi0,
    bandwidth = bandwidth,
    mode = mode
  )

  if (nrow(df) == 0L) {
    return(rep(NA_real_, length(y_grid)))
  }

  y <- df$excess
  w <- df$weight

  ok <- is.finite(y) & is.finite(w) & w > 0
  y <- y[ok]
  w <- w[ok]

  if (length(y) == 0L) {
    return(rep(NA_real_, length(y_grid)))
  }

  total_w <- sum(w)

  vapply(y_grid, function(y0) {
    sum(w[y > y0]) / total_w
  }, numeric(1))
}
# ============================================================
# 30. Empirical density helper at one angle
#
# Weighted kernel density on a fixed y-grid.

empirical_density_at_angle <- function(sample,
                                       phi0,
                                       y_grid,
                                       mode = c("interpolate", "window"),
                                       bandwidth = NULL,
                                       density_bw = "nrd0") {
  mode <- match.arg(mode)

  df <- extract_angle_distribution(
    sample = sample,
    phi0 = phi0,
    bandwidth = bandwidth,
    mode = mode
  )

  if (nrow(df) < 2L) {
    return(rep(NA_real_, length(y_grid)))
  }

  y <- df$excess
  w <- df$weight

  ok <- is.finite(y) & is.finite(w) & w > 0
  y <- y[ok]
  w <- w[ok]

  if (length(y) < 2L) {
    return(rep(NA_real_, length(y_grid)))
  }

  # Base R density supports weights in recent R versions
  dens <- density(
    x = y,
    weights = w / sum(w),
    bw = density_bw,
    from = min(y_grid),
    to = max(y_grid),
    n = length(y_grid)
  )

  approx(dens$x, dens$y, xout = y_grid, rule = 2)$y
}

# ============================================================
# 31. Build holistic discrepancy map
#
# ref_survival(phi, y) or ref_density(phi, y) must be supplied
# depending on mode.

build_holistic_diagnostic <- function(sample,
                                      phi_grid,
                                      y_grid,
                                      mode = c("log_survival_ratio",
                                               "survival_difference",
                                               "log_density_ratio"),
                                      extraction_mode = c("interpolate", "window"),
                                      extraction_bandwidth = NULL,
                                      ref_survival = NULL,
                                      ref_density = NULL,
                                      eps = 1e-8,
                                      density_bw = "nrd0") {
  mode <- match.arg(mode)
  extraction_mode <- match.arg(extraction_mode)

  Z <- matrix(NA_real_, nrow = length(phi_grid), ncol = length(y_grid))

  for (j in seq_along(phi_grid)) {
    phi0 <- phi_grid[j]

    if (mode %in% c("log_survival_ratio", "survival_difference")) {
      S_emp <- empirical_survival_at_angle(
        sample = sample,
        phi0 = phi0,
        y_grid = y_grid,
        mode = extraction_mode,
        bandwidth = extraction_bandwidth
      )

      if (is.null(ref_survival)) {
        stop("ref_survival must be supplied for survival-based modes.")
      }

      S_ref <- vapply(y_grid, function(y0) ref_survival(phi0, y0), numeric(1))

      Z[j,] <- switch(
        mode,
        log_survival_ratio = log((S_emp + eps) / (S_ref + eps)),
        survival_difference = S_emp - S_ref
      )
    }

    if (mode == "log_density_ratio") {
      f_emp <- empirical_density_at_angle(
        sample = sample,
        phi0 = phi0,
        y_grid = y_grid,
        mode = extraction_mode,
        bandwidth = extraction_bandwidth,
        density_bw = density_bw
      )

      if (is.null(ref_density)) {
        stop("ref_density must be supplied for density mode.")
      }

      f_ref <- vapply(y_grid, function(y0) ref_density(phi0, y0), numeric(1))
      Z[j,] <- log((f_emp + eps) / (f_ref + eps))
    }
  }

  new_holistic_diagnostic(
    phi_grid = phi_grid,
    y_grid = y_grid,
    value_matrix = Z,
    mode = mode,
    metadata = list(
      extraction_mode = extraction_mode,
      extraction_bandwidth = extraction_bandwidth
    )
  )
}
# ============================================================
# 32. Diverging palette

diagnostic_palette <- function(n = 201, center = c("white", "gray", "green")) {
  center <- match.arg(center)

  if (center == "white") {
    return(colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(n))
  }
  if (center == "gray") {
    return(colorRampPalette(c("#2166AC", "#DDDDDD", "#B2182B"))(n))
  }
  colorRampPalette(c("#2166AC", "#1A9850", "#B2182B"))(n)
}
# ============================================================
# 33. Square heatmap plot

plot.holistic_diagnostic <- function(x,
                                     style = c("square", "circular"),
                                     center = c("white", "gray", "green"),
                                     zlim = NULL,
                                     main = NULL,
                                     xlab = "Angle",
                                     ylab = "Excess",
                                     draw_contours = FALSE,
                                     ...) {
  style <- match.arg(style)
  center <- match.arg(center)

  if (style == "circular") {
    plot_holistic_circular(x, center = center, zlim = zlim, main = main, ...)
    return(invisible(x))
  }

  Z <- x$Z
  phi_grid <- x$phi_grid
  y_grid <- x$y_grid

  if (is.null(zlim)) {
    zmax <- max(abs(Z), na.rm = TRUE)
    zlim <- c(-zmax, zmax)
  }

  pal <- diagnostic_palette(201, center = center)
  main <- main %||% paste("Holistic diagnostic:", x$mode)

  image(
    x = phi_grid,
    y = y_grid,
    z = Z,
    col = pal,
    zlim = zlim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )

  abline(v = pretty(phi_grid), h = pretty(y_grid), col = "#00000010", lty = 1)

  if (draw_contours) {
    contour(
      x = phi_grid,
      y = y_grid,
      z = Z,
      levels = pretty(zlim, n = 8),
      add = TRUE,
      drawlabels = FALSE
    )
  }

  invisible(x)
}
# ============================================================
# 34. Circular plot
# Tile centers are mapped to polar coordinates.
# This is a visual approximation, but often very informative.

plot_holistic_circular <- function(x,
                                   center = c("white", "gray", "green"),
                                   zlim = NULL,
                                   main = NULL,
                                   point_cex = 1.2,
                                   ...) {
  center <- match.arg(center)

  Z <- x$Z
  phi_grid <- x$phi_grid
  y_grid <- x$y_grid

  if (is.null(zlim)) {
    zmax <- max(abs(Z), na.rm = TRUE)
    zlim <- c(-zmax, zmax)
  }

  pal <- diagnostic_palette(201, center = center)

  map_col <- function(z) {
    z <- pmax(zlim[1], pmin(zlim[2], z))
    idx <- 1 + floor((length(pal) - 1) * (z - zlim[1]) / diff(zlim))
    pal[idx]
  }

  # Expand grid
  grd <- expand.grid(phi = phi_grid, y = y_grid)
  grd$z <- as.vector(Z)

  # Simple polar mapping:
  # x = r cos(2*pi*phi), y = r sin(2*pi*phi)
  # If your phi is already in radians, replace 2*pi*phi by phi.
  theta <- 2 * pi * grd$phi
  r <- grd$y

  xcoord <- r * cos(theta)
  ycoord <- r * sin(theta)

  lim <- max(r, na.rm = TRUE) * 1.1

  plot(NA, NA,
       xlim = c(-lim, lim),
       ylim = c(-lim, lim),
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "",
       main = main %||% paste("Circular diagnostic:", x$mode))

  # faint radial circles
  rr <- pretty(y_grid)
  rr <- rr[rr > 0]
  for (rv in rr) {
    symbols(0, 0, circles = rv, inches = FALSE, add = TRUE, fg = "#00000020", bg = NA)
  }

  points(xcoord, ycoord, pch = 15, cex = point_cex, col = map_col(grd$z))

  invisible(x)
}

# ============================================================
# 45. Compare real vs simulated at one angle
#
# Overlaid weighted histogram + ECDF

plot_angle_comparison <- function(sample_emp,
                                  sample_ref,
                                  phi0,
                                  mode = c("interpolate", "window"),
                                  bandwidth = NULL,
                                  bins = 20,
                                  main_prefix = "Angle comparison") {
  mode <- match.arg(mode)

  d_emp <- extract_angle_distribution(
    sample = sample_emp,
    phi0 = phi0,
    bandwidth = bandwidth,
    mode = mode
  )

  d_ref <- extract_angle_distribution(
    sample = sample_ref,
    phi0 = phi0,
    bandwidth = bandwidth,
    mode = mode
  )

  if (nrow(d_emp) == 0L || nrow(d_ref) == 0L) {
    stop("No data available for comparison at requested angle.")
  }

  y1 <- d_emp$excess
  w1 <- d_emp$weight / sum(d_emp$weight)

  y2 <- d_ref$excess
  w2 <- d_ref$weight / sum(d_ref$weight)

  op <- par(mfrow = c(1, 2))
  on.exit(par(op), add = TRUE)

  # --- histogram
  rng <- range(c(y1, y2), na.rm = TRUE)
  breaks <- pretty(rng, n = bins)
  if (length(unique(breaks)) < 2L) {
    breaks <- seq(rng[1], rng[2] + 1e-6, length.out = 10)
  }

  mids <- head(breaks, -1) + diff(breaks) / 2

  id1 <- cut(y1, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  id2 <- cut(y2, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  h1 <- tapply(w1, id1, sum)
  h2 <- tapply(w2, id2, sum)

  nb <- length(breaks) - 1L
  hh1 <- rep(0, nb); hh2 <- rep(0, nb)
  if (!is.null(h1)) hh1[as.integer(names(h1))] <- as.numeric(h1)
  if (!is.null(h2)) hh2[as.integer(names(h2))] <- as.numeric(h2)

  ymax <- max(c(hh1, hh2), na.rm = TRUE)

  plot(mids, hh1, type = "h", lwd = 8,
       ylim = c(0, ymax),
       xlab = "Excess", ylab = "Weighted bin mass",
       main = paste0(main_prefix, "\nphi = ", signif(phi0, 4)))
  lines(mids, hh2, type = "h", lwd = 4, lty = 2)
  legend("topright",
         legend = c("Empirical", "Reference simulation"),
         lty = c(1, 2), lwd = c(8, 4), bty = "n")

  # --- ECDF
  F1 <- weighted_ecdf_fun(y1, w1)
  F2 <- weighted_ecdf_fun(y2, w2)
  xg <- seq(rng[1], rng[2], length.out = 300)

  plot(xg, F1(xg), type = "s", lwd = 2,
       xlab = "Excess", ylab = "F(x)",
       main = "Weighted ECDF")
  lines(xg, F2(xg), type = "s", lwd = 2, lty = 2)
  legend("bottomright",
         legend = c("Empirical", "Reference simulation"),
         lty = c(1, 2), lwd = 2, bty = "n")

  invisible(list(empirical = d_emp, reference = d_ref))
}

# ============================================================
# 46. Build a holistic comparison against simulated reference
#
# This uses the same diagnostic builder, but the reference is now
# empirical-from-simulation rather than analytic.
#
# Useful if you want a like-for-like comparison between:
#   - actual extracted event-level sample
#   - synthetic sample generated from the model

build_empirical_vs_simulated_map <- function(sample_emp,
                                             sample_ref,
                                             phi_grid,
                                             y_grid,
                                             extraction_mode = c("interpolate", "window"),
                                             extraction_bandwidth = NULL,
                                             eps = 1e-8) {
  extraction_mode <- match.arg(extraction_mode)

  Z <- matrix(NA_real_, nrow = length(phi_grid), ncol = length(y_grid))

  for (j in seq_along(phi_grid)) {
    phi0 <- phi_grid[j]

    S_emp <- empirical_survival_at_angle(
      sample = sample_emp,
      phi0 = phi0,
      y_grid = y_grid,
      mode = extraction_mode,
      bandwidth = extraction_bandwidth
    )

    S_ref <- empirical_survival_at_angle(
      sample = sample_ref,
      phi0 = phi0,
      y_grid = y_grid,
      mode = extraction_mode,
      bandwidth = extraction_bandwidth
    )

    Z[j,] <- log((S_emp + eps) / (S_ref + eps))
  }

  new_holistic_diagnostic(
    phi_grid = phi_grid,
    y_grid = y_grid,
    value_matrix = Z,
    mode = "log_survival_ratio",
    metadata = list(reference = "simulated_gp")
  )
}

test_holistic_diagnosis <- function(){
  # ============================================================
  # 36. Example holistic diagnostic

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

  phi_grid_diag <- seq(min(sample_obj$data$phi), max(sample_obj$data$phi), length.out = 60)
  y_grid_diag <- seq(0, max(sample_obj$data$excess, na.rm = TRUE), length.out = 80)

  # Example reference: mild angle-varying GP
  ref <- make_gp_reference(
    scale_fun = function(phi) 0.8 + 0.5 * (sin(2 * pi * phi)^2),
    shape_fun = function(phi) 0.1 + 0.05 * cos(2 * pi * phi)
  )

  hol_diag <- build_holistic_diagnostic(
    sample = sample_obj,
    phi_grid = phi_grid_diag,
    y_grid = y_grid_diag,
    mode = "log_survival_ratio",
    extraction_mode = "interpolate",
    ref_survival = ref$survival
  )

  print(hol_diag)
  plot(hol_diag, style = "square", center = "white", draw_contours = TRUE)
  plot(hol_diag, style = "circular", center = "white")


  # sample_obj <- new_envelope_sample(df_all)

  cal <- compute_threshold_calibration(
    sample = sample_obj,
    phi_grid = seq(min(sample_obj$data$phi), max(sample_obj$data$phi), length.out = 40),
    mode = "interpolate",
    target_rate = 0.2
  )

  print(cal)
  plot(cal)

  bias_tab <- summarize_threshold_bias(cal)
  head(bias_tab)
}


test_map_plot <- function() {

  # ============================================================
  # 47. Example workflow
  #
  # Use the same phi grid as your envelope data if possible.
  # Here I assume:
  #   sample_obj = empirical envelope sample from your real/extracted data

  # Example angle-varying GP parameter functions
  scale_fun_ex <- function(phi) 0.8 + 0.5 * (sin(2 * pi * phi)^2)
  shape_fun_ex <- function(phi) 0.10 + 0.05 * cos(2 * pi * phi)

  # Suppose your empirical sample already exists:
  # sample_obj <- new_envelope_sample(df_all)

  # Use the empirical angle grid for the simulation
  # (or choose a regular grid)
  phi_grid_sim <- sort(unique(sample_obj$data$phi))

  # Simulate, say, the same number of clusters as in the empirical sample
  n_clusters_emp <- length(unique(sample_obj$data$cluster_id))

  sample_ref <- simulate_gp_reference_sample_from_funs(
    phi_grid = phi_grid_sim,
    n_clusters = n_clusters_emp,
    scale_fun = scale_fun_ex,
    shape_fun = shape_fun_ex
  )

  print(sample_ref)

  # ------------------------------------------------------------
  # Single-angle comparison

  phi0 <- median(phi_grid_sim)

  plot_angle_comparison(
    sample_emp = sample_obj,
    sample_ref = sample_ref,
    phi0 = phi0,
    mode = "interpolate"
  )

  # ------------------------------------------------------------
  # Holistic square/circular comparison map
  y_grid_cmp <- seq(
    0,
    max(c(sample_obj$data$excess, sample_ref$data$excess), na.rm = TRUE),
    length.out = 80
  )

  hol_cmp <- build_empirical_vs_simulated_map(
    sample_emp = sample_obj,
    sample_ref = sample_ref,
    phi_grid = seq(min(phi_grid_sim), max(phi_grid_sim), length.out = 60),
    y_grid = y_grid_cmp,
    extraction_mode = "interpolate"
  )

  print(hol_cmp)
  plot(hol_cmp, style = "square", center = "white", draw_contours = TRUE)
  plot(hol_cmp, style = "circular", center = "white")
}