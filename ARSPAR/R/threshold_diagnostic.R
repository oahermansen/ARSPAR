# ============================================================
# 38. Event-level threshold calibration diagnostic
#
# Assumes the sample data contain one row per cluster-angle grid
# point, with:
#   phi      = angle
#   excess   = value relative to fitted threshold
#              (positive means exceedance)
#   weight   = cluster/angular quadrature weight

new_threshold_calibration <- function(phi_grid,
                                      rate,
                                      target_rate,
                                      n_eff = NULL,
                                      metadata = list()) {
  stopifnot(length(phi_grid) == length(rate))

  structure(
    list(
      phi_grid = phi_grid,
      rate = rate,
      target_rate = target_rate,
      n_eff = n_eff,
      metadata = metadata
    ),
    class = "threshold_calibration"
  )
}

print.threshold_calibration <- function(x, ...) {
  cat("threshold_calibration\n")
  cat("  target exceedance rate :", signif(x$target_rate, 6), "\n")
  cat("  angle grid size        :", length(x$phi_grid), "\n")
  cat("  empirical rate range   :",
      signif(min(x$rate, na.rm = TRUE), 6), "to",
      signif(max(x$rate, na.rm = TRUE), 6), "\n")
  cat("  mean abs deviation     :",
      signif(mean(abs(x$rate - x$target_rate), na.rm = TRUE), 6), "\n")

  worst <- order(abs(x$rate - x$target_rate), decreasing = TRUE)
  worst <- worst[seq_len(min(5, length(worst)))]

  cat("  worst angles:\n")
  for (i in worst) {
    cat("    phi = ", signif(x$phi_grid[i], 6),
        "  rate = ", signif(x$rate[i], 6),
        "  dev = ", signif(x$rate[i] - x$target_rate, 6), "\n", sep = "")
  }

  invisible(x)
}

# ============================================================
# 39. Weighted exceedance-rate-by-angle
# Uses the envelope sample directly.
# excess > 0 means event-level exceedance relative to threshold.

compute_threshold_calibration <- function(sample,
                                          phi_grid = NULL,
                                          mode = c("interpolate", "window"),
                                          bandwidth = NULL,
                                          target_rate = 0.2,
                                          progress = FALSE) {
  stopifnot(inherits(sample, "envelope_sample"))
  mode <- match.arg(mode)

  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("`progress` must be a single non-missing logical value.", call. = FALSE)
  }

  df <- sample$data
  angle_name <- sample$angle_name

  if (is.null(phi_grid)) {
    phi_grid <- sort(unique(df[[angle_name]]))
  }

  rate <- rep(NA_real_, length(phi_grid))
  n_eff <- rep(NA_real_, length(phi_grid))

  pb <- NULL
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = length(phi_grid), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (j in seq_along(phi_grid)) {
    phi0 <- phi_grid[j]

    d0 <- extract_angle_distribution(
      sample = sample,
      phi0 = phi0,
      bandwidth = bandwidth,
      mode = mode
    )

    if (nrow(d0) == 0L) next

    y <- d0$excess
    w <- d0$weight

    ok <- is.finite(y) & is.finite(w) & w > 0
    y <- y[ok]
    w <- w[ok]

    if (length(y) == 0L) next

    rate[j] <- sum(w[y > 0]) / sum(w)

    # Kish effective sample size
    n_eff[j] <- (sum(w)^2) / sum(w^2)

    if (isTRUE(progress)) {
      utils::setTxtProgressBar(pb, j)
    }
  }

  new_threshold_calibration(
    phi_grid = phi_grid,
    rate = rate,
    target_rate = target_rate,
    n_eff = n_eff,
    metadata = list(mode = mode, bandwidth = bandwidth)
  )
}

# ============================================================
# 40. Plot threshold calibration

plot.threshold_calibration <- function(x,
                                       show_effective_n = TRUE,
                                       main = "Event-level threshold calibration",
                                       xlab = "Angle",
                                       ylab = "Empirical exceedance rate",
                                       ...) {
  op <- par(mfrow = if (show_effective_n) c(1, 2) else c(1, 1))
  on.exit(par(op), add = TRUE)

  plot(
    x$phi_grid, x$rate,
    type = "b", pch = 19,
    ylim = range(c(x$rate, x$target_rate), na.rm = TRUE),
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
  abline(h = x$target_rate, lty = 2, lwd = 2)
  abline(h = mean(x$rate, na.rm = TRUE), lty = 3)

  if (show_effective_n) {
    plot(
      x$phi_grid, x$n_eff,
      type = "b", pch = 19,
      xlab = "Angle",
      ylab = "Effective sample size",
      main = "Support by angle"
    )
  }

  invisible(x)
}

# ============================================================
# 41. Threshold bias summary
#
# If empirical exceedance rate > target, threshold is effectively
# too low at that angle. If lower than target, threshold is too high.

summarize_threshold_bias <- function(obj) {
  stopifnot(inherits(obj, "threshold_calibration"))

  data.frame(
    phi = obj$phi_grid,
    empirical_rate = obj$rate,
    target_rate = obj$target_rate,
    discrepancy = obj$rate - obj$target_rate,
    interpretation = ifelse(
      obj$rate > obj$target_rate,
      "threshold effectively too low",
      ifelse(obj$rate < obj$target_rate,
             "threshold effectively too high",
             "aligned")
    ),
    n_eff = obj$n_eff
  )
}
