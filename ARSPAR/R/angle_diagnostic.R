
# ============================================================
# 23. Diagnostic object
new_angle_domain <- function(type = c("interval", "circular"),
                             min = NULL,
                             max = NULL,
                             period = NULL,
                             label = NULL) {
  type <- match.arg(type)

  if (type == "interval") {
    dom <- new_spar_angle_domain(
      type = "interval",
      lower = min,
      upper = max,
      label = label %||% "angle"
    )
  } else {
    if ((is.null(min) || is.null(max)) && is.finite(period)) {
      min <- 0
      max <- period
    }

    dom <- new_spar_angle_domain(
      type = "cyclical",
      lower = min,
      upper = max,
      label = label %||% "angle"
    )
  }

  structure(
    list(
      type = if (dom$type == "cyclical") "circular" else "interval",
      min = dom$lower,
      max = dom$upper,
      period = if (dom$type == "cyclical") dom$width else NA_real_,
      label = dom$label
    ),
    class = "angle_domain"
  )
}

angle_range <- function(phi, domain) {
  spar_angle_range(phi, as_spar_angle_domain(domain))
}

angle_in_range <- function(phi, range, domain, closed = TRUE) {
  spar_angle_in_range(phi, range, as_spar_angle_domain(domain), closed = closed)
}

angle_span <- function(phi, domain) {
  spar_angle_span(phi, as_spar_angle_domain(domain))
}

angle_grid <- function(n = 100, domain) {
  spar_angle_grid(n = n, domain = as_spar_angle_domain(domain))
}

angle_delta <- function(phi, phi0, domain) {
  spar_angle_delta(phi, phi0, as_spar_angle_domain(domain))
}

normalize_angle <- function(phi, domain) {
  spar_normalize_angle(phi, as_spar_angle_domain(domain))
}
angle_distance <- function(phi, phi0, domain) {
  spar_angle_distance(phi, phi0, as_spar_angle_domain(domain))
}
angle_gt <- function(phi, phi0, domain) {
  spar_angle_gt(phi, phi0, as_spar_angle_domain(domain))
}
angle_lt <- function(phi, phi0, domain) {
  spar_angle_lt(phi, phi0, as_spar_angle_domain(domain))
}

new_angle_diagnostic <- function(df, phi0, mode, bandwidth = NA_real_, metadata = list()) {
  structure(
    list(
      data = df,
      phi0 = phi0,
      mode = mode,
      bandwidth = bandwidth,
      metadata = metadata
    ),
    class = "angle_diagnostic"
  )
}

print.angle_diagnostic <- function(x, ...) {
  df <- x$data

  cat("angle_diagnostic\n")
  cat("  target angle        :", signif(x$phi0, 6), "\n")
  cat("  extraction mode     :", x$mode, "\n")
  if (is.finite(x$bandwidth)) {
    cat("  bandwidth           :", signif(x$bandwidth, 6), "\n")
  }

  if (nrow(df) == 0L) {
    cat("  observations        : 0\n")
    return(invisible(x))
  }

  y <- df$excess
  w <- df$weight

  qs <- weighted_quantile_safe(y, w, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))

  cat("  observations        :", nrow(df), "\n")
  cat("  positive mass count :", sum(y > 0, na.rm = TRUE), "\n")
  cat("  total weight        :", signif(sum(w, na.rm = TRUE), 6), "\n")
  cat("  weighted mean       :", signif(weighted_mean_safe(y, w), 6), "\n")
  cat("  weighted sd         :", signif(sqrt(weighted_var_safe(y, w)), 6), "\n")
  cat("  weighted quantiles  :\n")
  for (nm in names(qs)) {
    cat("    ", format(nm, width = 4), ": ", signif(qs[[nm]], 6), "\n", sep = "")
  }

  invisible(x)
}

# ============================================================
# 24. Main diagnostic constructor

diagnose_angle_distribution <- function(sample,
                                        phi0,
                                        bandwidth = NULL,
                                        mode = c("interpolate", "window")) {
  mode <- match.arg(mode)

  df <- extract_angle_distribution(
    sample = sample,
    phi0 = phi0,
    bandwidth = bandwidth,
    mode = mode
  )

  new_angle_diagnostic(
    df = df,
    phi0 = attr(df, "phi0"),
    mode = attr(df, "mode"),
    bandwidth = attr(df, "bandwidth")
  )
}

# ============================================================
# 25. Plot diagnostics
#
# Produces:
#   1. weighted histogram
#   2. weighted ECDF
#   3. mean excess over sub-thresholds of the extracted sample

plot.angle_diagnostic <- function(x, bins = 20, ...) {
  df <- x$data

  if (nrow(df) == 0L) {
    plot.new()
    title(main = "No data for requested angle")
    return(invisible(x))
  }

  y <- df$excess
  w <- df$weight

  op <- par(mfrow = c(1, 3))
  on.exit(par(op), add = TRUE)

  # 1. Weighted histogram
  breaks <- pretty(range(y, na.rm = TRUE), n = bins)
  if (length(unique(breaks)) < 2L) {
    breaks <- seq(min(y), max(y) + 1e-6, length.out = 5)
  }

  hist_info <- hist(y, breaks = breaks, plot = FALSE)
  mids <- hist_info$mids
  bin_id <- cut(y, breaks = hist_info$breaks, include.lowest = TRUE, labels = FALSE)
  nb <- length(hist_info$breaks) - 1

  bin_w <- numeric(nb)

  tmp <- tapply(w, bin_id, sum)

  if (!is.null(tmp)) {
    idx <- as.integer(names(tmp))
    bin_w[idx] <- tmp
  }

  barplot(
    height = bin_w,
    names.arg = signif(mids, 3),
    space = 0,
    xlab = "Excess",
    ylab = "Weighted count",
    las = 2,           # vertical labels
    cex.names = 0.7,   # smaller text
    main = paste0("Weighted histogram\nphi = ", signif(x$phi0, 4))
  )

  # 2. Weighted ECDF
  Fhat <- weighted_ecdf_fun(y, w)
  xgrid <- seq(min(y, na.rm = TRUE), max(y, na.rm = TRUE), length.out = 300)
  plot(
    xgrid, Fhat(xgrid), type = "s",
    xlab = "Excess",
    ylab = "F(x)",
    main = "Weighted ECDF"
  )

  # 3. Mean residual life style diagnostic
  u_grid <- unique(quantile(y, probs = seq(0.1, 0.9, by = 0.1), na.rm = TRUE))
  mrl <- vapply(u_grid, function(u0) {
    keep <- y > u0
    if (!any(keep)) return(NA_real_)
    weighted_mean_safe(y[keep] - u0, w[keep])
  }, numeric(1))

  plot(
    u_grid, mrl, type = "b", pch = 19,
    xlab = "Sub-threshold u",
    ylab = "E[Y-u | Y>u]",
    main = "Mean residual life"
  )

  invisible(x)
}

# ============================================================
# 26. Optional: quick comparison across several angles
#
# Returns diagnostics for multiple target angles and prints a
# compact summary table.

diagnose_multiple_angles <- function(sample,
                                     phi_vec,
                                     bandwidth = NULL,
                                     mode = c("interpolate", "window")) {
  mode <- match.arg(mode)

  diags <- lapply(phi_vec, function(phi0) {
    diagnose_angle_distribution(
      sample = sample,
      phi0 = phi0,
      bandwidth = bandwidth,
      mode = mode
    )
  })

  tab <- do.call(rbind, lapply(diags, function(d) {
    df <- d$data
    if (nrow(df) == 0L) {
      return(data.frame(
        phi = d$phi0,
        n = 0,
        mean = NA_real_,
        q50 = NA_real_,
        q90 = NA_real_
      ))
    }

    y <- df$excess
    w <- df$weight
    qs <- weighted_quantile_safe(y, w, probs = c(0.5, 0.9))

    data.frame(
      phi = d$phi0,
      n = nrow(df),
      mean = weighted_mean_safe(y, w),
      q50 = qs[1],
      q90 = qs[2]
    )
  }))

  rownames(tab) <- NULL
  list(diagnostics = diags, summary = tab)
}

# ============================================================
# 21. Weighted summary helpers

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0L) return(NA_real_)
  sum(w * x) / sum(w)
}

weighted_var_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) <= 1L) return(NA_real_)
  mu <- sum(w * x) / sum(w)
  sum(w * (x - mu)^2) / sum(w)
}

weighted_quantile_safe <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0L) {
    out <- rep(NA_real_, length(probs))
    names(out) <- paste0(probs * 100, "%")
    return(out)
  }

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)

  qfun <- function(p) {
    idx <- which(cw >= p)[1]
    x[idx]
  }

  out <- vapply(probs, qfun, numeric(1))
  names(out) <- paste0(probs * 100, "%")
  out
}

weighted_ecdf_fun <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]

  if (length(x) == 0L) {
    return(function(t) rep(NA_real_, length(t)))
  }

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)

  function(t) {
    vapply(t, function(tt) {
      idx <- max(which(x <= tt), 0L)
      if (idx == 0L) 0 else cw[idx]
    }, numeric(1))
  }
}

# ============================================================
# 22. Angle-specific extraction
#
# Two modes:
#   mode = "window"
#       use all rows with |phi - phi0| <= bandwidth
#   mode = "interpolate"
#       for each cluster, linearly interpolate the envelope at phi0
#
# The interpolation mode is often preferable because it gives one
# value per cluster at the target angle.

extract_angle_distribution <- function(sample,
                                       phi0,
                                       bandwidth = NULL,
                                       mode = c("interpolate", "window")) {
  stopifnot(inherits(sample, "envelope_sample"))
  mode <- match.arg(mode)

  df <- sample$data
  angle_name <- sample$angle_name
  response_name <- sample$response_name
  cluster_name <- sample$cluster_name
  weight_name <- sample$weight_name
  angle_domain <- if (is.null(sample$angle_domain)) {
    new_spar_angle_domain()
  } else {
    as_spar_angle_domain(sample$angle_domain)
  }
  phi0 <- as.numeric(phi0)
  if (length(phi0) != 1L || !is.finite(phi0)) {
    stop("`phi0` must be a single finite numeric value.", call. = FALSE)
  }

  if (mode == "window") {
    if (is.null(bandwidth) || bandwidth <= 0) {
      stop("For mode='window', provide a positive bandwidth.")
    }

    keep <- spar_angle_distance(df[[angle_name]], phi0, angle_domain) <= bandwidth
    out <- df[keep, , drop = FALSE]

    attr(out, "phi0") <- phi0
    attr(out, "bandwidth") <- bandwidth
    attr(out, "mode") <- mode
    return(out)
  }

  # mode = "interpolate"
  split_df <- split(df, df[[cluster_name]])
  rows <- vector("list", length(split_df))
  k <- 0L

  for (nm in names(split_df)) {
    dfi <- split_df[[nm]]
    phi <- dfi[[angle_name]]
    y   <- dfi[[response_name]]
    w   <- dfi[[weight_name]]

    ok <- is.finite(phi) & is.finite(y)
    phi <- phi[ok]
    y   <- y[ok]
    w   <- w[ok]

    if (length(phi) < 1L) next

    phi <- spar_normalize_angle(phi, angle_domain)
    phi0_use <- spar_normalize_angle(phi0, angle_domain)

    ord <- order(phi)
    phi <- phi[ord]
    y   <- y[ord]
    w   <- w[ord]

    angdist <- spar_angle_distance(phi, phi0_use, angle_domain)
    if (any(angdist < 1e-12)) {
      idx <- which.min(angdist)
      y0 <- y[idx]
    } else {
      if (length(phi) < 2L) next

      if (angle_domain$type == "interval") {
        if (phi0_use < min(phi) || phi0_use > max(phi)) {
          next
        }

        j <- max(which(phi < phi0_use))
        k2 <- min(which(phi > phi0_use))
        if (!is.finite(j) || !is.finite(k2) || j == k2) next

        t <- (phi0_use - phi[j]) / (phi[k2] - phi[j])
        y0 <- y[j] + t * (y[k2] - y[j])
      } else {
        phi_ext <- c(phi, phi[1] + angle_domain$width)
        y_ext <- c(y, y[1])

        phi0_ext <- phi0_use
        if (phi0_ext < phi[1]) {
          phi0_ext <- phi0_ext + angle_domain$width
        }

        j <- max(which(phi_ext <= phi0_ext))
        if (!is.finite(j) || j >= length(phi_ext)) next
        k2 <- j + 1L

        t <- (phi0_ext - phi_ext[j]) / (phi_ext[k2] - phi_ext[j])
        y0 <- y_ext[j] + t * (y_ext[k2] - y_ext[j])
      }
    }

    # one observation per cluster, equal cluster weight by default
    k <- k + 1L
    rows[[k]] <- data.frame(
      cluster_id = dfi[[cluster_name]][1],
      phi = phi0,
      excess = y0,
      weight = 1
    )
  }

  out <- if (k == 0L) {
    data.frame(cluster_id = integer(0), phi = numeric(0), excess = numeric(0), weight = numeric(0))
  } else {
    do.call(rbind, rows[seq_len(k)])
  }

  attr(out, "phi0") <- phi0
  attr(out, "bandwidth") <- NA_real_
  attr(out, "mode") <- mode
  out
}
