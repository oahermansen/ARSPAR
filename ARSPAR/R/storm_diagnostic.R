#' Summarize declustered storms by time interval
#'
#' @param x A `spar_representation` object.
#' @param interval One of `"year"`, `"n_samples"`, `"time_span"`.
#' @param n_samples Bin width in samples for `interval = "n_samples"`.
#' @param span_length Bin width in numeric time units for `interval = "time_span"`.
#' @param clusters Optional `clustered_excursion_spans` object.
#' @param gap_rule Declustering gap used when clusters are not available.
#' @param use_stored Logical; if `TRUE`, reuse stored clusters when available.
#'
#' @return Data frame with storm counts and span summaries per interval.
#' @export
spar_summarize_storm_intervals <- function(
    x,
    interval = c("year", "n_samples", "time_span"),
    n_samples = 1000L,
    span_length = NULL,
    clusters = NULL,
    gap_rule = 0,
    use_stored = TRUE
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  interval <- match.arg(interval)

  if (is.null(clusters) && isTRUE(use_stored)) {
    clusters <- x$excursions$clusters
  }
  if (is.null(clusters) || !inherits(clusters, "clustered_excursion_spans")) {
    dec <- spar_decluster_excursions(x, gap_rule = gap_rule, store = FALSE)
    clusters <- dec$clusters
  }

  spans <- clusters$spans
  if (nrow(spans) == 0L) {
    return(data.frame(
      bin_id = integer(0),
      bin_label = character(0),
      n_storms = integer(0),
      mean_total_span = numeric(0),
      q95_total_span = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  start_idx <- as.integer(spans$First)

  if (interval == "year") {
    if (is.null(x$data$time)) {
      stop("`x$data$time` is required for `interval = 'year'`.", call. = FALSE)
    }

    tt <- x$data$time
    if (inherits(tt, "Date") || inherits(tt, "POSIXt")) {
      tt_posix <- as.POSIXct(tt, tz = "UTC")
    } else {
      tt_posix <- suppressWarnings(as.POSIXct(tt, tz = "UTC"))
    }
    if (anyNA(tt_posix)) {
      stop("Could not parse `x$data$time` to datetime for yearly binning.", call. = FALSE)
    }

    bins <- format(tt_posix[start_idx], "%Y")
    bin_labels <- bins
  } else if (interval == "n_samples") {
    if (!is.numeric(n_samples) || length(n_samples) != 1L || is.na(n_samples) || n_samples < 1) {
      stop("`n_samples` must be a single positive integer.", call. = FALSE)
    }
    n_samples <- as.integer(n_samples)
    bin_id <- ((start_idx - 1L) %/% n_samples) + 1L
    left <- (bin_id - 1L) * n_samples + 1L
    right <- bin_id * n_samples
    bins <- as.character(bin_id)
    bin_labels <- sprintf("[%d,%d]", left, right)
  } else {
    if (is.null(x$data$time)) {
      stop("`x$data$time` is required for `interval = 'time_span'`.", call. = FALSE)
    }
    if (!is.numeric(span_length) || length(span_length) != 1L || is.na(span_length) || span_length <= 0) {
      stop("`span_length` must be a single positive numeric value.", call. = FALSE)
    }

    tt <- x$data$time
    if (inherits(tt, "Date") || inherits(tt, "POSIXt")) {
      tnum <- as.numeric(as.POSIXct(tt, tz = "UTC"))
      t_origin <- as.POSIXct(min(tt, na.rm = TRUE), tz = "UTC")
      t0 <- as.numeric(t_origin)
      bin_id <- floor((tnum[start_idx] - t0) / span_length) + 1L
      left <- t0 + (bin_id - 1L) * span_length
      right <- left + span_length
      bins <- as.character(bin_id)
      bin_labels <- sprintf("[%s,%s)", format(as.POSIXct(left, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d"), format(as.POSIXct(right, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d"))
    } else {
      tnum <- suppressWarnings(as.numeric(tt))
      if (anyNA(tnum)) {
        stop("`x$data$time` must be numeric/datetime-like for `interval = 'time_span'`.", call. = FALSE)
      }
      t0 <- min(tnum, na.rm = TRUE)
      bin_id <- floor((tnum[start_idx] - t0) / span_length) + 1L
      left <- t0 + (bin_id - 1L) * span_length
      right <- left + span_length
      bins <- as.character(bin_id)
      bin_labels <- sprintf("[%g,%g)", left, right)
    }
  }

  d <- data.frame(
    bin = bins,
    bin_label = bin_labels,
    total_span = as.numeric(spans$total_span),
    stringsAsFactors = FALSE
  )

  split_d <- split(d, d$bin)
  out <- lapply(split_d, function(dd) {
    data.frame(
      bin_id = as.integer(dd$bin[1]),
      bin_label = dd$bin_label[1],
      n_storms = nrow(dd),
      mean_total_span = mean(dd$total_span, na.rm = TRUE),
      q95_total_span = as.numeric(stats::quantile(dd$total_span, probs = 0.95, na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out <- out[order(out$bin_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Plot storm counts by interval
#'
#' @param x Optional `spar_representation` object.
#' @param summary_df Optional precomputed interval summary from
#'   `spar_summarize_storm_intervals()`.
#' @param interval Interval mode passed when `summary_df` is not supplied.
#' @param n_samples Bin width in samples for `interval = "n_samples"`.
#' @param span_length Bin width for `interval = "time_span"`.
#' @param gap_rule Declustering gap used when computing summary.
#' @param use_stored Logical; whether to reuse stored clusters.
#' @param main Plot title.
#'
#' @return The summary data frame (invisibly).
#' @export
spar_plot_storms_by_interval <- function(
    x = NULL,
    summary_df = NULL,
    interval = c("year", "n_samples", "time_span"),
    n_samples = 1000L,
    span_length = NULL,
    gap_rule = 0,
    use_stored = TRUE,
    main = NULL
) {
  interval <- match.arg(interval)

  if (is.null(summary_df)) {
    if (is.null(x)) {
      stop("Provide either `x` or `summary_df`.", call. = FALSE)
    }
    summary_df <- spar_summarize_storm_intervals(
      x = x,
      interval = interval,
      n_samples = n_samples,
      span_length = span_length,
      gap_rule = gap_rule,
      use_stored = use_stored
    )
  }

  if (!is.data.frame(summary_df) || !all(c("bin_label", "n_storms") %in% names(summary_df))) {
    stop("`summary_df` must be a summary returned by `spar_summarize_storm_intervals()`.", call. = FALSE)
  }

  if (is.null(main)) {
    main <- sprintf("Declustered storm count by %s", interval)
  }

  graphics::barplot(
    height = summary_df$n_storms,
    names.arg = summary_df$bin_label,
    las = 2,
    col = "#1f77b4",
    border = NA,
    ylab = "Number of storms",
    main = main
  )

  invisible(summary_df)
}

#' Extract support-based storm quantity by angle
#'
#' @param x A `spar_representation` object.
#' @param phi_grid Optional angular grid.
#' @param n_phi Grid size if `phi_grid` is `NULL`.
#' @param gap_rule Declustering gap.
#' @param use_stored Logical; whether to reuse stored clusters.
#' @param eps Positive threshold for support (`excess > eps`).
#' @param time_span Optional denominator for rates.
#' @param yearly_average Logical; if `TRUE`, append yearly-average support
#'   quantity using `x$data$time`.
#'
#' @return Data frame with support count and rates by angle.
#' @export
spar_storm_quantity_by_angle <- function(
    x,
    phi_grid = NULL,
    n_phi = 180L,
    gap_rule = 0,
    use_stored = TRUE,
    eps = 1e-12,
    time_span = NULL,
    yearly_average = FALSE
) {
  ad <- spar_compute_angular_density(
    x,
    time_span = time_span,
    phi_grid = phi_grid,
    n_phi = n_phi,
    eps = eps,
    source = "cluster_support",
    gap_rule = gap_rule,
    use_stored = use_stored,
    store = FALSE
  )

  out <- data.frame(
    phi = ad$phi_grid,
    support_count = ad$support_count,
    conditional_support = ad$conditional_support,
    declustered_rate = ad$declustered_rate,
    density = ad$density,
    stringsAsFactors = FALSE
  )

  if (isTRUE(yearly_average)) {
    if (is.null(x$data$time)) {
      stop("`x$data$time` is required when `yearly_average = TRUE`.", call. = FALSE)
    }

    tt <- x$data$time
    if (inherits(tt, "Date") || inherits(tt, "POSIXt")) {
      tt_posix <- as.POSIXct(tt, tz = "UTC")
    } else {
      tt_posix <- suppressWarnings(as.POSIXct(tt, tz = "UTC"))
    }
    if (anyNA(tt_posix)) {
      stop("Could not parse `x$data$time` to datetime for yearly averaging.", call. = FALSE)
    }

    years <- unique(format(tt_posix, "%Y"))
    n_years <- length(years)
    if (n_years < 1L) {
      stop("No valid year labels found for yearly averaging.", call. = FALSE)
    }

    out$n_years <- n_years
    out$support_count_yearly_avg <- out$support_count / n_years
  }

  out
}

#' Plot support-based storm quantity by angle
#'
#' @param x Optional `spar_representation` object.
#' @param quantity_df Optional precomputed output from
#'   `spar_storm_quantity_by_angle()`.
#' @param value One of `"support_count"`, `"declustered_rate"`,
#'   `"conditional_support"`, `"density"`, `"support_count_yearly_avg"`.
#' @param n_phi Grid size if computing internally.
#' @param gap_rule Declustering gap.
#' @param use_stored Logical; whether to reuse stored clusters.
#' @param eps Positive threshold for support.
#' @param main Plot title.
#'
#' @return The quantity data frame (invisibly).
#' @export
spar_plot_storm_quantity_by_angle <- function(
    x = NULL,
    quantity_df = NULL,
    value = c("support_count", "declustered_rate", "conditional_support", "density", "support_count_yearly_avg"),
    n_phi = 180L,
    gap_rule = 0,
    use_stored = TRUE,
    eps = 1e-12,
    yearly_average = FALSE,
    main = NULL
) {
  value <- match.arg(value)

  if (is.null(quantity_df)) {
    if (is.null(x)) {
      stop("Provide either `x` or `quantity_df`.", call. = FALSE)
    }
    quantity_df <- spar_storm_quantity_by_angle(
      x = x,
      n_phi = n_phi,
      gap_rule = gap_rule,
      use_stored = use_stored,
      eps = eps,
      yearly_average = yearly_average
    )
  }

  if (!is.data.frame(quantity_df) || !all(c("phi", value) %in% names(quantity_df))) {
    stop("`quantity_df` must contain `phi` and the selected `value` column.", call. = FALSE)
  }

  ylab <- switch(
    value,
    support_count = "Storm support count",
    declustered_rate = "Declustered storm rate",
    conditional_support = "Conditional support",
    density = "Angular density",
    support_count_yearly_avg = "Average storm support count per year"
  )

  if (is.null(main)) {
    main <- sprintf("Storm quantity by angle (%s)", value)
  }

  yvals <- quantity_df[[value]]
  y_rng <- range(yvals, na.rm = TRUE)
  if (!all(is.finite(y_rng))) {
    y_rng <- c(0, 1)
  }
  if (identical(value, "support_count_yearly_avg")) {
    y_rng[1] <- 0
  }

  graphics::plot(
    quantity_df$phi,
    quantity_df[[value]],
    type = "l",
    lwd = 2,
    col = "#d62728",
    xlab = "Angle",
    ylab = ylab,
    ylim = y_rng,
    main = main
  )

  invisible(quantity_df)
}

#' Compute upper-path excess quantiles by angle
#'
#' Uses all upper excursion path observations within an angular window around
#' each grid angle.
#'
#' @param x A `spar_representation` object.
#' @param phi_grid Optional angular grid.
#' @param n_phi Grid size if `phi_grid` is `NULL`.
#' @param bandwidth Angular half-window. If `NULL`, uses
#'   `bw_mult * median(diff(phi_grid))`.
#' @param bw_mult Multiplier for default bandwidth.
#' @param probs Quantile probabilities.
#' @param upper_paths Optional `excursion_upper_path_group`.
#' @param gap_rule Declustering gap used if upper paths must be built.
#' @param use_stored Logical; reuse stored upper paths when available.
#' @param eps Positive threshold (`excess > eps`) for usable values.
#'
#' @return Data frame with `phi`, `n_local`, quantile columns, and an attribute
#'   `insufficient_report` containing angle/quantile combinations with
#'   insufficient local sample size.
#' @export
spar_upper_path_excess_quantiles_by_angle <- function(
    x,
    phi_grid = NULL,
    n_phi = 180L,
    bandwidth = NULL,
    bw_mult = 2,
    probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
    upper_paths = NULL,
    gap_rule = 0,
    use_stored = TRUE,
    eps = 1e-12
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  probs <- as.numeric(probs)
  if (length(probs) == 0L || any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("`probs` must contain finite values strictly between 0 and 1.", call. = FALSE)
  }

  dom <- x$angular$domain
  if (!inherits(dom, "spar_angle_domain")) {
    dom <- NULL
  }

  if (is.null(phi_grid)) {
    n_phi <- as.integer(n_phi)
    if (!is.finite(n_phi) || length(n_phi) != 1L || n_phi < 2L) {
      stop("`n_phi` must be an integer >= 2.", call. = FALSE)
    }

    if (!is.null(dom)) {
      phi_grid <- spar_angle_grid(n = n_phi, domain = dom)
    } else if (!is.null(x$angular$phi)) {
      p <- as.numeric(x$angular$phi)
      p <- p[is.finite(p)]
      if (length(p) == 0L) {
        stop("Cannot infer `phi_grid`; no finite angles available.", call. = FALSE)
      }
      phi_grid <- seq(min(p), max(p), length.out = n_phi)
    } else {
      stop("Cannot infer `phi_grid`; provide `phi_grid`.", call. = FALSE)
    }
  } else {
    phi_grid <- as.numeric(phi_grid)
    if (length(phi_grid) < 2L || any(!is.finite(phi_grid))) {
      stop("`phi_grid` must be a numeric vector with at least two finite values.", call. = FALSE)
    }
  }

  if (is.null(dom)) {
    dom <- new_spar_angle_domain(type = "interval", lower = min(phi_grid), upper = max(phi_grid))
  }
  dom <- as_spar_angle_domain(dom)

  if (is.null(bandwidth)) {
    d <- diff(sort(unique(as.numeric(phi_grid))))
    d <- d[is.finite(d) & d > 0]
    step <- if (length(d) > 0L) median(d) else 1e-6
    bandwidth <- bw_mult * step
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || is.na(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a single positive numeric value.", call. = FALSE)
  }

  if (is.null(upper_paths) && isTRUE(use_stored)) {
    upper_paths <- x$excursions$upper_paths
  }
  if (is.null(upper_paths) || !inherits(upper_paths, "excursion_upper_path_group")) {
    upper_paths <- spar_build_upper_excursion_paths(
      x,
      path_group = NULL,
      gap_rule = gap_rule,
      store = FALSE
    )
  }

  if (length(upper_paths$paths) == 0L) {
    out <- data.frame(phi = phi_grid, n_local = 0L, stringsAsFactors = FALSE)
    for (p in probs) {
      nm <- sprintf("q%02d", as.integer(round(100 * p)))
      out[[nm]] <- NA_real_
    }
    attr(out, "insufficient_report") <- data.frame(
      phi = numeric(0),
      quantile = character(0),
      n_local = integer(0),
      required_n = integer(0),
      stringsAsFactors = FALSE
    )
    return(out)
  }

  phi_all <- unlist(lapply(upper_paths$paths, function(p) as.numeric(p$phi)), use.names = FALSE)
  y_all <- unlist(lapply(upper_paths$paths, function(p) as.numeric(p$excess)), use.names = FALSE)

  ok <- is.finite(phi_all) & is.finite(y_all) & (y_all > eps)
  phi_all <- phi_all[ok]
  y_all <- y_all[ok]

  q_names <- sprintf("q%02d", as.integer(round(100 * probs)))
  if (any(duplicated(q_names))) {
    stop("`probs` generates duplicate quantile names; use distinct probabilities.", call. = FALSE)
  }

  q_mat <- matrix(NA_real_, nrow = length(phi_grid), ncol = length(probs))
  colnames(q_mat) <- q_names
  n_local <- integer(length(phi_grid))
  report_rows <- vector("list", length(phi_grid) * length(probs))
  rr <- 0L

  for (j in seq_along(phi_grid)) {
    phi0 <- phi_grid[j]
    dphi <- spar_angle_distance(phi_all, phi0, dom)
    keep <- is.finite(dphi) & (dphi <= bandwidth)
    yj <- y_all[keep]
    n_local[j] <- length(yj)

    for (k in seq_along(probs)) {
      p <- probs[k]
      req_n <- max(1L, as.integer(ceiling(1 / min(p, 1 - p))))
      if (n_local[j] < req_n) {
        rr <- rr + 1L
        report_rows[[rr]] <- data.frame(
          phi = phi0,
          quantile = sprintf("q%02d", as.integer(round(100 * p))),
          n_local = n_local[j],
          required_n = req_n,
          stringsAsFactors = FALSE
        )
        q_mat[j, k] <- NA_real_
      } else {
        q_mat[j, k] <- as.numeric(stats::quantile(yj, probs = p, names = FALSE, na.rm = TRUE, type = 8))
      }
    }
  }

  insufficient_report <- if (rr > 0L) do.call(rbind, report_rows[seq_len(rr)]) else data.frame(
    phi = numeric(0),
    quantile = character(0),
    n_local = integer(0),
    required_n = integer(0),
    stringsAsFactors = FALSE
  )

  out <- data.frame(phi = phi_grid, n_local = n_local, q_mat, stringsAsFactors = FALSE)

  if (dom$type == "cyclical" && length(phi_grid) >= 2L) {
    is_wrapped_endpoint <- is.finite(dom$width) &&
      isTRUE(all.equal(phi_grid[length(phi_grid)] - phi_grid[1], dom$width, tolerance = 1e-8))
    if (is_wrapped_endpoint) {
      out[nrow(out), setdiff(names(out), "phi")] <- out[1, setdiff(names(out), "phi")]
    }
  }

  attr(out, "insufficient_report") <- insufficient_report
  attr(out, "bandwidth") <- bandwidth
  attr(out, "bw_mult") <- bw_mult
  if (nrow(insufficient_report) > 0L) {
    warning(sprintf(
      "Insufficient local data for %d angle-quantile combinations. See attr(result, 'insufficient_report').",
      nrow(insufficient_report)
    ), call. = FALSE)
  }

  out
}

#' Build support and upper-path quantile diagnostic table by angle
#'
#' @param x A `spar_representation` object.
#' @param phi_grid Optional angular grid.
#' @param n_phi Grid size if `phi_grid` is `NULL`.
#' @param gap_rule Declustering gap.
#' @param use_stored Logical; whether to reuse stored clusters/upper paths.
#' @param eps Positive threshold (`excess > eps`).
#' @param probs Quantile probabilities.
#' @param bandwidth Angular half-window for local quantiles.
#' @param bw_mult Bandwidth multiplier when `bandwidth = NULL`.
#'
#' @return Data frame with yearly average storm support and upper-path excess
#'   quantiles by angle. Attribute `insufficient_report` contains local-sample
#'   insufficiency details.
#' @export
spar_storm_excess_dual_axis_by_angle <- function(
    x,
    phi_grid = NULL,
    n_phi = 180L,
    gap_rule = 0,
    use_stored = TRUE,
    eps = 1e-12,
    probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
    bandwidth = NULL,
    bw_mult = 2
) {
  support_df <- spar_storm_quantity_by_angle(
    x = x,
    phi_grid = phi_grid,
    n_phi = n_phi,
    gap_rule = gap_rule,
    use_stored = use_stored,
    eps = eps,
    yearly_average = TRUE
  )

  q_df <- spar_upper_path_excess_quantiles_by_angle(
    x = x,
    phi_grid = support_df$phi,
    n_phi = n_phi,
    bandwidth = bandwidth,
    bw_mult = bw_mult,
    probs = probs,
    gap_rule = gap_rule,
    use_stored = use_stored,
    eps = eps
  )

  out <- merge(support_df, q_df, by = "phi", all = TRUE, sort = TRUE)
  attr(out, "insufficient_report") <- attr(q_df, "insufficient_report")
  out$bandwidth <- as.numeric(attr(q_df, "bandwidth"))
  out$bw_mult <- as.numeric(attr(q_df, "bw_mult"))
  out
}

#' Plot yearly storm support and upper-path excess quantiles by angle
#'
#' @param x Optional `spar_representation` object.
#' @param diag_df Optional data frame from
#'   `spar_storm_excess_dual_axis_by_angle()`.
#' @param n_phi Grid size if computing internally.
#' @param gap_rule Declustering gap.
#' @param probs Quantile probabilities when computing internally.
#' @param bandwidth Angular half-window for local quantiles.
#' @param bw_mult Bandwidth multiplier when `bandwidth = NULL`.
#' @param support_col Column name for left-axis storm quantity.
#' @param main Plot title.
#'
#' @return The diagnostic table (invisibly).
#' @export
spar_plot_storm_excess_dual_axis <- function(
    x = NULL,
    diag_df = NULL,
    n_phi = 180L,
    gap_rule = 0,
    probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
    bandwidth = NULL,
    bw_mult = 2,
    support_col = "support_count_yearly_avg",
    main = "Yearly storm support and excess quantiles by angle"
) {
  if (is.null(diag_df)) {
    if (is.null(x)) {
      stop("Provide either `x` or `diag_df`.", call. = FALSE)
    }
    diag_df <- spar_storm_excess_dual_axis_by_angle(
      x = x,
      n_phi = n_phi,
      gap_rule = gap_rule,
      probs = probs,
      bandwidth = bandwidth,
      bw_mult = bw_mult
    )
  }

  if (!is.data.frame(diag_df) || !all(c("phi", support_col) %in% names(diag_df))) {
    stop("`diag_df` must contain `phi` and selected `support_col`.", call. = FALSE)
  }

  q_cols <- grep("^q[0-9]{2}$", names(diag_df), value = TRUE)
  if (length(q_cols) == 0L) {
    stop("`diag_df` does not contain quantile columns (`q..`).", call. = FALSE)
  }

  y_left <- diag_df[[support_col]]
  y_left_rng <- range(y_left, na.rm = TRUE)
  if (!all(is.finite(y_left_rng))) {
    y_left_rng <- c(0, 1)
  }
  y_left_rng[1] <- 0

  xvals <- as.numeric(diag_df$phi)
  x_rng <- range(xvals, na.rm = TRUE)
  x_minor <- seq(from = floor(x_rng[1] / 0.1) * 0.1, to = ceiling(x_rng[2] / 0.1) * 0.1, by = 0.1)
  x_major <- seq(from = floor(x_rng[1] / 0.2) * 0.2, to = ceiling(x_rng[2] / 0.2) * 0.2, by = 0.2)

  graphics::plot(
    diag_df$phi,
    y_left,
    type = "l",
    lwd = 2.5,
    col = "#1f77b4",
    xlab = "Angle",
    ylab = "Average yearly storm support count",
    xlim = x_rng,
    ylim = y_left_rng,
    main = main,
    axes = FALSE
  )

  y_left_grid <- pretty(y_left_rng, n = 8)
  y_left_grid <- y_left_grid[y_left_grid >= y_left_rng[1] & y_left_grid <= y_left_rng[2]]
  graphics::abline(v = x_minor, col = "#f0f0f0", lwd = 1)
  graphics::abline(v = x_major, col = "#d9d9d9", lwd = 1.2)
  graphics::abline(h = y_left_grid, col = "#ececec", lwd = 1)
  graphics::axis(1, at = x_major)
  graphics::axis(2)
  graphics::box()

  q_vals <- unlist(diag_df[q_cols], use.names = FALSE)
  y_right_rng <- range(q_vals, na.rm = TRUE)
  if (!all(is.finite(y_right_rng))) {
    y_right_rng <- c(0, 1)
  }

  graphics::par(new = TRUE)
  graphics::plot(
    diag_df$phi,
    diag_df[[q_cols[1]]],
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = "",
    xlim = x_rng,
    ylim = y_right_rng
  )

  q_cols_plot <- q_cols[vapply(q_cols, function(nm) any(is.finite(diag_df[[nm]])), logical(1))]
  q_cols_skip <- setdiff(q_cols, q_cols_plot)
  q_colors <- grDevices::hcl.colors(max(length(q_cols_plot), 1L), palette = "Dark 3")

  if (length(q_cols_plot) > 0L) {
    for (i in seq_along(q_cols_plot)) {
      graphics::lines(diag_df$phi, diag_df[[q_cols_plot[i]]], lwd = 1.8, col = q_colors[i])
    }
  }

  graphics::axis(4)
  graphics::mtext("Upper-path excess quantiles", side = 4, line = 3)

  q_labels <- if (length(q_cols_plot) > 0L) {
    paste0("Excess ", sub("^q", "Q", q_cols_plot), "%")
  } else {
    character(0)
  }

  legend_labels <- c("Yearly storm support", q_labels)
  legend_cols <- c("#1f77b4", if (length(q_cols_plot) > 0L) q_colors else character(0))
  legend_lwd <- c(2.5, if (length(q_cols_plot) > 0L) rep(1.8, length(q_cols_plot)) else numeric(0))
  graphics::legend("topright", legend = legend_labels, col = legend_cols, lty = 1, lwd = legend_lwd, bty = "n")

  if (length(q_cols_skip) > 0L) {
    warning(sprintf(
      "Quantile lines with all-NA values skipped in plot: %s",
      paste(q_cols_skip, collapse = ", ")
    ), call. = FALSE)
  }

  invisible(diag_df)
}
