#' Construct angular density diagnostic object
#'
#' @param phi_grid Numeric angle grid.
#' @param support_count Number of supporting clusters per angle.
#' @param conditional_support Conditional support proportion per angle.
#' @param declustered_rate Time-normalized declustered support rate per angle.
#' @param density Normalized angular density integrating to one.
#' @param global_rate Global declustered event rate.
#' @param time_span Time-span denominator used for rates.
#' @param metadata Optional metadata.
#'
#' @return An object of class `"angular_density_diagnostic"`.
#'
#' @export
new_angular_density_diagnostic <- function(
    phi_grid,
    support_count,
    conditional_support,
    declustered_rate,
    density,
    global_rate,
    time_span,
    metadata = list()
) {
  stopifnot(
    is.numeric(phi_grid),
    is.numeric(support_count),
    is.numeric(conditional_support),
    is.numeric(declustered_rate),
    is.numeric(density),
    length(phi_grid) == length(support_count),
    length(phi_grid) == length(conditional_support),
    length(phi_grid) == length(declustered_rate),
    length(phi_grid) == length(density),
    is.numeric(global_rate), length(global_rate) == 1L,
    is.numeric(time_span), length(time_span) == 1L
  )

  structure(
    list(
      phi_grid = phi_grid,
      support_count = support_count,
      conditional_support = conditional_support,
      declustered_rate = declustered_rate,
      density = density,
      global_rate = global_rate,
      time_span = time_span,
      metadata = metadata
    ),
    class = "angular_density_diagnostic"
  )
}

#' Print angular density diagnostic
#'
#' @param x An `angular_density_diagnostic` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.angular_density_diagnostic <- function(x, ...) {
  cat("angular_density_diagnostic\n")
  cat("  angle grid size           :", length(x$phi_grid), "\n")
  cat("  time span                 :", signif(x$time_span, 6), "\n")
  cat("  global declustered rate   :", signif(x$global_rate, 6), "\n")
  cat("  support count range       :", min(x$support_count, na.rm = TRUE), "to", max(x$support_count, na.rm = TRUE), "\n")
  cat("  conditional support range :",
      signif(min(x$conditional_support, na.rm = TRUE), 6), "to",
      signif(max(x$conditional_support, na.rm = TRUE), 6), "\n")
  cat("  declustered rate range    :",
      signif(min(x$declustered_rate, na.rm = TRUE), 6), "to",
      signif(max(x$declustered_rate, na.rm = TRUE), 6), "\n")
  invisible(x)
}

#' Plot angular density diagnostic
#'
#' @param x An `angular_density_diagnostic` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
plot.angular_density_diagnostic <- function(x, ...) {
  op <- par(mfrow = c(1, 2))
  on.exit(par(op), add = TRUE)

  plot(
    x$phi_grid,
    x$declustered_rate,
    type = "l",
    lwd = 2,
    xlab = "Angle",
    ylab = "Declustered rate",
    main = "Declustered rate by angle"
  )
  abline(h = x$global_rate, lty = 2, lwd = 2)

  plot(
    x$phi_grid,
    x$density,
    type = "l",
    lwd = 2,
    xlab = "Angle",
    ylab = "Angular density",
    main = "Normalized angular density"
  )

  invisible(x)
}
