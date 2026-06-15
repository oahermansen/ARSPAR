spar_pick_gpd_backend <- function(preferred = c("ismev", "POT", "evir"), require_available = FALSE) {
  preferred <- as.character(preferred)
  preferred <- preferred[nzchar(preferred)]
  if (length(preferred) == 0L) {
    stop("`preferred` must contain at least one backend name.", call. = FALSE)
  }

  for (nm in preferred) {
    if (requireNamespace(nm, quietly = TRUE)) {
      return(list(name = nm, available = TRUE, note = NULL))
    }
  }

  if (isTRUE(require_available)) {
    stop(sprintf("Install one of these packages for GPD fits: %s", paste(preferred, collapse = ", ")), call. = FALSE)
  }

  list(name = "none", available = FALSE, note = sprintf("Install one of: %s", paste(preferred, collapse = ", ")))
}

spar_first_non_null <- function(...) {
  xs <- list(...)
  for (x in xs) {
    if (!is.null(x)) return(x)
  }
  NULL
}

spar_extract_gpd_scale_shape <- function(par) {
  if (is.null(par)) return(c(scale = NA_real_, shape = NA_real_))
  par <- as.numeric(par)
  nm <- names(par)

  if (!is.null(nm)) {
    i_scale <- grep("scale|sigma|beta", nm, ignore.case = TRUE)
    i_shape <- grep("shape|xi", nm, ignore.case = TRUE)
    if (length(i_scale) >= 1L && length(i_shape) >= 1L) {
      return(c(scale = par[i_scale[1]], shape = par[i_shape[1]]))
    }
  }

  if (length(par) >= 2L) c(scale = par[1], shape = par[2]) else c(scale = NA_real_, shape = NA_real_)
}

spar_fit_gpd_constant <- function(y, min_n = 30L, backend = NULL) {
  y <- as.numeric(y)
  y <- y[is.finite(y) & y > 0]

  if (is.null(backend)) {
    backend <- spar_pick_gpd_backend(require_available = FALSE)$name
  }

  out <- list(
    n = length(y),
    scale = NA_real_,
    shape = NA_real_,
    fit = NULL,
    converged = FALSE,
    engine = backend
  )

  if (length(y) < min_n) return(out)
  if (!backend %in% c("ismev", "POT", "evir")) return(out)

  fit <- tryCatch(
    {
      if (identical(backend, "ismev")) {
        ismev::gpd.fit(y, threshold = 0, show = FALSE)
      } else if (identical(backend, "POT")) {
        POT::fitgpd(y, threshold = 0, est = "mle")
      } else {
        evir::gpd(y, threshold = 0)
      }
    },
    error = function(e) {
      warning(sprintf("%s GPD fit failed for sample of size %d: %s", backend, length(y), conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) return(out)

  par <- spar_first_non_null(fit$mle, fit$param, fit$par.ests, fit$estimate)
  ss <- spar_extract_gpd_scale_shape(par)

  out$scale <- as.numeric(ss["scale"])
  out$shape <- as.numeric(ss["shape"])
  out$fit <- fit
  out$converged <- is.finite(out$scale) && is.finite(out$shape)
  out
}

spar_gp_survival <- function(y, sigma, xi) {
  y <- as.numeric(y)
  out <- rep(NA_real_, length(y))

  ok <- is.finite(y) & y >= 0 & is.finite(sigma) & sigma > 0 & is.finite(xi)
  if (!any(ok)) return(out)

  if (abs(xi) < 1e-10) {
    out[ok] <- exp(-y[ok] / sigma)
    return(out)
  }

  z <- 1 + xi * y[ok] / sigma
  out[ok] <- 0
  ok2 <- z > 0
  out_idx <- which(ok)
  out[out_idx[ok2]] <- z[ok2]^(-1 / xi)
  out
}
