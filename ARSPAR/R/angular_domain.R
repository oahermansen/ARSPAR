#' Create an angular domain specification
#'
#' Constructs a `spar_angle_domain` object describing how angular coordinates
#' should be interpreted.
#'
#' An angular domain specifies the support and geometry of angle values, such as:
#' \itemize{
#'   \item whether the domain is cyclical or interval-based,
#'   \item the lower and upper bounds of the support,
#'   \item whether the support is closed on the left and/or right,
#'   \item the measurement units,
#'   \item and an optional label for display or bookkeeping.
#' }
#'
#' This object is intended to be stored in the `angular$domain` slot of a
#' `spar_representation` once an angular-radial representation has been defined.
#' It provides the contextual information needed for domain-aware operations such
#' as normalization, ordering, wrapped differences, and angular distances.
#'
#' @param type Character string specifying the angular geometry. Must be one of
#'   `"cyclical"` or `"interval"`.
#' @param lower Single numeric value giving the lower bound of the angular
#'   support.
#' @param upper Single numeric value giving the upper bound of the angular
#'   support. Must be strictly greater than `lower`.
#' @param closed_left Logical; whether the support is closed on the left.
#' @param closed_right Logical; whether the support is closed on the right.
#' @param units Character string describing the measurement units. Must be one of
#'   `"radians"` or `"degrees"`.
#' @param label Optional character label for the domain.
#'
#' @details
#' A cyclical domain indicates that angles should be interpreted modulo the
#' domain width, typically using wrapped notions of difference and distance.
#'
#' An interval domain indicates that angles should be interpreted on a linear
#' interval without wrapping.
#'
#' The constructor stores the derived width `upper - lower` for convenience.
#'
#' @return An object of class `"spar_angle_domain"`.
#'
#' @examples
#' dom1 <- new_spar_angle_domain()
#'
#' dom2 <- new_spar_angle_domain(
#'   type = "interval",
#'   lower = -pi / 2,
#'   upper = pi / 2,
#'   label = "upper-half support"
#' )
#'
#' dom1
#' dom2
#'
#' @seealso [print.spar_angle_domain()]
#'
#' @export
new_spar_angle_domain <- function(
    type = c("cyclical", "interval"),
    lower = 0,
    upper = 2 * pi,
    closed_left = TRUE,
    closed_right = FALSE,
    units = c("radians", "degrees"),
    label = NULL
) {
  type <- match.arg(type)
  units <- match.arg(units)
  
  if (!is.numeric(lower) || length(lower) != 1L || is.na(lower)) {
    stop("`lower` must be a single non-missing numeric value.", call. = FALSE)
  }
  
  if (!is.numeric(upper) || length(upper) != 1L || is.na(upper)) {
    stop("`upper` must be a single non-missing numeric value.", call. = FALSE)
  }
  
  if (upper <= lower) {
    stop("`upper` must be strictly greater than `lower`.", call. = FALSE)
  }
  
  if (!is.logical(closed_left) || length(closed_left) != 1L || is.na(closed_left)) {
    stop("`closed_left` must be a single non-missing logical value.", call. = FALSE)
  }
  
  if (!is.logical(closed_right) || length(closed_right) != 1L || is.na(closed_right)) {
    stop("`closed_right` must be a single non-missing logical value.", call. = FALSE)
  }
  
  if (!is.null(label)) {
    if (!is.character(label) || length(label) != 1L || is.na(label)) {
      stop("`label` must be `NULL` or a single non-missing character string.", call. = FALSE)
    }
  }
  
  structure(
    list(
      type = type,
      lower = lower,
      upper = upper,
      width = upper - lower,
      closed_left = closed_left,
      closed_right = closed_right,
      units = units,
      label = label
    ),
    class = "spar_angle_domain"
  )
}

validate_spar_angle_domain <- function(x) {
  if (!inherits(x, "spar_angle_domain")) {
    stop("`x` must be a `spar_angle_domain` object.", call. = FALSE)
  }

  if (!is.character(x$type) || length(x$type) != 1L || !x$type %in% c("cyclical", "interval")) {
    stop("Invalid angle domain type.", call. = FALSE)
  }

  if (!is.numeric(x$lower) || length(x$lower) != 1L || !is.finite(x$lower)) {
    stop("Invalid angle domain lower bound.", call. = FALSE)
  }

  if (!is.numeric(x$upper) || length(x$upper) != 1L || !is.finite(x$upper)) {
    stop("Invalid angle domain upper bound.", call. = FALSE)
  }

  if (!is.numeric(x$width) || length(x$width) != 1L || !is.finite(x$width) || x$width <= 0) {
    stop("Invalid angle domain width.", call. = FALSE)
  }

  if (!identical(unname(x$upper - x$lower), unname(x$width))) {
    stop("Angle domain width does not match upper-lower.", call. = FALSE)
  }

  if (!is.logical(x$closed_left) || length(x$closed_left) != 1L || is.na(x$closed_left)) {
    stop("Invalid angle domain left-closure flag.", call. = FALSE)
  }

  if (!is.logical(x$closed_right) || length(x$closed_right) != 1L || is.na(x$closed_right)) {
    stop("Invalid angle domain right-closure flag.", call. = FALSE)
  }

  invisible(x)
}

as_spar_angle_domain <- function(domain) {
  if (inherits(domain, "spar_angle_domain")) {
    return(validate_spar_angle_domain(domain))
  }

  if (inherits(domain, "angle_domain")) {
    type <- if (identical(domain$type, "circular")) "cyclical" else "interval"

    if (identical(type, "cyclical")) {
      lower <- domain$min
      upper <- domain$max
      if ((!is.finite(lower) || !is.finite(upper)) && is.finite(domain$period)) {
        lower <- 0
        upper <- domain$period
      }
    } else {
      lower <- domain$min
      upper <- domain$max
    }

    return(new_spar_angle_domain(
      type = type,
      lower = lower,
      upper = upper,
      closed_left = TRUE,
      closed_right = FALSE,
      units = "radians",
      label = domain$label %||% NULL
    ))
  }

  stop("`domain` must be an angle-domain object.", call. = FALSE)
}

spar_normalize_angle <- function(phi, domain) {
  domain <- as_spar_angle_domain(domain)
  phi <- as.numeric(phi)

  if (domain$type == "interval") {
    return(phi)
  }

  ((phi - domain$lower) %% domain$width) + domain$lower
}

spar_angle_delta <- function(phi, phi0, domain) {
  domain <- as_spar_angle_domain(domain)

  phi <- as.numeric(phi)
  phi0 <- as.numeric(phi0)

  if (domain$type == "interval") {
    return(phi - phi0)
  }

  ((phi - phi0 + domain$width / 2) %% domain$width) - domain$width / 2
}

spar_angle_distance <- function(phi, phi0, domain) {
  abs(spar_angle_delta(phi, phi0, domain))
}

spar_angle_gt <- function(phi, phi0, domain) {
  spar_angle_delta(phi, phi0, domain) > 0
}

spar_angle_lt <- function(phi, phi0, domain) {
  spar_angle_delta(phi, phi0, domain) < 0
}

spar_angle_positive_distance <- function(from, to, domain){
  if(domain$type == "cyclical"){
    spar_angle_delta(to, from, domain) %% 4
  }
}

spar_angle_in_range <- function(phi, range, domain, closed = TRUE) {
  domain <- as_spar_angle_domain(domain)

  if (!is.numeric(range) || length(range) != 2L || anyNA(range)) {
    stop("`range` must be a numeric length-2 vector.", call. = FALSE)
  }

  phi <- as.numeric(phi)
  start <- as.numeric(range[1])
  end <- as.numeric(range[2])

  if (domain$type == "interval") {
    if (isTRUE(closed)) {
      return(phi >= start & phi <= end)
    }
    return(phi > start & phi < end)
  }
  
  dphi <- spar_angle_positive_distance(phi, end, domain)
  dstart <- spar_angle_positive_distance(start, end, domain)
  if(isTRUE(closed)){
    return(dphi <= dstart)
  }
  return(dphi < dstart)
  
  d_start <- spar_angle_delta(phi, start, domain)
  d_end <- spar_angle_delta(phi, end, domain)

  if (isTRUE(closed)) {
    if (start <= end) {
      return(d_start >= 0 & d_end <= 0)
    }
    return(d_start >= 0 | d_end <= 0)
  }

  if (start < end) {
    return(d_start > 0 & d_end < 0)
  }
  d_start > 0 | d_end < 0
}

spar_angle_range_from_transitions <- function(phi, domain){
  domain <- as_spar_angle_domain(domain)
  if(domain$type != "cyclical"){
    stop("Only defined for cyclical domains")
  }
  phi <- as.numeric(phi)
  phi <- phi[is.finite(phi)]
  if (length(phi) == 0L) {
    return(c(NA_real_, NA_real_))
  }
  if (length(phi) == 1L) {
    p <- spar_normalize_angle(phi, domain)
    return(c(p, p))
  }

  pair_ranges <- lapply(seq_len(length(phi) - 1L), function(i) {
    spar_angle_range(phi[c(i, i + 1L)], domain)
  })

  return(spar_angle_merge_overlapping_range_list(pair_ranges, domain))
  if (length(merged) == 0L) {
    return(c(NA_real_, NA_real_))
  }
  if (length(merged) == 1L) {
    return(as.numeric(merged[[1]]))
  }

  spans <- vapply(merged, function(r) spar_angle_span(r, domain), numeric(1))
  as.numeric(merged[[which.min(spans)]])
}

spar_angle_merge_overlapping_range_list <- function(ranges, domain){
  res <- c()
  for(i in 1:length(ranges)){
    if(i == 1){
      res <- ranges[[i]]
    }else{
      merged <- spar_angle_merge_overlapping_ranges(res, ranges[[i]], domain)
      res <- merged
    }
  }
  return(res)
}

spar_angle_merge_range_list <- function(ranges, domain, inner = FALSE){
  domain <- as_spar_angle_domain(domain)
  if (length(ranges) == 0L) {
    return(list())
  }
  normalize_range <- function(r) {
    if (!is.numeric(r) || length(r) != 2L || anyNA(r) || any(!is.finite(r))) return(NULL)
    as.numeric(r)
  }

  rs <- lapply(ranges, normalize_range)
  rs <- rs[!vapply(rs, is.null, logical(1))]
  if (length(rs) == 0L) return(list())

  overlaps_or_contains <- function(a, b) {
    if (!is.numeric(a) || !is.numeric(b) || length(a) != 2L || length(b) != 2L || anyNA(a) || anyNA(b)) {
      return(FALSE)
    }
    any(spar_angle_in_range(a, b, domain)) ||
      any(spar_angle_in_range(b, a, domain))
  }

  changed <- TRUE
  while (changed) {
    changed <- FALSE
    out <- list()
    for (r in rs) {
      merged_r <- r
      keep_out <- list()
      did_merge <- FALSE
      for (k in seq_along(out)) {
        o <- out[[k]]
        if (overlaps_or_contains(o, merged_r)) {
          merged_r <- spar_angle_merge_overlapping_ranges(o, merged_r, domain)
          did_merge <- TRUE
          changed <- TRUE
        } else {
          keep_out[[length(keep_out) + 1L]] <- o
        }
      }
      keep_out[[length(keep_out) + 1L]] <- merged_r
      out <- keep_out
      if (did_merge) changed <- TRUE
    }
    rs <- out
  }

  rs
}

spar_angle_merge_overlapping_ranges <- function(range1, range2, domain){
  domain <- as_spar_angle_domain(domain)
  if(domain$type != "cyclical"){
    stop("Only defined for cyclical domains")
  }
  if(!is.numeric(range1) || !is.numeric(range2) || length(range1) != 2L || length(range2) != 2L || anyNA(range1) || anyNA(range2)){
    stop("Only defined for range arguments.")
  }
  check1 <- which(!spar_angle_in_range(range1, range2, domain))
  check2 <- which(!spar_angle_in_range(range2, range1, domain))
  if(length(check1) > 1 && length(check2) > 1){
    stop("Ranges are not overlapping")
  }
  if(length(check1) > 1 || length(check2) > 1){
    out <- if (length(check1) > 1) range1 else range2
    return(as.numeric(out))
  }
  
  outside <- range1[check1]
  dist <- spar_angle_distance(outside, range1[-check1], domain)
  dc <- c(spar_angle_delta(outside, range2[1], domain), spar_angle_delta(outside, range2[2], domain))
  distcheck <- dist >= abs(dc)
  if(!all(distcheck)){
    range2[distcheck] <- outside
  }else{
    ind <- which.min(abs(dc))
    range2[ind] <- outside
  }
  as.numeric(range2)
}

spar_angle_range <- function(phi, domain) {
  domain <- as_spar_angle_domain(domain)
  phi <- phi[is.finite(phi)]
  n <- length(phi)

  if (n == 0L) {
    return(c(NA_real_, NA_real_))
  }

  if (domain$type == "interval") {
    return(c(min(phi), max(phi)))
  }

  phi <- sort(spar_normalize_angle(phi, domain))
  if (n == 1L) {
    return(c(phi, phi))
  }

  ext <- c(phi, phi[1] + domain$width)
  gaps <- diff(ext)
  i_gap <- which.max(gaps)

  start <- phi[(i_gap %% n) + 1L]
  end <- phi[i_gap]

  c(start, end)
}

spar_angle_span <- function(phi, domain) {
  domain <- as_spar_angle_domain(domain)
  r <- spar_angle_range(phi, domain)

  if (anyNA(r)) {
    return(NA_real_)
  }

  if (domain$type == "interval") {
    return(r[2] - r[1])
  }

  if (r[2] >= r[1]) {
    return(r[2] - r[1])
  }
  domain$width - (r[1] - r[2])
}

#' Compute shortest directed angular span between two angles
#'
#' Returns the shortest directed arc from `phi_from` to `phi_to` under the
#' supplied domain. For cyclical domains, the arc may cross the seam and is then
#' returned as two contiguous segments in plotting coordinates.
#'
#' @param phi_from Start angle.
#' @param phi_to End angle.
#' @param domain Angular domain object.
#'
#' @return A list containing normalized endpoints, directed delta, seam-crossing
#'   flag, and a `segments` data frame with columns `start` and `end`.
#'
#' @export
spar_smallest_angular_span <- function(phi_from, phi_to, domain) {
  domain <- as_spar_angle_domain(domain)

  phi_from <- as.numeric(phi_from)
  phi_to <- as.numeric(phi_to)

  if (length(phi_from) != 1L || length(phi_to) != 1L || !is.finite(phi_from) || !is.finite(phi_to)) {
    stop("`phi_from` and `phi_to` must be single finite numeric values.", call. = FALSE)
  }

  if (domain$type == "interval") {
    return(list(
      from = phi_from,
      to = phi_to,
      delta = phi_to - phi_from,
      crosses_seam = FALSE,
      segments = data.frame(start = phi_from, end = phi_to)
    ))
  }

  w <- domain$width
  lo <- domain$lower
  hi <- domain$upper

  from_n <- spar_normalize_angle(phi_from, domain)
  delta <- spar_angle_delta(phi_to, from_n, domain)
  to_u <- from_n + delta

  if (!is.finite(delta) || abs(delta) < .Machine$double.eps^0.5) {
    return(list(
      from = from_n,
      to = from_n,
      delta = 0,
      crosses_seam = FALSE,
      segments = data.frame(start = from_n, end = from_n)
    ))
  }

  if (to_u >= lo && to_u <= hi) {
    return(list(
      from = from_n,
      to = to_u,
      delta = delta,
      crosses_seam = FALSE,
      segments = data.frame(start = from_n, end = to_u)
    ))
  }

  if (to_u > hi) {
    seg <- data.frame(start = c(from_n, lo), end = c(hi, to_u - w))
  } else {
    seg <- data.frame(start = c(from_n, hi), end = c(lo, to_u + w))
  }

  list(
    from = from_n,
    to = spar_normalize_angle(phi_to, domain),
    delta = delta,
    crosses_seam = TRUE,
    segments = seg
  )
}

spar_angle_grid <- function(n = 100L, domain, respect_closure = FALSE) {
  domain <- as_spar_angle_domain(domain)
  if (!is.numeric(n) || length(n) != 1L || n < 1L || is.na(n)) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!respect_closure) {
    if (n == 1L) {
        return(domain$lower)
    }
    return(seq(domain$lower, domain$upper, length.out = n))
  }else {
    o_left <- as.numeric(!domain$closed_left)
    o_right <- as.numeric(!domain$closed_right)
    
    if (n == 1L) {
      if (domain$closed_left) {
        return(domain$lower)
      } else {
        return(domain$lower + domain$width/2)
      }
    }
    n <- n + (!o_left + !o_right)
    out <- seq(domain$lower, domain$upper, length.out = n)
    lower_ind <- 1 + o_left
    upper_ind <- n - o_right
    
    return(out[lower_ind:upper_ind])
  }
}

#' Print an angular domain specification
#'
#' Prints a compact summary of a `spar_angle_domain` object.
#'
#' @param x A `spar_angle_domain` object.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#'
#' @export
print.spar_angle_domain <- function(x, ...) {
  validate_spar_angle_domain(x)
  cat("<spar_angle_domain>\n")
  cat("  type:  ", x$type, "\n", sep = "")
  cat("  range: ", x$lower, " to ", x$upper, "\n", sep = "")
  cat("  width: ", x$width, "\n", sep = "")
  cat("  units: ", x$units, "\n", sep = "")
  cat(
    "  closure: ",
    if (isTRUE(x$closed_left)) "[" else "(",
    "lower, upper",
    if (isTRUE(x$closed_right)) "]" else ")",
    "\n",
    sep = ""
  )
  if (!is.null(x$label)) {
    cat("  label: ", x$label, "\n", sep = "")
  }
  invisible(x)
}
