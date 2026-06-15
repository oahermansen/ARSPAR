spar_apply_schema_names <- function(x, data, space = c("original", "transformed", "angular")) {

  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }

  space <- match.arg(space)

  if (is.null(data)) {
    return(NULL)
  }

  if (!is.matrix(data) && !is.data.frame(data)) {
    data <- as.matrix(data)
  }

  names_expected <- spar_schema_names(x, space)

  if (!is.null(names_expected)) {

    if (length(names_expected) != ncol(data)) {
      stop(
        sprintf(
          "Schema names for '%s' space have length %d but data has %d columns.",
          space,
          length(names_expected),
          ncol(data)
        ),
        call. = FALSE
      )
    }

    colnames(data) <- names_expected
  }

  data
}

spar_schema_meta_names <- function(x) {

  sch <- x$schema

  list(
    observation_id = sch$observation_id_name %||% "observation_id",
    time = sch$time_name %||% "time"
  )
}
spar_slice_obs <- function(x, i) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.matrix(x) || is.data.frame(x)) {
    return(x[i, , drop = FALSE])
  }

  if (is.atomic(x)) {
    return(x[i])
  }

  if (is.list(x)) {
    return(lapply(x, function(v) {
      if (is.null(v)) {
        NULL
      } else if (is.matrix(v) || is.data.frame(v)) {
        v[i, , drop = FALSE]
      } else {
        v[i]
      }
    }))
  }

  x
}
