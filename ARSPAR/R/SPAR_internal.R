#' Build the internal active observation matrix
#'
#' Constructs the internal observation-aligned matrix used as the basis for
#' active SPAR data views.
#'
#' This function collects available observation-level components from a
#' `spar_representation` object and column-binds them into a single matrix.
#' The resulting matrix uses **internal field names** rather than schema
#' names and is intended strictly for internal use.
#'
#' The internal naming convention is stable and algorithm-oriented. Examples
#' include:
#' \itemize{
#'   \item original-space variables
#'   \item transformed variables
#'   \item `"R"` for the radial coordinate
#'   \item `"phi"` for the angular coordinate
#'   \item `"threshold"` for per-observation thresholds
#'   \item `"excess"` for excess values
#' }
#'
#' This matrix is the foundation used by [spar_internal_view()] and
#' [spar_active_view()].
#'
#' @param x A `spar_representation` object.
#' @param include_time Logical; if `TRUE`, include the time variable in the
#'   resulting matrix.
#' @param include_observation_id Logical; if `TRUE`, include the observation ID
#'   column in the resulting matrix.
#'
#' @return A numeric matrix where each row corresponds to an observation and
#' each column corresponds to an internally named field.
#'
#' @details
#' The current implementation builds the matrix by column-binding components
#' stored in the representation object. In the future this function may instead
#' return a matrix backed by an optimized workspace layout.
#'
#' @seealso [spar_internal_view()], [spar_active_view()]
#'
#' @keywords internal
spar_build_active_matrix <- function(
    x,
    include_time = FALSE,
    include_observation_id = FALSE
) {
  
  parts <- list()
  
  # original
  parts[[length(parts) + 1L]] <- x$data$X_original
  
  # transformed
  if (!is.null(x$transform$transformed_data)) {
    parts[[length(parts) + 1L]] <- x$transform$transformed_data
  }
  
  # angular
  if (!is.null(x$angular$R) && !is.null(x$angular$phi)) {
    parts[[length(parts) + 1L]] <- cbind(
      R = x$angular$R,
      phi = x$angular$phi
    )
  }
  
  # threshold
  if (!is.null(x$threshold$per_observation)) {
    parts[[length(parts) + 1L]] <- matrix(
      x$threshold$per_observation,
      ncol = 1,
      dimnames = list(NULL, "threshold")
    )
  }
  
  # excess
  if (!is.null(x$excess$value)) {
    parts[[length(parts) + 1L]] <- matrix(
      x$excess$value,
      ncol = 1,
      dimnames = list(NULL, "excess")
    )
  }
  
  mat <- do.call(cbind, parts)
  
  mat
}
#' Construct the internal observation view
#'
#' Creates a `spar_obs_view` object representing the **internal observation
#' view** of a `spar_representation`.
#'
#' The internal view exposes a matrix-like interface to observation-level data
#' using **internal field names**. This interface is intended for algorithms,
#' diagnostics, and other internal computations where stable naming is required.
#'
#' Internally, this function uses [spar_build_active_matrix()] to assemble the
#' underlying observation matrix.
#'
#' @param x A `spar_representation` object.
#' @param include_time Logical; if `TRUE`, include the observation time column
#'   in the view.
#' @param include_observation_id Logical; if `TRUE`, include the observation ID
#'   column in the view.
#'
#' @return A `spar_obs_view` object containing an internally named observation
#' matrix in the `data` field.
#'
#' @details
#' The internal view is intended for use by internal SPAR algorithms. Column
#' names follow the internal naming convention rather than the user-facing
#' schema.
#'
#' For example, angular coordinates will appear as `"R"` and `"phi"` even if
#' the user-facing schema defines different names.
#'
#' This function is typically used as the first step of internal analysis
#' routines:
#'
#' \preformatted{
#' V <- spar_internal_view(x)
#' M <- as.matrix(V)
#' }
#'
#' @seealso [spar_build_active_matrix()], [spar_active_view()]
#'
#' @keywords internal
spar_internal_view <- function(
    x,
    include_time = FALSE,
    include_observation_id = FALSE
) {
  
  mat <- spar_build_active_matrix(
    x,
    include_time = include_time,
    include_observation_id = include_observation_id
  )
  
  structure(
    list(data = mat),
    class = "spar_obs_view"
  )
  
}

#' Construct the user-facing active observation view
#'
#' Creates a `spar_obs_view` object representing the **active observation view**
#' of a `spar_representation`.
#'
#' This view provides a matrix-like interface to observation-level data using
#' **schema-defined names**, making it suitable for interactive use and
#' diagnostics.
#'
#' Internally, the active view is constructed by first building the internal
#' observation matrix via [spar_internal_view()] and then projecting column
#' names through the representation's schema.
#'
#' @param x A `spar_representation` object.
#' @param include_time Logical; if `TRUE`, include the observation time column
#'   in the view.
#' @param include_observation_id Logical; if `TRUE`, include the observation ID
#'   column in the view.
#'
#' @return A `spar_obs_view` object containing the active observation matrix
#' with **schema-projected column names**.
#'
#' @details
#' The active observation view is intended for user-facing data exploration and
#' diagnostic workflows. Column names follow the naming schema stored in the
#' representation object.
#'
#' For example, if the schema specifies custom names for angular coordinates,
#' the active view will use those names rather than the internal `"R"` and
#' `"phi"` identifiers.
#'
#' The semantics of the returned object are matrix-like:
#'
#' \itemize{
#'   \item rows correspond to observations
#'   \item columns correspond to available observation-level fields
#' }
#'
#' Subsetting operations on the resulting object are supported via the
#' `[.spar_obs_view` method.
#'
#' @examples
#' # Assuming `obj` is a spar_representation
#' # V <- spar_active_view(obj)
#' # as.matrix(V)
#'
#' @seealso [spar_internal_view()], [spar_build_active_matrix()]
#'
#' @export
spar_active_view <- function(
    x,
    include_time = FALSE,
    include_observation_id = FALSE
) {
  
  view <- spar_internal_view(
    x,
    include_time = include_time,
    include_observation_id = include_observation_id
  )
  
  mat <- view$data
  
  # apply schema projection
  colnames(mat) <- spar_project_schema_names(
    x,
    colnames(mat)
  )
  
  view$data <- mat
  
  view
}

#' Coerce a spar_obs_view to a matrix
#'
#' Returns the underlying joined observation matrix stored in the
#' `data` field of a `spar_obs_view`.
#'
#' @param x A `spar_obs_view` object.
#' @param ... Unused.
#'
#' @return A numeric matrix containing the joined observation data.
#'
#' @keywords internal
#' @export
as.matrix.spar_obs_view <- function(x, ...) {
  x$data
}

#' Get dimensions of a spar_obs_view.
#'
#' Returns the dimensions of the underlying joined observation matrix stored in the
#' `data` field of a `spar_obs_view`.
#'
#' @param x A `spar_obs_view` object.
#' @param ... Unused.
#'
#' @return A numeric vector containing the dimension sizes.
#'
#' @keywords internal
#' @export
dim.spar_obs_view <- function(x) {
  dim(x$data)
}

#' Dispatch colnames on the underlying matrix
#' 
#' @param x A `spar_obs_view` object.
#' @param ... Unused.
#'
#' @return A collection of column names.
#'
#' @keywords internal
#' @export
colnames.spar_obs_view <- function(x) {
  colnames(x$data)
}

#' Dispatch rownames on the underlying matrix
#' 
#' @param x A `spar_obs_view` object.
#' @param ... Unused.
#'
#' @return A collection of row names.
#'
#' @keywords internal
#' @export
rownames.spar_obs_view <- function(x) {
  rownames(x$data)
}

#' Print information about the spar_obs_view.
#' 
#' @param x A `spar_obs_view` object.
#' @param ... Unused.
#'
#'
#' @keywords internal
#' @export
print.spar_obs_view <- function(x, ...) {
  cat("<spar_obs_view>\n")
  print(utils::head(x$data, 10L))
  if (nrow(x$data) > 10) {
    cat("... with", nrow(x$data) - 10, "more rows\n")
  }
  invisible(x)
}

#' Process a row selector
#'
#' S3 generic for normalizing row selectors used when subsetting active SPAR
#' observation views.
#'
#' @param x Context object, typically a matrix-like object supplying the row
#'   dimension.
#' @param i Row selector.
#' @param n Number of rows available.
#' @param ... Additional arguments for methods.
#'
#' @return An integer vector of row positions.
#'
#' @keywords internal
#' @noRd
spar_process_rows <- function(x, i, n = nrow(x), ...) {
  UseMethod("spar_process_rows", i)
}


#' Process a column selector
#'
#' S3 generic for normalizing column selectors used when subsetting active SPAR
#' observation views.
#'
#' @param x Context object, typically a matrix-like object supplying the column
#'   names.
#' @param j Column selector.
#' @param names Column names available for matching.
#' @param ... Additional arguments for methods.
#'
#' @return An integer vector of column positions.
#'
#' @keywords internal
#' @noRd
spar_process_cols <- function(x, j, names = colnames(x), ...) {
  UseMethod("spar_process_cols", j)
}

#' Resolve a row selector to integer positions
#'
#' Normalizes a row selector into integer row positions for a matrix-like active
#' observation view.
#'
#' This function is a thin wrapper around [spar_process_rows()] that handles
#' missing and `NULL` selectors by returning all available row positions.
#'
#' It is intended as the central entry point for row-index interpretation in
#' internal SPAR subsetting code.
#'
#' @param x A matrix-like object providing the row dimension.
#' @param i Row selector. May be missing, `NULL`, numeric, or logical. Future
#'   custom row-index classes may be supported through
#'   [spar_process_rows()].
#' @param n Number of rows available. Defaults to `nrow(x)`.
#'
#' @return An integer vector of row positions.
#'
#' @details
#' This function exists so that the interpretation of row selectors is
#' centralized in one place. It is intended to support future extensions such
#' as time-index classes, observation-index classes, cluster-index classes, or
#' span-based row selectors.
#'
#' @examples
#' m <- matrix(1:12, nrow = 4)
#' spar_resolve_rows(m, NULL)
#' spar_resolve_rows(m, 1:2)
#' spar_resolve_rows(m, c(TRUE, FALSE, TRUE, FALSE))
#'
#' @seealso [spar_process_rows()], [spar_resolve_cols()]
#'
#' @keywords internal
#' @export
spar_resolve_rows <- function(x, i, n = nrow(x)) {
  if (missing(i) || is.null(i)) {
    return(seq_len(n))
  }
  
  spar_process_rows(x, i, n = n)
}

#' Resolve a column selector to integer positions
#'
#' Normalizes a column selector into integer column positions for a matrix-like
#' active observation view.
#'
#' This function is a thin wrapper around [spar_process_cols()] that handles
#' missing and `NULL` selectors by returning all available column positions.
#'
#' It is intended as the central entry point for column-index interpretation in
#' internal SPAR subsetting code.
#'
#' @param x A matrix-like object providing column names.
#' @param j Column selector. May be missing, `NULL`, numeric, logical, or
#'   character. Future custom column-index classes may be supported through
#'   [spar_process_cols()].
#' @param names Column names available for matching. Defaults to `colnames(x)`.
#'
#' @return An integer vector of column positions.
#'
#' @details
#' Character column selectors are matched against the column names of the active
#' observation view.
#'
#' This function exists so that the interpretation of column selectors is
#' centralized in one place. It is intended to support future extensions such
#' as layout-index classes, span selectors, grouped field selectors, or other
#' custom column-selection mechanisms.
#'
#' @examples
#' m <- matrix(1:12, nrow = 4)
#' colnames(m) <- c("a", "b", "c")
#'
#' spar_resolve_cols(m, NULL)
#' spar_resolve_cols(m, 1:2)
#' spar_resolve_cols(m, c(TRUE, FALSE, TRUE))
#' spar_resolve_cols(m, c("a", "c"))
#'
#' @seealso [spar_process_cols()], [spar_resolve_rows()]
#'
#' @keywords internal
#' @export
spar_resolve_cols <- function(x, j, names = colnames(x)) {
  if (missing(j) || is.null(j)) {
    return(seq_along(names))
  }
  
  spar_process_cols(x, j, names = names)
}

#' Subset the active observation view of a SPAR representation
#'
#' Subsets the active joined observation-aligned view of a
#' `spar_representation`.
#'
#' This method does **not** subset the nested internal list structure of the
#' representation object itself. Instead, it constructs an observation-aligned
#' active view and then applies row and column indexing to that view.
#'
#' By default, the active view uses the user-facing naming schema. This keeps
#' subsetting convenient for interactive use:
#'
#' \preformatted{
#' x[1:10, ]
#' x[, c("R", "phi")]
#' }
#'
#' Internal code can request the internal naming convention by setting
#' `internal = TRUE`:
#'
#' \preformatted{
#' x[1:10, c("R", "phi", "threshold"), internal = TRUE]
#' }
#'
#' Row and column selectors are normalized through [spar_resolve_rows()] and
#' [spar_resolve_cols()]. This provides a central extension point for future
#' custom index classes such as time indices, observation indices, layout
#' indices, span indices, or cluster indices.
#'
#' @param x A `spar_representation` object.
#' @param i Row selector. May be missing, `NULL`, numeric, or logical. Future
#'   index classes may also be supported through the row-index processing
#'   system.
#' @param j Column selector. May be missing, `NULL`, numeric, logical, or
#'   character. Character selectors are matched against the column names of the
#'   chosen active view.
#' @param drop Logical; passed through to matrix subsetting. If `TRUE`, a
#'   dropped result may be returned as a bare vector or scalar. If `FALSE`, the
#'   result is returned as a `spar_obs_view`.
#' @param internal Logical; if `TRUE`, subset the internal observation view
#'   using stable internal field names. If `FALSE` (default), subset the
#'   user-facing active view using schema-projected names.
#'
#' @return
#' If `drop = FALSE`, a `spar_obs_view` object representing the selected
#' portion of the chosen active observation view.
#'
#' If `drop = TRUE`, the dropped result of subsetting the underlying active
#' view matrix.
#'
#' @details
#' The semantics of `x[i, j]` for a `spar_representation` are intentionally
#' matrix-like:
#'
#' \itemize{
#'   \item rows correspond to observations,
#'   \item columns correspond to fields in the active joined observation view.
#' }
#'
#' The exact contents of the active observation view are determined by
#' [spar_active_view()] or [spar_internal_view()], depending on the value of
#' `internal`.
#'
#' This method is designed so that the underlying implementation of the active
#' view can later be optimized, for example by using a packed contiguous
#' workspace layout, without changing the external subsetting interface.
#'
#' @examples
#' # Assuming `obj` is a spar_representation and an active view can be built:
#' # obj[1:10, ]
#' # obj[, c("R", "phi")]
#' # obj[, c("R", "phi", "threshold"), internal = TRUE]
#' # obj[1:5, 1:3, drop = TRUE]
#'
#' @seealso [spar_active_view()], [spar_internal_view()],
#'   [spar_resolve_rows()], [spar_resolve_cols()]
#'
#' @export
`[.spar_representation` <- function(x, i, j, drop = FALSE, internal = FALSE) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  
  view <- if (isTRUE(internal)) {
    spar_internal_view(x)
  } else {
    spar_active_view(x)
  }
  
  if (!inherits(view, "spar_obs_view")) {
    stop(
      "`spar_active_view(x)` must return a `spar_obs_view` object.",
      call. = FALSE
    )
  }
  
  mat <- view$data
  
  if (!is.matrix(mat)) {
    stop(
      "The active observation view must contain a matrix in `view$data`.",
      call. = FALSE
    )
  }
  
  ii <- if (missing(i)) {
    spar_resolve_rows(mat, NULL, n = nrow(mat))
  } else {
    spar_resolve_rows(mat, i, n = nrow(mat))
  }
  
  jj <- if (missing(j)) {
    spar_resolve_cols(mat, NULL, names = colnames(mat))
  } else {
    spar_resolve_cols(mat, j, names = colnames(mat))
  }
  
  out <- mat[ii, jj, drop = drop]
  
  if (isTRUE(drop)) {
    return(out)
  }
  
  structure(
    list(
      data = out
    ),
    class = "spar_obs_view"
  )
}


#' @keywords internal
#' @noRd
spar_process_cols.default <- function(x, j, names = colnames(x), ...) {
  stop(
    sprintf(
      "Unsupported column index type: %s",
      paste(class(j), collapse = "/")
    ),
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
spar_process_cols.NULL <- function(x, j, names = colnames(x), ...) {
  seq_along(names)
}

#' @keywords internal
#' @noRd
spar_process_cols.numeric <- function(x, j, names = colnames(x), ...) {
  j <- as.integer(j)
  p <- length(names)
  
  if (anyNA(j)) {
    stop("Column index contains NA values.", call. = FALSE)
  }
  
  if (any(j < 1L | j > p)) {
    stop("Column index out of bounds.", call. = FALSE)
  }
  
  j
}

#' @keywords internal
#' @noRd
spar_process_cols.logical <- function(x, j, names = colnames(x), ...) {
  p <- length(names)
  
  if (length(j) != p) {
    stop(
      "Logical column index must have length equal to the number of columns.",
      call. = FALSE
    )
  }
  
  which(j)
}

#' @keywords internal
#' @noRd
spar_process_cols.character <- function(x, j, names = colnames(x), ...) {
  idx <- match(j, names)
  
  if (anyNA(idx)) {
    bad <- unique(j[is.na(idx)])
    stop(
      sprintf(
        "Unknown column name(s): %s",
        paste(bad, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  idx
}


#' @keywords internal
#' @noRd
spar_process_rows.default <- function(x, i, n = nrow(x), ...) {
  stop(
    sprintf(
      "Unsupported row index type: %s",
      paste(class(i), collapse = "/")
    ),
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
spar_process_rows.NULL <- function(x, i, n = nrow(x), ...) {
  seq_len(n)
}

#' @keywords internal
#' @noRd
spar_process_rows.numeric <- function(x, i, n = nrow(x), ...) {
  i <- as.integer(i)
  
  if (anyNA(i)) {
    stop("Row index contains NA values.", call. = FALSE)
  }
  
  if (any(i < 1L | i > n)) {
    stop("Row index out of bounds.", call. = FALSE)
  }
  
  i
}

#' @keywords internal
#' @noRd
spar_process_rows.logical <- function(x, i, n = nrow(x), ...) {
  if (length(i) != n) {
    stop(
      "Logical row index must have length equal to the number of rows.",
      call. = FALSE
    )
  }
  
  which(i)
}
