#' Build a mutate-style SPAR transform step from expressions
#'
#' Captures named mutate-style expressions and compiles them into a
#' `spar_transform` object.
#'
#' @param x A `spar_representation` object.
#' @param ... Named expressions defining output columns for the transform step.
#' @param .name Single non-empty character string naming the transform step.
#'
#' @return An object of class `"spar_transform"`.
#'
#' @details
#' The compile-time context uses transformed data when available in
#' `x$transform$transformed_data`; otherwise it falls back to original data.
#' Prior-step references are only available if matching cached step outputs are
#' explicitly provided elsewhere; this function builds a single step object.
#'
#' @examples
#' # step <- spar_transform_mutate(
#' #   x,
#' #   Hs = Hs - mean(Hs),
#' #   Tm = Tm - mean(Tm),
#' #   .name = "center"
#' # )
#'
#' @export
spar_transform_mutate <- function(
    x,
    ...,
    .name
) {
  if (!inherits(x, "spar_representation")) {
    stop("`x` must be a `spar_representation` object.", call. = FALSE)
  }
  
  quos <- rlang::enquos(..., .named = TRUE)
  
  if (length(quos) == 0L) {
    stop("At least one named expression must be supplied.", call. = FALSE)
  }
  
  if (any(names(quos) == "")) {
    stop("All transform expressions must be named.", call. = FALSE)
  }
  
  if (!is.character(.name) || length(.name) != 1L || is.na(.name) || !nzchar(.name)) {
    stop("`.name` must be a single non-empty character string.", call. = FALSE)
  }
  
  Xt <- x$transform$transformed_data
  current_df <- if (is.null(Xt)) {
    as.data.frame(x$data$X_original, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    as.data.frame(Xt, check.names = FALSE, stringsAsFactors = FALSE)
  }
  
  original_df <- as.data.frame(
    x$data$X_original,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  ctx <- new_spar_transform_context(
    current = current_df,
    original = original_df,
    steps = list()
  )
  
  exprs <- lapply(quos, rlang::quo_get_expr)
  
  spar_compile_transform_step(
    exprs = exprs,
    name = .name,
    context = ctx,
    data_names = names(ctx$current),
    step_names = names(ctx$steps)
  )
}

#' Create a low-level SPAR transform step
#'
#' Constructs a `spar_transform` object representing a single transform step in
#' a SPAR transformation chain.
#'
#' A transform step stores:
#' \itemize{
#'   \item a step name,
#'   \item output column names,
#'   \item compiled transform expressions,
#'   \item a forward transform function,
#'   \item an optional inverse function,
#'   \item and step metadata.
#' }
#'
#' @param name Single non-empty character string naming the transform step.
#' @param output_names Character vector of output column names.
#' @param exprs Named list of `spar_transform_expr` objects.
#' @param forward Function implementing the forward transform.
#' @param inverse Optional inverse function. If `NULL`, the step is treated as
#'   non-invertible.
#' @param metadata Named list of additional metadata.
#'
#' @return An object of class `"spar_transform"`.
#'
#' @keywords internal
new_spar_transform <- function(
    name,
    output_names,
    exprs,
    forward,
    inverse = NULL,
    metadata = list()
) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }
  
  if (!is.character(output_names) || length(output_names) == 0L || anyNA(output_names)) {
    stop("`output_names` must be a non-empty character vector.", call. = FALSE)
  }
  
  if (!is.list(exprs) || length(exprs) == 0L) {
    stop("`exprs` must be a non-empty list.", call. = FALSE)
  }
  
  if (!is.function(forward)) {
    stop("`forward` must be a function.", call. = FALSE)
  }
  
  if (!is.null(inverse) && !is.function(inverse)) {
    stop("`inverse` must be NULL or a function.", call. = FALSE)
  }
  
  structure(
    list(
      name = name,
      output_names = output_names,
      exprs = exprs,
      forward = forward,
      inverse = inverse,
      metadata = metadata,
      invertible = !is.null(inverse)
    ),
    class = "spar_transform"
  )
}

#' Validate a SPAR transform step
#'
#' Checks that a `spar_transform` object has the expected structure.
#'
#' @param x A `spar_transform` object.
#'
#' @return `TRUE`, invisibly, if validation succeeds.
#'
#' @keywords internal
validate_spar_transform <- function(x) {
  if (!inherits(x, "spar_transform")) {
    stop("Object is not a `spar_transform`.", call. = FALSE)
  }
  
  if (!is.character(x$name) || length(x$name) != 1L) {
    stop("Invalid transform name.", call. = FALSE)
  }
  
  if (!is.character(x$output_names) || length(x$output_names) == 0L) {
    stop("Invalid transform output names.", call. = FALSE)
  }
  
  if (!is.list(x$exprs) || length(x$exprs) == 0L) {
    stop("Invalid transform expressions.", call. = FALSE)
  }
  
  if (!is.function(x$forward)) {
    stop("Invalid transform forward function.", call. = FALSE)
  }
  
  if (!is.null(x$inverse) && !is.function(x$inverse)) {
    stop("Invalid transform inverse function.", call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Compile a mutate-style SPAR transform step
#'
#' Compiles a named list of mutate-style expressions into a single
#' `spar_transform` object.
#'
#' @param exprs A named list of expressions.
#' @param name Single non-empty character string naming the transform step.
#' @param context A `spar_transform_context` object used during compilation.
#' @param data_names Optional character vector of current-state field names.
#' @param step_names Optional character vector of available prior transform step
#'   names.
#'
#' @return An object of class `"spar_transform"`.
#'
#' @keywords internal
spar_compile_transform_step <- function(
    exprs,
    name,
    context,
    data_names = names(context$current),
    step_names = names(context$steps)
) {
  if (!inherits(context, "spar_transform_context")) {
    stop("`context` must be a `spar_transform_context`.", call. = FALSE)
  }
  
  validate_spar_transform_context(context)
  
  if (is.null(exprs) || !is.list(exprs) || length(exprs) == 0L) {
    stop("`exprs` must be a non-empty named list of expressions.", call. = FALSE)
  }
  
  if (is.null(names(exprs)) || any(names(exprs) == "")) {
    stop("`exprs` must be a named list of expressions.", call. = FALSE)
  }
  
  compiled <- lapply(names(exprs), function(nm) {
    spar_compile_transform_expr(
      expr = exprs[[nm]],
      name = nm,
      context = context,
      data_names = data_names,
      step_names = step_names
    )
  })
  names(compiled) <- names(exprs)

  output_names <- names(exprs)
  uses_data_output <- identical(output_names, ".data")

  if (isTRUE(uses_data_output)) {
    if (length(compiled) != 1L) {
      stop("The special `.data` output must be defined as a single expression.", call. = FALSE)
    }

    output_names <- data_names
  }

  forward <- function(current, original = NULL, steps = list()) {
    if (is.null(original)) {
      original <- current
    }
    
    current_df <- as.data.frame(current, check.names = FALSE, stringsAsFactors = FALSE)
    original_df <- as.data.frame(original, check.names = FALSE, stringsAsFactors = FALSE)
    
    step_ctx <- new_spar_transform_context(
      current = current_df,
      original = original_df,
      steps = steps
    )
    
    n_obs <- nrow(current_df)
    if (isTRUE(uses_data_output)) {
      value <- spar_eval_transform_expr(
        x = compiled[[1L]],
        context = step_ctx,
        data_names = names(step_ctx$current),
        step_names = names(step_ctx$steps)
      )

      out_df <- as.data.frame(value, check.names = FALSE, stringsAsFactors = FALSE)

      if (nrow(out_df) != n_obs) {
        stop(
          sprintf(
            "Transform step '%s' `.data` output produced %d rows; expected %d.",
            name,
            nrow(out_df),
            n_obs
          ),
          call. = FALSE
        )
      }

      if (ncol(out_df) != length(output_names)) {
        stop(
          sprintf(
            "Transform step '%s' `.data` output produced %d columns; expected %d.",
            name,
            ncol(out_df),
            length(output_names)
          ),
          call. = FALSE
        )
      }

      colnames(out_df) <- output_names
      out_mat <- as.matrix(out_df)
    } else {
      out <- vector("list", length(compiled))
      names(out) <- output_names

      for (nm in output_names) {
        value <- spar_eval_transform_expr(
          x = compiled[[nm]],
          context = step_ctx,
          data_names = names(step_ctx$current),
          step_names = names(step_ctx$steps)
        )

        if (is.matrix(value) || is.data.frame(value)) {
          if (ncol(value) != 1L) {
            stop(
              sprintf(
                "Transform step '%s' output '%s' must evaluate to a scalar, vector, or one-column object.",
                name,
                nm
              ),
              call. = FALSE
            )
          }
          value <- value[, 1L, drop = TRUE]
        }

        if (length(value) == 1L && n_obs > 1L) {
          value <- rep(value, n_obs)
        }

        if (length(value) != n_obs) {
          stop(
            sprintf(
              "Transform step '%s' output '%s' has length %d; expected %d.",
              name,
              nm,
              length(value),
              n_obs
            ),
            call. = FALSE
          )
        }

        out[[nm]] <- value
      }

      out_df <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
      out_mat <- as.matrix(out_df)
    }
    
    if (!is.numeric(out_mat)) {
      stop("Compiled transform step produced non-numeric output.", call. = FALSE)
    }

    if (nrow(out_mat) != n_obs) {
      stop(
        sprintf(
          "Transform step '%s' produced %d rows; expected %d.",
          name,
          nrow(out_mat),
          n_obs
        ),
        call. = FALSE
      )
    }
    
    colnames(out_mat) <- output_names

    if (!identical(colnames(out_mat), output_names)) {
      stop(
        sprintf("Transform step '%s' produced unexpected output columns.", name),
        call. = FALSE
      )
    }

    out_mat
  }
  
  runtime_dependencies <- lapply(compiled, function(z) z$runtime_dependencies)
  
  metadata <- list(
    kind = "mutate",
    runtime_dependencies = runtime_dependencies,
    needs_step_cache = any(vapply(
      runtime_dependencies,
      function(z) isTRUE(z$needs_step_cache),
      logical(1)
    ))
  )
  
  tr <- new_spar_transform(
    name = name,
    output_names = output_names,
    exprs = compiled,
    forward = forward,
    inverse = NULL,
    metadata = metadata
  )
  
  validate_spar_transform(tr)
  tr
}

#' Print a SPAR transform step
#'
#' Prints a concise summary of a compiled `spar_transform` step.
#'
#' @param x A `spar_transform` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.spar_transform <- function(x, ...) {
  validate_spar_transform(x)

  cat("<spar_transform>", x$name, "\n")
  cat("  outputs:", paste(x$output_names, collapse = ", "), "\n")

  rt <- x$metadata$runtime_dependencies %||% list()
  needs_cache <- isTRUE(x$metadata$needs_step_cache)
  cat("  needs step cache:", if (needs_cache) "yes" else "no", "\n")

  if (length(rt) > 0L) {
    step_refs <- unique(unlist(lapply(rt, function(z) names(z$steps %||% list()))))
    if (length(step_refs) > 0L) {
      cat("  references steps:", paste(step_refs, collapse = ", "), "\n")
    }
  }

  invisible(x)
}

#' Summarize diagnostics for a transform step
#'
#' Builds a compact diagnostics object describing dependencies and captured
#' summary products for each output expression in a compiled `spar_transform`.
#'
#' @param x A `spar_transform` object.
#'
#' @return A list with class `"spar_transform_diagnostics"`.
#'
#' @export
spar_transform_diagnostics <- function(x) {
  validate_spar_transform(x)

  per_output <- lapply(x$exprs, function(expr_obj) {
    src <- expr_obj$source_dependencies %||% list()
    rt <- expr_obj$runtime_dependencies %||% list()

    list(
      output = expr_obj$name,
      source_current = src$current %||% character(0),
      source_original = src$original %||% character(0),
      source_steps = src$steps %||% list(),
      runtime_current = rt$current %||% character(0),
      runtime_original = rt$original %||% character(0),
      runtime_steps = rt$steps %||% list(),
      captured_count = length(expr_obj$captured %||% list())
    )
  })

  structure(
    list(
      name = x$name,
      outputs = x$output_names,
      needs_step_cache = isTRUE(x$metadata$needs_step_cache),
      per_output = per_output
    ),
    class = "spar_transform_diagnostics"
  )
}

#' Print transform-step diagnostics
#'
#' @param x A `spar_transform_diagnostics` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.spar_transform_diagnostics <- function(x, ...) {
  cat("<spar_transform_diagnostics>", x$name, "\n")
  cat("  outputs:", paste(x$outputs, collapse = ", "), "\n")
  cat("  needs step cache:", if (isTRUE(x$needs_step_cache)) "yes" else "no", "\n")

  for (item in x$per_output) {
    step_refs <- names(item$runtime_steps)
    cat(
      sprintf(
        "  - %s: captured=%d runtime_steps=%s\n",
        item$output,
        item$captured_count,
        if (length(step_refs) == 0L) "none" else paste(step_refs, collapse = ",")
      )
    )
  }

  invisible(x)
}

#' Evaluate a compiled SPAR transform step
#'
#' Evaluates a single compiled `spar_transform` step against a current state,
#' original state, and optional prior-step cache.
#'
#' @param x A `spar_transform` object.
#' @param current Current-state data, coercible to a data frame or matrix.
#' @param original Original-space data, coercible to a data frame or matrix.
#'   If `NULL`, `current` is used.
#' @param steps Named list of prior step outputs available for step-qualified
#'   references.
#'
#' @return A numeric matrix containing the output of the transform step.
#'
#' @keywords internal
#' @export
spar_eval_transform_step <- function(
    x,
    current,
    original = NULL,
    steps = list()
) {
  if (!inherits(x, "spar_transform")) {
    stop("`x` must be a `spar_transform`.", call. = FALSE)
  }
  
  if (is.null(original)) {
    original <- current
  }
  
  x$forward(
    current = current,
    original = original,
    steps = steps
  )
}

#' Run an ordered chain of SPAR transform steps
#'
#' Applies an ordered list of compiled `spar_transform` steps to original data.
#'
#' Each step is evaluated in sequence, with its output becoming the current
#' state for the next step. Prior step outputs are cached only to the extent
#' required by later runtime step-qualified references.
#'
#' @param X_original Original-space data as a numeric matrix or matrix-like
#'   object.
#' @param steps Ordered list of `spar_transform` objects.
#'
#' @return A named list with components:
#' \describe{
#'   \item{current}{The final transformed matrix.}
#'   \item{step_cache}{A named list of intermediate step outputs retained for
#'   runtime step-qualified references.}
#'   \item{cache_requirements}{A named list describing which fields were deemed
#'   necessary to cache for each step.}
#' }
#'
#' @details
#' This implementation uses runtime dependency metadata from compiled transform
#' expressions. Step outputs needed only for compile-time summary products do
#' not need to remain cached at runtime.
#'
#' @keywords internal
#' @export
spar_run_transform_steps <- function(
    X_original,
    steps
) {
  X_current <- X_original
  
  if (!is.matrix(X_current)) {
    X_current <- as.matrix(X_current)
  }
  
  if (!is.numeric(X_current)) {
    stop("`X_original` must be numeric or coercible to a numeric matrix.", call. = FALSE)
  }
  
  if (length(steps) == 0L) {
    return(list(
      current = X_current,
      step_cache = list(),
      cache_requirements = list()
    ))
  }
  
  if (!is.list(steps)) {
    stop("`steps` must be a list.", call. = FALSE)
  }
  
  for (k in seq_along(steps)) {
    if (!inherits(steps[[k]], "spar_transform")) {
      stop("All elements of `steps` must be `spar_transform` objects.", call. = FALSE)
    }
  }
  
  cache_requirements <- spar_step_cache_requirements(steps)
  
  X_orig <- X_current
  step_cache <- list()
  
  for (k in seq_along(steps)) {
    step <- steps[[k]]
    
    X_current <- spar_eval_transform_step(
      x = step,
      current = X_current,
      original = X_orig,
      steps = step_cache
    )
    
    req_fields <- cache_requirements[[step$name]]
    step_cache[[step$name]] <- spar_trim_step_cache(X_current, req_fields)
  }
  
  list(
    current = X_current,
    step_cache = step_cache,
    cache_requirements = cache_requirements
  )
}
