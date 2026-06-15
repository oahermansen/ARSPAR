#' Create a low-level SPAR transform context
#'
#' Constructs a `spar_transform_context` object describing the namespaces
#' available during transform compilation and execution.
#'
#' A transform context contains:
#' \itemize{
#'   \item `current`: the data immediately prior to the current transform step,
#'   \item `original`: the original-space data,
#'   \item `steps`: a named list of prior transform step outputs.
#' }
#'
#' This object is used to resolve SPAR transform references such as:
#' \itemize{
#'   \item `Hs` for current-state references,
#'   \item `.Hs` for original-space references,
#'   \item `center.Hs` for prior-step references.
#' }
#'
#' @param current A data frame-like object representing the current state.
#' @param original A data frame-like object representing the original-space
#'   state.
#' @param steps A named list of prior transform step outputs. Each element
#'   should be data frame-like.
#'
#' @return An object of class `"spar_transform_context"`.
#'
#' @keywords internal
new_spar_transform_context <- function(
    current,
    original,
    steps = list()
) {
  structure(
    list(
      current = current,
      original = original,
      steps = steps
    ),
    class = "spar_transform_context"
  )
}

#' Construct a SPAR transform chain
#'
#' Create a `spar_transform_chain` object representing an ordered sequence
#' of compiled transform steps. A chain provides a structured container
#' for composing multiple `spar_transform` objects and managing execution,
#' dependency analysis, and cache planning across the full transform pipeline.
#'
#' Each step in the chain must be a compiled `spar_transform` with a unique,
#' non-empty name. The order of steps determines the forward evaluation order
#' when the chain is executed.
#'
#' At construction time the chain stores only the provided steps and optional
#' metadata. Cache requirements and dependency analysis can be computed later
#' using helper functions such as `spar_chain_cache_requirements()` or
#' `spar_finalize_transform_chain()`.
#'
#' @param steps A list of objects inheriting from `"spar_transform"`.
#'   Each element must have a unique, non-empty `name`.
#' @param name Optional character scalar giving a descriptive name for the chain.
#' @param cache_requirements Optional precomputed cache requirement structure.
#'   Typically left `NULL` and populated later.
#' @param analysis Optional structure storing dependency analysis or compilation
#'   metadata for the chain.
#' @param inverse Optional inverse function for the full transform chain.
#'
#' @return
#' An object of class `"spar_transform_chain"` containing:
#'
#' * `name` – optional chain name
#' * `steps` – ordered list of `spar_transform` steps
#' * `cache_requirements` – optional cache analysis
#' * `analysis` – optional dependency metadata
#'
#' @details
#' The chain object is intended to support both forward execution and later
#' extension to reverse transforms. For this reason cache requirements and
#' dependency metadata are stored separately from the step list.
#'
#' @seealso
#' `spar_add_transform_step()`, `spar_run_transform_chain()`,
#' `spar_chain_cache_requirements()`
#'
#' @export
new_spar_transform_chain <- function(
    steps = list(),
    name = NULL,
    cache_requirements = NULL,
    analysis = NULL,
    inverse = NULL
) {
  stopifnot(is.list(steps))

  if (length(steps) > 0L) {
    ok <- vapply(steps, inherits, logical(1), "spar_transform")
    if (!all(ok)) {
      stop("All 'steps' must inherit from 'spar_transform'.", call. = FALSE)
    }

    step_names <- vapply(steps, function(s) s$name %||% "", character(1))
    if (any(step_names == "")) {
      stop("All transform steps must have non-empty names.", call. = FALSE)
    }
    if (anyDuplicated(step_names)) {
      stop("Transform step names must be unique within a chain.", call. = FALSE)
    }
  }

  if (!is.null(inverse) && !is.function(inverse)) {
    stop("`inverse` must be NULL or a function.", call. = FALSE)
  }

  structure(
    list(
      name = name,
      steps = steps,
      cache_requirements = cache_requirements,
      analysis = analysis,
      inverse = inverse
    ),
    class = "spar_transform_chain"
  )
}

#' Validate a SPAR transform chain
#'
#' Validate structural invariants of a `spar_transform_chain` object.
#' This function checks that the object has the correct class and that
#' the step list contains valid `spar_transform` objects with unique
#' non-empty names.
#'
#' The function is primarily intended for internal use but may also be
#' useful in debugging or development contexts.
#'
#' @param x Object expected to inherit from `"spar_transform_chain"`.
#'
#' @return
#' Invisibly returns `x` if validation succeeds.
#'
#' @details
#' The following invariants are enforced:
#'
#' * `x` inherits from `"spar_transform_chain"`
#' * `steps` is a list
#' * each element of `steps` inherits from `"spar_transform"`
#' * each step has a non-empty name
#' * step names are unique within the chain
#'
#' @export
validate_spar_transform_chain <- function(x) {
  if (!inherits(x, "spar_transform_chain")) {
    stop("Expected a 'spar_transform_chain' object.", call. = FALSE)
  }

  if (!is.list(x$steps)) {
    stop("'steps' must be a list.", call. = FALSE)
  }

  if (!is.null(x$inverse) && !is.function(x$inverse)) {
    stop("'inverse' must be NULL or a function.", call. = FALSE)
  }

  if (length(x$steps) > 0L) {
    ok <- vapply(x$steps, inherits, logical(1), "spar_transform")
    if (!all(ok)) {
      stop("All elements of 'steps' must inherit from 'spar_transform'.", call. = FALSE)
    }

    step_names <- vapply(x$steps, function(s) s$name %||% "", character(1))
    if (any(step_names == "")) {
      stop("Every step must have a non-empty name.", call. = FALSE)
    }
    if (anyDuplicated(step_names)) {
      stop("Step names must be unique.", call. = FALSE)
    }
  }

  invisible(x)
}

#' Print a SPAR transform chain
#'
#' Print a concise summary of a `spar_transform_chain`, including its name,
#' the number of steps, and the ordered list of step names.
#'
#' @param x A `spar_transform_chain` object.
#' @param ... Additional arguments passed to other methods (unused).
#'
#' @return
#' Invisibly returns `x`.
#'
#' @details
#' The printed summary includes the chain name (if present), the number
#' of steps, and the step ordering. If cache requirements have been
#' computed, the summary indicates that cache information is available.
#'
#' @export
print.spar_transform_chain <- function(x, ...) {
  validate_spar_transform_chain(x)

  cat("<spar_transform_chain>")
  if (!is.null(x$name)) cat(" ", x$name, sep = "")
  cat("\n")

  n <- length(x$steps)
  cat("  steps:", n, "\n")

  if (n > 0L) {
    nm <- vapply(x$steps, `[[`, character(1), "name")
    for (i in seq_along(nm)) {
      cat(sprintf("    %d. %s\n", i, nm[[i]]))
    }
  }

  if (!is.null(x$cache_requirements)) {
    cat("  cache requirements: available\n")
  }

  invisible(x)
}

#' Extract step names from a transform chain
#'
#' Return the ordered vector of step names contained in a
#' `spar_transform_chain`.
#'
#' @param chain A `spar_transform_chain` object.
#'
#' @return
#' A character vector of step names in forward evaluation order.
#'
#' @seealso
#' `spar_add_transform_step()`
#'
#' @export
spar_chain_step_names <- function(chain) {
  validate_spar_transform_chain(chain)
  vapply(chain$steps, `[[`, character(1), "name")
}

#' Add a step to a transform chain
#'
#' Append a compiled `spar_transform` step to an existing
#' `spar_transform_chain`. The step must have a unique name within the
#' chain.
#'
#' Adding a step invalidates any previously computed cache requirement
#' or dependency analysis metadata stored on the chain.
#'
#' @param chain A `spar_transform_chain` object.
#' @param step A compiled `spar_transform` to append.
#'
#' @return
#' The updated `spar_transform_chain`.
#'
#' @details
#' The new step is appended to the end of the chain and will therefore
#' be evaluated after all existing steps during forward execution.
#'
#' Step names must remain unique within the chain to preserve the
#' integrity of step-qualified references (`step.field`) used in
#' transform expressions.
#'
#' @seealso
#' `new_spar_transform_chain()`, `spar_chain_step_names()`
#'
#' @export
spar_add_transform_step <- function(chain, step) {
  validate_spar_transform_chain(chain)

  if (!inherits(step, "spar_transform")) {
    stop("'step' must inherit from 'spar_transform'.", call. = FALSE)
  }

  nm <- step$name %||% ""
  if (nm == "") {
    stop("Transform step must have a non-empty name.", call. = FALSE)
  }

  if (nm %in% spar_chain_step_names(chain)) {
    stop(sprintf("Step '%s' already exists in chain.", nm), call. = FALSE)
  }

  chain$steps[[length(chain$steps) + 1L]] <- step
  chain$cache_requirements <- NULL
  chain$analysis <- NULL
  validate_spar_transform_chain(chain)
}

#' Merge cache requirement structures
#'
#' Combine two cache requirement structures by taking the union of
#' required fields for each step.
#'
#' @param x First cache requirement structure.
#' @param y Second cache requirement structure.
#'
#' @return
#' A cache requirement structure containing the union of required
#' cached outputs for each step.
#'
#' @details
#' Cache requirements are represented as named lists where each element
#' corresponds to a step name and contains a character vector of fields
#' that must remain available in the step cache.
#'
#' When merging requirements, the resulting requirement for each step
#' is the union of the required fields from both inputs.
#'
#' This function is primarily used to combine forward and reverse
#' dependency requirements when constructing bidirectional transform
#' pipelines.
#'
#' @export
spar_merge_cache_requirements <- function(x, y) {
  nms <- union(names(x), names(y))
  out <- setNames(vector("list", length(nms)), nms)

  for (nm in nms) {
    lhs <- x[[nm]] %||% character(0)
    rhs <- y[[nm]] %||% character(0)
    out[[nm]] <- union(lhs, rhs)
  }

  out
}

#' Compute cache requirements for a transform chain
#'
#' Compute the cache requirements implied by the steps of a
#' `spar_transform_chain`.
#'
#' Cache requirements determine which intermediate step outputs must
#' remain available in the step cache during execution so that later
#' steps can evaluate their expressions.
#'
#' @param chain A `spar_transform_chain`.
#' @param direction Character scalar specifying which cache requirements
#'   to return. One of `"forward"` or `"combined"`.
#'
#' @return
#' A named list mapping step names to the fields that must remain cached.
#'
#' @details
#' Forward cache requirements are derived from runtime dependencies
#' recorded during transform step compilation. These dependencies reflect
#' which prior step outputs are referenced by later steps.
#'
#' The `"combined"` direction is reserved for future bidirectional
#' transform support and represents the union of forward and reverse
#' cache requirements.
#'
#' @seealso
#' `spar_step_cache_requirements()`, `spar_merge_cache_requirements()`
#'
#' @export
spar_chain_cache_requirements <- function(chain, direction = c("forward", "combined")) {
  direction <- match.arg(direction)
  validate_spar_transform_chain(chain)

  forward <- spar_step_cache_requirements(chain$steps)
  reverse <- NULL
  combined <- forward

  if (!is.null(reverse)) {
    combined <- spar_merge_cache_requirements(forward, reverse)
  }

  switch(
    direction,
    forward = forward,
    combined = combined
  )
}

#' Execute a transform chain
#'
#' Execute all steps in a `spar_transform_chain` in forward order,
#' starting from an original data matrix.
#'
#' @param chain A `spar_transform_chain`.
#' @param X_original Numeric matrix representing the original data
#'   representation.
#'
#' @return
#' A list containing the result of executing the transform chain,
#' including the final representation and any cached step outputs.
#'
#' @details
#' The chain is executed sequentially using the step execution
#' infrastructure provided by `spar_run_transform_steps()`. Step cache
#' trimming is applied according to runtime dependency analysis so that
#' only required intermediate outputs are retained.
#'
#' @seealso
#' `spar_run_transform_steps()`
#'
#' @export
spar_run_transform_chain <- function(chain, X_original) {
  validate_spar_transform_chain(chain)

  spar_run_transform_steps(
    X_original = X_original,
    steps = chain$steps
  )
}

#' Finalize a transform chain
#'
#' Compute and attach cache requirement metadata to a
#' `spar_transform_chain`.
#'
#' @param chain A `spar_transform_chain`.
#'
#' @return
#' The updated chain with cache requirement information populated.
#'
#' @details
#' This function computes forward cache requirements using
#' `spar_chain_cache_requirements()` and stores the result in the
#' `cache_requirements` field of the chain object.
#'
#' The stored structure contains:
#'
#' * `forward` – forward execution cache requirements
#' * `reverse` – placeholder for reverse transform requirements
#' * `combined` – union of forward and reverse requirements
#'
#' Reverse requirements are currently `NULL` but the structure is
#' included to maintain a stable API once reverse transforms are
#' implemented.
#'
#' @export
spar_finalize_transform_chain <- function(chain) {
  validate_spar_transform_chain(chain)

  chain$cache_requirements <- list(
    forward = spar_chain_cache_requirements(chain, "forward"),
    reverse = NULL,
    combined = spar_chain_cache_requirements(chain, "combined")
  )

  chain
}

#' Compile a sequence of transform steps into a chain
#'
#' Compile a list of transform step specifications into a fully
#' constructed `spar_transform_chain`.
#'
#' Each step specification is compiled in sequence using a progressively
#' enriched transform context so that later steps can reference outputs
#' from earlier steps via step-qualified references (`step.field`).
#'
#' @param step_specs A list of step specifications. Each specification
#'   must contain:
#'   * `name` – step name
#'   * `exprs` – named list of transform expressions
#' @param X_original Numeric matrix representing the original data.
#' @param name Optional chain name.
#'
#' @return
#' A finalized `spar_transform_chain` containing compiled transform
#' steps and cache requirement metadata.
#'
#' @details
#' Compilation proceeds sequentially:
#'
#' 1. A compile-time context is created from the current representation,
#'    original data, and outputs of previously compiled steps.
#' 2. The step is compiled using `spar_compile_transform_step()`.
#' 3. The step is evaluated to update the current representation.
#' 4. The resulting output is stored in the step cache for use by later
#'    steps.
#'
#' This process ensures that step-qualified references such as
#' `step.field` can be resolved during compilation.
#'
#' @seealso
#' `spar_compile_transform_step()`, `spar_run_transform_chain()`
#'
#' @export
spar_compile_transform_chain <- function(step_specs, X_original, name = NULL) {
  stopifnot(is.list(step_specs))

  X_original_df <- as.data.frame(X_original, check.names = FALSE)
  current_df <- X_original_df
  current_mat <- as.matrix(X_original)

  step_cache <- list()
  compiled_steps <- vector("list", length(step_specs))

  for (i in seq_along(step_specs)) {
    spec <- step_specs[[i]]

    step_name <- spec$name %||% NULL
    exprs <- spec$exprs %||% NULL

    if (is.null(step_name) || !nzchar(step_name)) {
      stop(sprintf("Step spec %d is missing a valid 'name'.", i), call. = FALSE)
    }
    if (is.null(exprs) || !is.list(exprs) || length(exprs) == 0L) {
      stop(sprintf("Step spec '%s' must provide non-empty 'exprs'.", step_name), call. = FALSE)
    }

    ctx <- new_spar_transform_context(
      current = current_df,
      original = X_original_df,
      steps = step_cache
    )

    step <- spar_compile_transform_step(
      exprs = exprs,
      name = step_name,
      context = ctx
    )

    compiled_steps[[i]] <- step

    current_mat <- spar_eval_transform_step(
      x = step,
      current = current_mat,
      original = X_original,
      steps = step_cache
    )

    current_df <- as.data.frame(current_mat, check.names = FALSE)
    step_cache[[step_name]] <- current_df
  }

  chain <- new_spar_transform_chain(
    steps = compiled_steps,
    name = name
  )

  spar_finalize_transform_chain(chain)
}

#' Build a mutate-style transform step specification
#'
#' Captures named mutate expressions into a step specification suitable for
#' [spar_compile_transform_chain()] or [spar_build_transform_chain()].
#'
#' @param ... Named expressions defining output columns.
#' @param .name Single non-empty character string naming the step.
#'
#' @return A list with class `"spar_transform_step_spec"`.
#'
#' @export
spar_step_mutate <- function(..., .name) {
  quos <- rlang::enquos(..., .named = TRUE)

  if (length(quos) == 0L) {
    stop("At least one named expression must be supplied.", call. = FALSE)
  }

  if (any(names(quos) == "")) {
    stop("All step expressions must be named.", call. = FALSE)
  }

  if (!is.character(.name) || length(.name) != 1L || is.na(.name) || !nzchar(.name)) {
    stop("`.name` must be a single non-empty character string.", call. = FALSE)
  }

  structure(
    list(
      name = .name,
      exprs = lapply(quos, rlang::quo_get_expr)
    ),
    class = "spar_transform_step_spec"
  )
}

spar_as_transform_step_spec <- function(x, idx = NULL) {
  if (!is.list(x)) {
    stop("Each step specification must be a list.", call. = FALSE)
  }

  step_name <- x$name %||% ""
  exprs <- x$exprs %||% NULL

  if (!is.character(step_name) || length(step_name) != 1L || !nzchar(step_name)) {
    tag <- if (is.null(idx)) "" else paste0(" at position ", idx)
    stop(sprintf("Invalid step name%s.", tag), call. = FALSE)
  }

  if (is.null(exprs) || !is.list(exprs) || length(exprs) == 0L) {
    stop(sprintf("Step '%s' must provide non-empty 'exprs'.", step_name), call. = FALSE)
  }

  if (is.null(names(exprs)) || any(names(exprs) == "")) {
    stop(sprintf("Step '%s' must provide named expressions.", step_name), call. = FALSE)
  }

  list(name = step_name, exprs = exprs)
}

#' Build and compile a transform chain ergonomically
#'
#' Accepts mutate-style step specifications and compiles them into a finalized
#' `spar_transform_chain`.
#'
#' @param X_original Original numeric data matrix.
#' @param ... Step specs, usually created by [spar_step_mutate()]. You may also
#'   pass a single list of step specs.
#' @param name Optional chain name.
#'
#' @return A finalized `spar_transform_chain`.
#'
#' @export
spar_build_transform_chain <- function(X_original, ..., name = NULL) {
  dots <- list(...)

  is_step_spec_like <- function(z) {
    is.list(z) && !is.null(z$name) && !is.null(z$exprs)
  }

  if (
    length(dots) == 1L &&
    is.list(dots[[1]]) &&
    !inherits(dots[[1]], "spar_transform_step_spec") &&
    !is_step_spec_like(dots[[1]])
  ) {
    candidates <- dots[[1]]
  } else {
    candidates <- dots
  }

  if (length(candidates) == 0L) {
    stop("At least one step specification is required.", call. = FALSE)
  }

  step_specs <- lapply(seq_along(candidates), function(i) {
    spar_as_transform_step_spec(candidates[[i]], idx = i)
  })

  step_names <- vapply(step_specs, function(z) z$name, character(1))
  if (anyDuplicated(step_names)) {
    stop("Step specification names must be unique.", call. = FALSE)
  }

  spar_compile_transform_chain(
    step_specs = step_specs,
    X_original = X_original,
    name = name
  )
}

#' Validate a SPAR transform context
#'
#' Checks that a `spar_transform_context` has the expected structure.
#'
#' @param x A `spar_transform_context` object.
#'
#' @return `TRUE`, invisibly, if validation succeeds.
#'
#' @keywords internal
validate_spar_transform_context <- function(x) {
  if (!inherits(x, "spar_transform_context")) {
    stop("Object is not a `spar_transform_context`.", call. = FALSE)
  }

  if (!is.list(x$current) && !is.data.frame(x$current)) {
    stop("`current` must be data frame-like.", call. = FALSE)
  }

  if (!is.list(x$original) && !is.data.frame(x$original)) {
    stop("`original` must be data frame-like.", call. = FALSE)
  }

  if (!is.list(x$steps)) {
    stop("`steps` must be a list.", call. = FALSE)
  }

  if (length(x$steps) > 0L) {
    nm <- names(x$steps)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("`steps` must be a named list.", call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Compile a single transform expression
#'
#' Compiles one mutate-style transform expression into a `spar_transform_expr`
#' object.
#'
#' Compilation performs:
#' \itemize{
#'   \item dependency analysis of the original expression,
#'   \item summary-product capture,
#'   \item expression rewriting,
#'   \item dependency analysis of the compiled expression,
#'   \item metadata assembly.
#' }
#'
#' The original expression dependencies describe source provenance, while the
#' compiled expression dependencies describe runtime requirements after summary
#' products have been captured into parameters.
#'
#' @param expr An R expression defining the transformed output column.
#' @param name Name of the output column produced by the expression.
#' @param context A `spar_transform_context` object.
#' @param data_names Optional character vector of current-state field names.
#'   Defaults to `names(context$current)`.
#' @param step_names Optional character vector of available prior transform step
#'   names. Defaults to `names(context$steps)`.
#'
#' @return An object of class `"spar_transform_expr"`.
#'
#' @keywords internal
spar_compile_transform_expr <- function(
    expr,
    name,
    context,
    data_names = names(context$current),
    step_names = names(context$steps)
) {
  if (!inherits(context, "spar_transform_context")) {
    stop("`context` must be a `spar_transform_context`.", call. = FALSE)
  }

  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }

  validate_spar_transform_context(context)

  source_dependencies <- spar_transform_dependencies(
    expr = expr,
    data_names = data_names,
    step_names = step_names
  )

  captured <- spar_capture_summary_products(
    expr = expr,
    context = context,
    data_names = data_names,
    step_names = step_names
  )

  runtime_dependencies <- spar_transform_dependencies(
    expr = captured$expr,
    data_names = data_names,
    step_names = step_names
  )

  structure(
    list(
      name = name,
      original_expr = expr,
      compiled_expr = captured$expr,
      params = captured$params,
      source_dependencies = source_dependencies,
      runtime_dependencies = runtime_dependencies,
      captured = captured$captured
    ),
    class = "spar_transform_expr"
  )
}

#' Compute runtime step-cache requirements for a transform chain
#'
#' Inspects an ordered list of `spar_transform` objects and determines which
#' prior-step fields are required at runtime by later steps.
#'
#' Only runtime dependencies are considered, meaning dependencies that remain
#' after summary-safe products have been captured into parameters.
#'
#' @param steps Ordered list of `spar_transform` objects.
#'
#' @return A named list. Each element is a character vector of field names that
#' must be retained for the corresponding step at runtime.
#'
#' @details
#' If a step is not referenced by any later compiled transform expression, its
#' runtime cache requirement is an empty character vector.
#'
#' @keywords internal
#' @export
spar_step_cache_requirements <- function(steps) {
  if (!is.list(steps)) {
    stop("`steps` must be a list.", call. = FALSE)
  }

  if (length(steps) == 0L) {
    return(list())
  }

  for (k in seq_along(steps)) {
    if (!inherits(steps[[k]], "spar_transform")) {
      stop("All elements of `steps` must be `spar_transform` objects.", call. = FALSE)
    }
  }

  step_names <- vapply(steps, function(s) s$name, character(1))

  reqs <- stats::setNames(vector("list", length(step_names)), step_names)
  for (nm in step_names) {
    reqs[[nm]] <- character(0)
  }

  for (k in seq_along(steps)) {
    exprs <- steps[[k]]$exprs

    for (j in seq_along(exprs)) {
      rt_deps <- exprs[[j]]$runtime_dependencies
      step_deps <- rt_deps$steps

      if (length(step_deps) == 0L) {
        next
      }

      for (dep_step in names(step_deps)) {
        reqs[[dep_step]] <- unique(c(reqs[[dep_step]], step_deps[[dep_step]]))
      }
    }
  }

  reqs
}

#' Trim a step output to required runtime cache fields
#'
#' Reduces a step output matrix or data frame to the subset of columns required
#' for future runtime references.
#'
#' @param x A matrix or data frame representing one step output.
#' @param fields Character vector of field names required for future runtime
#'   references.
#'
#' @return A data frame with only the required fields. If `fields` is empty, an
#'   empty named list is returned.
#'
#' @keywords internal
#' @export
spar_trim_step_cache <- function(x, fields) {
  if (length(fields) == 0L) {
    return(list())
  }

  x_df <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)

  missing_fields <- setdiff(fields, names(x_df))
  if (length(missing_fields) > 0L) {
    stop(
      sprintf(
        "Required cached field(s) not found in step output: %s",
        paste(missing_fields, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  x_df[, fields, drop = FALSE]
}

#' Evaluate a compiled transform expression
#'
#' Evaluates a `spar_transform_expr` against a transform context and returns the
#' resulting vector.
#'
#' Any remaining transform references in the compiled expression are resolved
#' against the supplied context. Fitted parameters stored in the compiled
#' expression are also made available during evaluation.
#'
#' @param x A `spar_transform_expr` object.
#' @param context A `spar_transform_context` object.
#' @param data_names Optional character vector of current-state field names.
#'   Defaults to `names(context$current)`.
#' @param step_names Optional character vector of available prior transform step
#'   names. Defaults to `names(context$steps)`.
#'
#' @return The evaluated result of the compiled expression, typically a vector
#' of length equal to the number of observations.
#'
#' @examples
#' ctx <- new_spar_transform_context(
#'   current = data.frame(Hs = 1:3, check.names = FALSE),
#'   original = data.frame(Hs = 4:6, check.names = FALSE)
#' )
#'
#' cx <- spar_compile_transform_expr(
#'   expr = quote(Hs - mean(Hs)),
#'   name = "Hs",
#'   context = ctx
#' )
#'
#' spar_eval_transform_expr(cx, ctx)
#'
#' @keywords internal
spar_eval_transform_expr <- function(
    x,
    context,
    data_names = names(context$current),
    step_names = names(context$steps)
) {
  if (!inherits(x, "spar_transform_expr")) {
    stop("`x` must be a `spar_transform_expr`.", call. = FALSE)
  }

  if (!inherits(context, "spar_transform_context")) {
    stop("`context` must be a `spar_transform_context`.", call. = FALSE)
  }

  bound <- spar_bind_transform_refs(
    expr = x$compiled_expr,
    context = context,
    data_names = data_names,
    step_names = step_names
  )

  eval_env <- rlang::env(parent = baseenv(), !!!x$params, !!!bound$bindings)
  rlang::eval_bare(bound$expr, env = eval_env)
}

#' Collect transform references from an expression
#'
#' Recursively traverses an expression and collects SPAR transform references.
#'
#' Function names are ignored. Only symbols in argument positions are collected.
#'
#' @param expr An R expression.
#' @param data_names Optional character vector of field names available in the
#'   current transform input context.
#' @param step_names Optional character vector of available prior transform step
#'   names.
#'
#' @return A list with components:
#' \describe{
#'   \item{refs}{List of parsed references.}
#'   \item{current}{Character vector of current-state fields referenced.}
#'   \item{original}{Character vector of original-space fields referenced.}
#'   \item{steps}{Named list mapping prior step names to referenced fields.}
#'   \item{data}{Logical scalar indicating use of the `.data` reference.}
#' }
#'
#' @examples
#' expr <- quote(Hs - mean(center.Tm.2) + .Hs)
#' spar_collect_transform_refs(
#'   expr,
#'   data_names = c("Hs", "Tm.2"),
#'   step_names = "center"
#' )
#'
#' @keywords internal
#' @export
spar_collect_transform_refs <- function(
    expr,
    data_names = NULL,
    step_names = NULL
) {
  refs <- list()

  recurse <- function(node, in_fun_position = FALSE) {
    if (rlang::is_symbol(node)) {
      nm <- rlang::as_string(node)
      if (!in_fun_position) {
        refs[[length(refs) + 1L]] <<- spar_parse_transform_ref(
          name = nm,
          data_names = data_names,
          step_names = step_names
        )
      }
      return(invisible(NULL))
    }

    if (rlang::is_call(node)) {
      args <- as.list(node)

      if (length(args) >= 1L) {
        recurse(args[[1L]], in_fun_position = TRUE)
      }

      if (length(args) >= 2L) {
        for (k in 2:length(args)) {
          recurse(args[[k]], in_fun_position = FALSE)
        }
      }

      return(invisible(NULL))
    }

    invisible(NULL)
  }

  recurse(expr)

  current_refs <- Filter(function(z) identical(z$type, "current"), refs)
  original_refs <- Filter(function(z) identical(z$type, "original"), refs)
  step_refs <- Filter(function(z) identical(z$type, "step"), refs)

  current <- unique(vapply(current_refs, function(z) z$field, character(1)))
  original <- unique(vapply(original_refs, function(z) z$field, character(1)))
  uses_data_ref <- length(Filter(function(z) identical(z$type, "data"), refs)) > 0L

  if (isTRUE(uses_data_ref) && !is.null(data_names)) {
    current <- unique(c(current, data_names))
  }

  steps <- list()
  if (length(step_refs) > 0L) {
    step_ids <- unique(vapply(step_refs, function(z) z$step, character(1)))
    steps <- stats::setNames(vector("list", length(step_ids)), step_ids)

    for (nm in step_ids) {
      steps[[nm]] <- unique(vapply(
        Filter(function(z) identical(z$step, nm), step_refs),
        function(z) z$field,
        character(1)
      ))
    }
  }

  list(refs = refs, current = current, original = original, steps = steps, data = uses_data_ref)
}

#' Return recognized summary-safe transform functions
#'
#' Returns the names of functions that are treated as dataset-level summary
#' functions during transform expression capture.
#'
#' Summary-safe calls may be evaluated once and stored as fitted parameters
#' instead of requiring rowwise products at runtime.
#'
#' @return Character vector of function names.
#'
#' @examples
#' spar_transform_summary_funs()
#'
#' @keywords internal
#' @export
spar_transform_summary_funs <- function() {
  c(
    "mean", "sd", "var", "min", "max", "median",
    "IQR", "mad", "sum", "prod",
    "quantile", "range",
    "colMeans", "colSums"
  )
}

#' Check whether an expression is a summary-safe transform call
#'
#' Determines whether an expression is a call to a recognized summary-safe
#' function.
#'
#' @param expr An R expression.
#'
#' @return Logical scalar.
#'
#' @examples
#' spar_is_summary_call(quote(mean(Hs)))
#' spar_is_summary_call(quote(Hs + 1))
#'
#' @keywords internal
#' @export
spar_is_summary_call <- function(expr) {
  if (!rlang::is_call(expr)) {
    return(FALSE)
  }

  fn <- rlang::call_name(expr)
  !is.null(fn) && fn %in% spar_transform_summary_funs()
}

#' Build dependency metadata for a transform expression
#'
#' Creates a compact dependency summary for a transform expression based on the
#' SPAR transform reference language.
#'
#' This helper is intended to be stored in transform metadata so that execution
#' code can determine whether a transform chain can be evaluated by simple
#' pipe-forward execution or whether intermediate prior-step outputs need to be
#' cached.
#'
#' @param expr An R expression.
#'
#' @return A named list with components:
#' \describe{
#'   \item{current}{Character vector of current-state fields referenced.}
#'   \item{original}{Character vector of original-space fields referenced.}
#'   \item{steps}{Named list of prior-step field dependencies.}
#'   \item{data}{Logical scalar indicating whether `.data` is referenced.}
#'   \item{needs_step_cache}{Logical scalar indicating whether prior-step
#'   references are present.}
#' }
#'
#' @examples
#' spar_transform_dependencies(
#'   quote(Hs + center.Tm.2 + .Hs),
#'   data_names = c("Hs", "Tm.2"),
#'   step_names = "center"
#' )
#'
#' @seealso [spar_collect_transform_refs()]
#'
#' @keywords internal
#' @export
spar_transform_dependencies <- function(
    expr,
    data_names = NULL,
    step_names = NULL
) {
  deps <- spar_collect_transform_refs(
    expr = expr,
    data_names = data_names,
    step_names = step_names
  )

  list(
    current = deps$current,
    original = deps$original,
    steps = deps$steps,
    data = isTRUE(deps$data),
    needs_step_cache = length(deps$steps) > 0L
  )
}

#' Parse a SPAR transform reference
#'
#' Parses a symbol name according to the SPAR transform reference language.
#'
#' Parsing is context-aware when `data_names` and `step_names` are supplied.
#' This is important because user field names may themselves contain periods.
#'
#' The parsing rules are:
#' \itemize{
#'   \item `.data` is treated as a reference to the full current data,
#'   \item names beginning with `.` are treated as original-data references,
#'   \item exact matches in `data_names` are treated as current-state fields,
#'   \item names of the form `step.field` are treated as step-qualified
#'     references only when `step` is a known step name,
#'   \item otherwise the full name is treated as a current-state field name.
#' }
#'
#' @param name A single character string.
#' @param data_names Optional character vector of field names available in the
#'   current transform input context.
#' @param step_names Optional character vector of available prior transform step
#'   names.
#'
#' @return A named list with components `type`, `step`, `field`, and `raw`.
#'
#' @examples
#' spar_parse_transform_ref("Hs")
#' spar_parse_transform_ref(".Hs")
#' spar_parse_transform_ref("center.Hs", step_names = "center")
#' spar_parse_transform_ref("Tm.2", data_names = "Tm.2", step_names = "Tm")
#'
#' @keywords internal
#' @export
spar_parse_transform_ref <- function(name, data_names = NULL, step_names = NULL) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }

  if (identical(name, ".data")) {
    return(list(type = "data", step = NULL, field = ".data", raw = name))
  }

  if (!is.null(data_names) && name %in% data_names) {
    return(list(type = "current", step = NULL, field = name, raw = name))
  }

  if (startsWith(name, ".")) {
    field <- substring(name, 2L)
    if (!nzchar(field)) {
      stop("Invalid original-data reference.", call. = FALSE)
    }
    return(list(type = "original", step = NULL, field = field, raw = name))
  }

  if (grepl(".", name, fixed = TRUE)) {
    pos <- regexpr(".", name, fixed = TRUE)[1L]
    prefix <- substring(name, 1L, pos - 1L)
    suffix <- substring(name, pos + 1L)

    if (!is.null(step_names) && prefix %in% step_names && nzchar(suffix)) {
      return(list(type = "step", step = prefix, field = suffix, raw = name))
    }
  }

  list(type = "current", step = NULL, field = name, raw = name)
}

#' Warn about transform step names colliding with dotted field-name prefixes
#'
#' Emits a warning when a proposed transform step name matches the prefix of an
#' existing dotted field name.
#'
#' @param step_name Proposed transform step name.
#' @param field_names Character vector of existing field names.
#'
#' @return Invisibly `TRUE`.
#'
#' @keywords internal
spar_warn_step_name_prefix_collision <- function(step_name, field_names) {
  dotted <- field_names[grepl(".", field_names, fixed = TRUE)]
  if (length(dotted) == 0L) {
    return(invisible(TRUE))
  }

  prefixes <- sub("\\..*$", "", dotted)
  hits <- dotted[prefixes == step_name]

  if (length(hits) > 0L) {
    warning(
      sprintf(
        "Transform step name '%s' matches the prefix of dotted field name(s): %s. ",
        step_name,
        paste(unique(hits), collapse = ", ")
      ),
      "This is supported, but may be confusing because SPAR uses 'step.field' reference syntax.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Check for dotted field names
#'
#' Warns when names contain `.` and returns metadata describing them.
#'
#' @param names Character vector of field names.
#' @param what Description of the naming context, e.g. `"original"` or
#'   `"transformed"`.
#'
#' @return A named list with dotted-name metadata.
#'
#' @keywords internal
spar_check_dotted_names <- function(names, what = "field") {
  dotted <- names[grepl(".", names, fixed = TRUE)]

  if (length(dotted) > 0L) {
    warning(
      sprintf(
        "%s names containing '.' are supported but discouraged in SPAR: %s. ",
        tools::toTitleCase(what),
        paste(dotted, collapse = ", ")
      ),
      "The '.' character is also used in transform references such as 'step.field'.",
      call. = FALSE
    )
  }

  list(
    has_dotted_names = length(dotted) > 0L,
    dotted_names = unique(dotted)
  )
}

#' Resolve a parsed transform reference against a transform context
#'
#' Resolves a parsed SPAR transform reference against a transform context.
#'
#' The context is expected to provide:
#' \itemize{
#'   \item `current`: current-state data,
#'   \item `original`: original-space data,
#'   \item `steps`: named list of prior step outputs.
#' }
#'
#' @param ref A parsed transform reference as returned by
#'   [spar_parse_transform_ref()].
#' @param context A named list containing `current`, `original`, and `steps`.
#'
#' @return The resolved field vector or object.
#'
#' @examples
#' ctx <- list(
#'   current = data.frame(Hs = 1:3, check.names = FALSE),
#'   original = data.frame(Hs = 4:6, check.names = FALSE),
#'   steps = list(center = data.frame("Tm.2" = 7:9, check.names = FALSE))
#' )
#'
#' spar_resolve_transform_ref(
#'   spar_parse_transform_ref("Hs", data_names = "Hs"),
#'   ctx
#' )
#'
#' @keywords internal
#' @export
spar_resolve_transform_ref <- function(ref, context) {
  if (!is.list(ref) || is.null(ref$type) || is.null(ref$field)) {
    stop("`ref` must be a parsed transform reference.", call. = FALSE)
  }

  if (!is.list(context)) {
    stop("`context` must be a list.", call. = FALSE)
  }

  get_field <- function(obj, field, what) {
    if (is.null(obj)) {
      stop(sprintf("Transform context '%s' is not available.", what), call. = FALSE)
    }

    if (!field %in% names(obj)) {
      stop(
        sprintf("Field '%s' not found in transform context '%s'.", field, what),
        call. = FALSE
      )
    }

    obj[[field]]
  }

  switch(
    ref$type,
    current = get_field(context$current, ref$field, "current"),
    original = get_field(context$original, ref$field, "original"),
    data = as.data.frame(context$current, check.names = FALSE, stringsAsFactors = FALSE),
    step = {
      if (is.null(ref$step) || !nzchar(ref$step)) {
        stop("Invalid step-qualified transform reference.", call. = FALSE)
      }

      if (is.null(context$steps) || !ref$step %in% names(context$steps)) {
        stop(
          sprintf("Referenced transform step '%s' is not available in context.", ref$step),
          call. = FALSE
        )
      }

      get_field(context$steps[[ref$step]], ref$field, ref$step)
    },
    stop("Unknown transform reference type.", call. = FALSE)
  )
}

#' Bind transform references in an expression
#'
#' Rewrites SPAR transform references in an expression into temporary symbols and
#' returns the corresponding bound values.
#'
#' Internal compiler parameter symbols beginning with `.param_` are left
#' unchanged.
#'
#' @param expr An R expression.
#' @param context A transform context list containing `current`, `original`, and
#'   `steps`.
#' @param data_names Optional character vector of field names available in the
#'   current transform input context.
#' @param step_names Optional character vector of available prior transform step
#'   names.
#'
#' @return A named list with components:
#' \describe{
#'   \item{expr}{The rewritten expression.}
#'   \item{bindings}{A named list of temporary bindings needed to evaluate the
#'   rewritten expression.}
#' }
#'
#' @keywords internal
#' @export
spar_bind_transform_refs <- function(
    expr,
    context,
    data_names = NULL,
    step_names = NULL
) {
  bindings <- list()
  counter <- 0L

  rewrite <- function(node, in_fun_position = FALSE) {
    if (rlang::is_symbol(node)) {
      if (in_fun_position) {
        return(node)
      }

      nm <- rlang::as_string(node)

      # Internal compiler parameters are already resolvable in the evaluation
      # environment and should not be treated as transform references.
      if (startsWith(nm, ".param_")) {
        return(node)
      }

      ref <- spar_parse_transform_ref(
        name = nm,
        data_names = data_names,
        step_names = step_names
      )

      counter <<- counter + 1L
      bind_name <- paste0(".ref_", counter)

      bindings[[bind_name]] <<- spar_resolve_transform_ref(ref, context)
      return(rlang::sym(bind_name))
    }

    if (rlang::is_call(node)) {
      args <- as.list(node)
      out_args <- vector("list", length(args))

      if (length(args) >= 1L) {
        out_args[[1L]] <- rewrite(args[[1L]], in_fun_position = TRUE)
      }

      if (length(args) >= 2L) {
        for (k in 2:length(args)) {
          out_args[[k]] <- rewrite(args[[k]], in_fun_position = FALSE)
        }
      }

      return(as.call(out_args))
    }

    node
  }

  list(
    expr = rewrite(expr),
    bindings = bindings
  )
}

#' Capture summary products from a transform expression
#'
#' Recursively traverses a transform expression and replaces recognized
#' summary-safe calls with parameter symbols, storing the evaluated summary
#' products in `params`.
#'
#' In this implementation, summary-safe calls are first rewritten so that any
#' SPAR transform references are resolved against the supplied transform
#' context. This allows calls such as `mean(.Hs)` and `mean(center.Tm.2)` to be
#' evaluated correctly.
#'
#' Captured summary products store dependency metadata using the name
#' `dependencies`. Raw parsed references are stored separately as
#' `parsed_refs`.
#'
#' @param expr An R expression.
#' @param context A transform context list containing `current`, `original`, and
#'   `steps`.
#' @param data_names Optional character vector of current input field names.
#' @param step_names Optional character vector of available prior transform step
#'   names.
#' @param params Named list of already captured parameters. If `NULL`, an empty
#'   list is created.
#'
#' @return A named list with components:
#' \describe{
#'   \item{expr}{The rewritten expression.}
#'   \item{params}{Updated parameter list.}
#'   \item{captured}{A list describing captured summary products.}
#' }
#'
#' Each element of `captured` is a list with components:
#' \describe{
#'   \item{call}{The original captured summary call.}
#'   \item{rewritten_call}{The rewritten call after binding transform references.}
#'   \item{param}{The parameter name used to store the captured summary value.}
#'   \item{dependencies}{Dependency summary for the captured product.}
#'   \item{parsed_refs}{Raw parsed transform references used by the captured
#'   product.}
#' }
#'
#' @examples
#' ctx <- list(
#'   current = data.frame(Hs = c(1, 2, 3), "Tm.2" = c(4, 5, 6), check.names = FALSE),
#'   original = data.frame(Hs = c(10, 20, 30), "Tm.2" = c(40, 50, 60), check.names = FALSE),
#'   steps = list(
#'     center = data.frame(Hs = c(-1, 0, 1), "Tm.2" = c(-1, 0, 1), check.names = FALSE)
#'   )
#' )
#'
#' spar_capture_summary_products(
#'   quote(Hs - mean(center.Tm.2)),
#'   context = ctx,
#'   data_names = c("Hs", "Tm.2"),
#'   step_names = "center"
#' )
#'
#' @keywords internal
#' @export
spar_capture_summary_products <- function(
    expr,
    context,
    data_names = NULL,
    step_names = NULL,
    params = NULL
) {
  if (is.null(params)) {
    params <- list()
  }

  if (!is.list(params)) {
    stop("`params` must be NULL or a list.", call. = FALSE)
  }

  if (is.null(names(params))) {
    names(params) <- rep.int("", length(params))
  }

  captured <- list()

  next_param_name <- function(params) {
    existing <- names(params)
    existing <- existing[nzchar(existing)]

    base <- ".param_fit_"

    if (length(existing) == 0L) {
      return(paste0(base, "1"))
    }

    hit <- grepl("^\\.param_fit_[0-9]+$", existing)
    if (!any(hit)) {
      return(paste0(base, "1"))
    }

    idx <- as.integer(sub("^\\.param_fit_", "", existing[hit]))
    paste0(base, max(idx, na.rm = TRUE) + 1L)
  }

  capture_call <- function(node, params) {
    if (spar_is_summary_call(node)) {
      bound <- spar_bind_transform_refs(
        expr = node,
        context = context,
        data_names = data_names,
        step_names = step_names
      )

      eval_env <- rlang::env(parent = baseenv(), !!!bound$bindings)
      value <- rlang::eval_bare(bound$expr, env = eval_env)

      param_name <- next_param_name(params)
      params[[param_name]] <- value

      deps <- spar_collect_transform_refs(
        expr = node,
        data_names = data_names,
        step_names = step_names
      )

      dependencies <- list(
        current = deps$current,
        original = deps$original,
        steps = deps$steps,
        needs_step_cache = length(deps$steps) > 0L
      )

      captured[[length(captured) + 1L]] <<- list(
        call = node,
        rewritten_call = bound$expr,
        param = param_name,
        dependencies = dependencies,
        parsed_refs = deps$refs
      )

      return(list(
        expr = rlang::sym(param_name),
        params = params
      ))
    }

    if (!rlang::is_call(node)) {
      return(list(expr = node, params = params))
    }

    args <- as.list(node)
    out_args <- vector("list", length(args))
    out_args[[1L]] <- args[[1L]]

    if (length(args) >= 2L) {
      for (k in 2:length(args)) {
        res <- capture_call(args[[k]], params = params)
        out_args[[k]] <- res$expr
        params <- res$params
      }
    }

    list(expr = as.call(out_args), params = params)
  }

  out <- capture_call(expr, params = params)

  list(
    expr = out$expr,
    params = out$params,
    captured = captured
  )
}
