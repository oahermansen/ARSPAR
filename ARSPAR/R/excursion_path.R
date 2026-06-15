
# ============================================================
# 1. Radial / angular transform interface
#
# Replace these with your actual gauge and pseudo-angle map.
# Example below: 2D L1 pseudo-angle
#   g(x)     = x1 + x2
#   phi(x)   = x1 / (x1 + x2)

L1_radial <- function(X) {
  # X: matrix n x d
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  abs(X[, 1]) + abs(X[, 2])
}

L1_angle <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  s <- abs(X[, 1]) + abs(X[, 2])
  a <- sign(X[, 2])
  a[a==0] <- 1
  out <- (a*(1-X[, 1]/s))
  out[s == 0] <- NA_real_
  out
}

transform_points <- function(X, radial_fun = L1_radial, angle_fun = L1_angle, u) {
  R <- radial_fun(X)
  phi <- angle_fun(X)
  Y <- pmax(R - u, 0)
  list(R = R, phi = phi, Y = Y)
}

new_exceedance_index <- function(idx, t = idx, R = NULL,
                                 phi = NULL, u = NULL,
                                 excess = NULL,
                                 metadata = list()){
  structure(
    list(
      idx = idx,
      t = t,
      R = R,
      phi = phi,
      u = u,
      excess = excess,
      metadata = metadata
    ),
    class = "exceedance_index"
  )
}

extract_exceedance_index <- function(...){
  UseMethod("extract_exceedance_index")
}

extract_exceedance_index.default <- function(...){
  stop("Not implemented.")
}
extract_exceedance_index.data.frame <- function(data, mapping = ep_map(R=R, phi=phi, u=u, t=t, excess=excess)){
  datamap <- resolve_mapping(mapping)
  stopifnot(all(c("R", "phi", "u", "t", "excess") %in% datamap))
  mapped_df <- resolve_mapped_df(data, mapping)
  stopifnot(!is.null(mapped_df$R), !is.null(mapped_df$excess), !is.null(mapped_df$u),
            !is.null(mapped_df$phi), length(mapped_df$R) == length(mapped_df$excess),
            length(mapped_df$u)==1 || length(mapped_df$u) == length(mapped_df$R),
            length(mapped_df$t)==length(mapped_df$R), length(mapped_df$phi)== length(mapped_df$R))
  keep <- mapped_df$excess > 0
  new_exceedance_index(idx=mapped_df$t[keep], t=mapped_df$t[keep], R=mapped_df$R[keep],
                       phi=mapped_df$phi[keep], u=mapped_df$u[keep], excess=mapped_df$excess[keep])
}

extract_exceedance_spans <- function(exc_index, angle_domain = NULL) {
  dom <- NULL
  if (!is.null(angle_domain)) {
    dom <- as_spar_angle_domain(angle_domain)
  }

  if (length(exc_index$t) == 0L) {
    return(new_exceedance_spans(data.frame(
      span_id = integer(0),
      First = numeric(0),
      Last = numeric(0),
      n_points = integer(0),
      lphi = numeric(0),
      uphi = numeric(0)
    )))
  }

  id <- 1
  span_id <- c()
  First <- c()
  Last <- c()
  n_points <- c()
  lphi <- c()
  uphi <- c()
  first <- 0
  prev <- 0
  cnt <- 0
  span_start_pos <- 1L
  pos <- 0L
  phi_all <- as.numeric(exc_index$phi)
  for (t in exc_index$t) {
    pos <- pos + 1L
    if(first == 0) {
      first <- t
      prev <- t
      span_start_pos <- pos
    }
    if(t - prev > 1) {
      span_end_pos <- pos - 1L
      phi_seg <- phi_all[span_start_pos:span_end_pos]
      if (is.null(dom)) {
        lr <- c(NA_real_, NA_real_)
      } else if (dom$type == "cyclical") {
        lr <- spar_angle_range_from_transitions(phi_seg, dom)
      } else {
        lr <- spar_angle_range(phi_seg, dom)
      }
      First <- c(First, first)
      Last <- c(Last, prev)
      span_id <- c(span_id, id)
      n_points <- c(n_points, cnt)
      lphi <- c(lphi, as.numeric(lr[1]))
      uphi <- c(uphi, as.numeric(lr[2]))

      first <- t
      id <- id + 1
      cnt <- 0
      span_start_pos <- pos
    }

    prev <- t
    cnt <- cnt + 1
  }

  First <- c(First, first)
  Last <- c(Last, prev)
  span_id <- c(span_id, id)
  n_points <- c(n_points, cnt)
  phi_seg <- phi_all[span_start_pos:pos]
  if (is.null(dom)) {
    lr <- c(NA_real_, NA_real_)
  } else if (dom$type == "cyclical") {
    lr <- spar_angle_range_from_transitions(phi_seg, dom)
  } else {
    lr <- spar_angle_range(phi_seg, dom)
  }
  lphi <- c(lphi, as.numeric(lr[1]))
  uphi <- c(uphi, as.numeric(lr[2]))

  new_exceedance_spans(data.frame(span_id, First, Last, n_points, lphi, uphi))
}

new_exceedance_spans <- function(spans,
                                 t = NULL,
                                 metadata = list()) {
  stopifnot(is.data.frame(spans))
  stopifnot(all(c("First", "Last") %in% names(spans)))

  structure(
    list(
      spans = spans,
      t = t,
      metadata = metadata
    ),
    class = "exceedance_spans"
  )
}

extract_clustered_excursion_spans <- function(exc_spans, gap_rule, angle_domain = NULL) {
  spans <- exc_spans$spans

  if (nrow(spans) == 0L) {
    return(new_clustered_excursion_spans(
      data.frame(
        cluster_id = integer(0),
        First = numeric(0),
        Last = numeric(0),
        max_gap = numeric(0),
        total_gap = numeric(0),
        n_spans = integer(0),
        n_points = integer(0),
        total_span = numeric(0)
      ),
      gap_rule = gap_rule
    ))
  }

  spans$gap <- c(0, spans$First[-1] - spans$Last[-nrow(spans)] - 1)
  is_new_cluster <- c(TRUE, spans$gap[-1] > gap_rule)
  spans$cluster_id <- cumsum(is_new_cluster)

  cluster_groups <- split(spans, spans$cluster_id)

  max_gap <- c()
  total_gap <- c()
  n_spans <- c()
  n_points <- c()
  First <- c()
  Last <- c()
  total_span <- c()
  cluster_id <- c()
  for(i in seq_along(cluster_groups)){
    cluster <- cluster_groups[[i]]
    span_range <- min(cluster$span_id):max(cluster$span_id)
    max_gap <- c(max_gap,max(c(cluster$gap[-1], 0)))
    total_gap <- c(total_gap, sum(cluster$gap[-1]))
    n_spans <- c(n_spans, length(span_range))
    n_points <- c(n_points, sum(cluster$n_points))
    First <- c(First, cluster$First[1])
    Last <- c(Last, cluster$Last[nrow(cluster)])
    total_span <- c(total_span, cluster$Last[nrow(cluster)]-cluster$First[1]+1)
    cluster_id <- c(cluster_id, i)
  }
  span <- data.frame(cluster_id, First, Last, max_gap, total_gap, n_spans, n_points, total_span)

  support_ranges <- NULL
  if (!is.null(angle_domain) && all(c("lphi", "uphi") %in% names(spans))) {
    dom <- as_spar_angle_domain(angle_domain)
    rows <- list()
    row_id <- 0L
    for (i in seq_along(cluster_groups)) {
      cl <- cluster_groups[[i]]
      ranges <- lapply(seq_len(nrow(cl)), function(k) as.numeric(c(cl$lphi[k], cl$uphi[k])))
      ranges <- Filter(function(r) length(r) == 2L && all(is.finite(r)), ranges)
      if (length(ranges) == 0L) next

      merged <- if (dom$type == "cyclical") spar_angle_merge_range_list(ranges, dom) else ranges
      for (j in seq_along(merged)) {
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          cluster_id = i,
          support_id = j,
          lphi = as.numeric(merged[[j]][1]),
          uphi = as.numeric(merged[[j]][2])
        )
      }
    }
    if (length(rows) > 0L) {
      support_ranges <- do.call(rbind, rows)
    }
  }

  new_clustered_excursion_spans(
    span,
    gap_rule,
    metadata = list(support_ranges = support_ranges)
  )
}

new_clustered_excursion_spans <- function(spans,
                                          gap_rule,
                                          metadata = list()) {
  stopifnot(is.data.frame(spans))
  stopifnot(all(c("cluster_id", "First", "Last") %in% names(spans)))

  structure(
    list(
      spans = spans,
      gap_rule = gap_rule,
      metadata = metadata
    ),
    class = "clustered_excursion_spans"
  )
}
build_excursion_path_group <- function(data,
                                       clustered_spans,
                                       mapping = NULL,
                                       schema = default_excursion_schema(),
                                       u = NULL,
                                       t = NULL,
                                       radial_fun = L1_radial,
                                       angle_fun = L1_angle,
                                       metadata = list(),
                                       keep_prev_next = TRUE,
                                       ...) {
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(!is.null(clustered_spans))
  spans <- clustered_spans$spans
  stopifnot(is.data.frame(spans))
  stopifnot(all(c("cluster_id", "First", "Last") %in% names(spans)))

  # Resolve schema once for the whole group
  map_schema <- resolve_mapping(mapping)
  schema <- resolve_schema(schema = schema, mapping = map_schema)

  n_clusters <- nrow(spans)
  paths <- vector("list", n_clusters)

  for (i in seq_len(n_clusters)) {
    first_i <- spans$First[i]
    last_i  <- spans$Last[i]

    stopifnot(first_i >= 1, last_i <= nrow(data), first_i <= last_i)

    # --------------------------------------------------------
    # Slice cluster rows from original data
    # --------------------------------------------------------
    cluster_slice <- data[first_i:last_i, , drop = FALSE]

    # --------------------------------------------------------
    # Optional neighboring points for closure/interpolation
    # --------------------------------------------------------
    X_prev <- NULL
    X_next <- NULL

    if (keep_prev_next) {
      if (first_i > 1) {
        if (is.data.frame(data)) {
          X_prev <- as.numeric(data[first_i - 1, schema$X, drop = TRUE])
        } else {
          X_prev <- as.numeric(data[first_i - 1, schema$X, drop = FALSE])
        }
      }

      if (last_i < nrow(data)) {
        if (is.data.frame(data)) {
          X_next <- as.numeric(data[last_i + 1, schema$X, drop = TRUE])
        } else {
          X_next <- as.numeric(data[last_i + 1, schema$X, drop = FALSE])
        }
      }
    }

    # --------------------------------------------------------
    # Path-specific metadata
    # --------------------------------------------------------
    path_metadata <- c(
      list(
        cluster_id = spans$cluster_id[i],
        First = first_i,
        Last = last_i,
        X_prev = X_prev,
        X_next = X_next
      )
    )

    # Carry over any extra span summary columns into metadata
    extra_cols <- setdiff(names(spans), c("cluster_id", "First", "Last"))
    if (length(extra_cols) > 0) {
      for (nm in extra_cols) {
        path_metadata[[nm]] <- spans[[nm]][i]
      }
    }

    # --------------------------------------------------------
    # Build excursion_path object
    # --------------------------------------------------------
    path_i <- new_excursion_path(
      cluster_slice,
      mapping = mapping,
      schema = schema,
      u = u,
      t = t,
      radial_fun = radial_fun,
      angle_fun = angle_fun,
      metadata = path_metadata,
      ...
    )

    paths[[i]] <- path_i
  }

  # Name paths by cluster id
  names(paths) <- paste0("cluster_", spans$cluster_id)

  # Return grouped object
  new_excursion_path_group(
    paths = paths,
    spans = clustered_spans,
    schema = schema,
    metadata = c(
      metadata,
      list(
        u = u
      )
    )
  )
}
new_excursion_path_group <- function(paths,
                                     spans,
                                     schema = default_excursion_schema(),
                                     metadata = list()) {
  stopifnot(is.list(paths))
  stopifnot(inherits(spans, "clustered_excursion_spans"))

  structure(
    list(
      paths = paths,
      spans = spans,
      schema = schema,
      metadata = metadata
    ),
    class = "excursion_path_group"
  )
}
length.excursion_path_group <- function(x) {
  length(x$paths)
}
# ============================================================
# 2. Path object
#
# Stores:
#   - original-space points
#   - threshold u
#   - transform functions

ep_map <- function(...) {
  mapping <- rlang::enquos(...)
  class(mapping) <- "excursion_mapping"
  mapping
}

resolve_mapping <- function(mapping = NULL) {
  if (is.null(mapping)) return(NULL)

  out <- list()

  for (nm in names(mapping)) {
    expr <- rlang::get_expr(mapping[[nm]])

    # simple symbol: phi = phi
    if (rlang::is_symbol(expr)) {
      out[[nm]] <- rlang::as_string(expr)
      next
    }

    # string literal: phi = "phi"
    if (rlang::is_string(expr)) {
      out[[nm]] <- rlang::as_string(expr)
      next
    }

    # vector of strings already materialized
    if (is.character(expr)) {
      out[[nm]] <- expr
      next
    }

    # c(hs, tm2)
    if (rlang::is_call(expr, "c")) {
      args <- as.list(expr)[-1]

      out[[nm]] <- vapply(args, function(a) {
        if (rlang::is_symbol(a)) {
          rlang::as_string(a)
        } else if (rlang::is_string(a)) {
          rlang::as_string(a)
        } else {
          stop(
            "Unsupported element inside c(...) for mapping '", nm, "': ",
            rlang::expr_text(a)
          )
        }
      }, character(1))

      next
    }

    if (rlang::is_call(expr, "data.frame")) {
      args <- as.list(expr)[-1]  # drop function name

      arg_names <- names(args)
      is_named <- !is.null(arg_names) & nzchar(arg_names)

      positional <- args[!is_named]
      named <- args[is_named]

      out[[nm]] <- vapply(positional, rlang::as_name, character(1))

      if (length(named) > 0) {
        out[[paste0(nm, "_meta")]] <- lapply(named, function(x) {
          if (rlang::is_atomic(x)) {
            x
          } else {
            rlang::expr_text(x)
          }
        })
      }

      next
    }

    # matrix(a, b, ncol = 2)
    if (rlang::is_call(expr, "matrix")) {
      args <- as.list(expr)[-1]  # drop function name

      arg_names <- names(args)
      is_named <- !is.null(arg_names) & nzchar(arg_names)

      positional <- args[!is_named]
      named <- args[is_named]

      out[[nm]] <- vapply(positional, rlang::as_name, character(1))

      if (length(named) > 0) {
        out[[paste0(nm, "_meta")]] <- lapply(named, function(x) {
          if (rlang::is_atomic(x)) {
            x
          } else {
            rlang::expr_text(x)
          }
        })
      }

      next
    }

    stop(
      "Unsupported mapping expression for '", nm, "': ",
      rlang::expr_text(expr)
    )
  }

  out
}

resolve_mapped_object <- function(df, mapping) {
  if (is.null(mapping)) return(list())

  out <- list()

  for (nm in names(mapping)) {
    val <- rlang::eval_tidy(mapping[[nm]], df)

    if (is.character(val)) {
      out[[nm]] <- val
    } else if (is.numeric(val) || is.vector(val)) {
      out[[nm]] <- val
    } else {
      out[[nm]] <- val
    }
  }
  out
}

resolve_mapped_df <- function(df, mapping) {
  if (is.null(mapping)) return(data.frame())

  out <- list()

  for (nm in names(mapping)) {
    val <- rlang::eval_tidy(mapping[[nm]], df)

    if (is.character(val)) {
      out[[nm]] <- val
    } else if (is.numeric(val) || is.vector(val)) {
      out[[nm]] <- val
    } else {
      out[[nm]] <- val
    }
  }
  data.frame(out)
}

resolve_schema <- function(schema, mapping, defaults = default_excursion_schema()) {
  mapping <- resolve_mapping(mapping)

  out <- merge_schema_mapping(defaults, schema)
  out <- merge_schema_mapping(out, mapping)

  out
}

merge_schema_mapping <- function(schema, mapping) {
  if (is.null(mapping)) return(schema)

  out <- schema
  for (nm in names(mapping)) {
    map_value <- mapping[[nm]]
    out[[nm]] <- if (is.null(map_value)) schema[[nm]] else map_value
  }

  out
}

default_excursion_schema <- function() {
  list(
    X = c(1,2),
    radial_fun = "radial_fun",
    angle_fun = "angle_fun",
    u = "u",
    t = "t"
  )
}

eval_excess <- function(u, R, schema) {
  UseMethod("eval_excess")
}

eval_excess.data.frame <- function(u, R, schema) {

}

eval_excess.numeric <- function(u, R, schema) {
  dims <- dim(u)
  if (is.null(dims)) {
    if (length(u) == 1 || length(u) == length(R)) return (pmax(R-u,0))
    else stop ("Undefined behavior, length of thresholds not equal to one or the length of radial values. To provide ")
  } else {
    if(dims[2] == 1) {
      if (dims[1] == 1 || dims[1] == length(R)) return (pmax(R-u[,1],0))
      else stop ("Undefined behavior, length of thresholds not equal to one or the length of radial values. To provide angular map include an angle column.")
    } else if (dims[2] >= length(schema$phi) + 1){
      c_names <- colnames(u)
      if (dims[2] == length(schema$phi) + 1){
        phi_cols <- colnames()
        if(nrow(u) <= length(R)){
          stop ("Undefined behavior, matrix u with row count not equal to the observations not yet supported.")
          if (dims[2] == 2) {
            if(schema$phi %in% c_names) return (pmax(R-u[u[,c_names==schema$phi],c_names!=schema$phi]))
          }
        }
        else { ## nrow(u) == length(R)
          if (dims[2] == 2) {
            if(schema$phi %in% c_names) return (pmax(R-u[,c_names!=schema$phi],0))
            if(schema$R %in% c_names) return (pmax(R-u[,c_names==schema$R], 0))
          } else if(all(schema$phi %in% c_names)) return (pmax(R-u[,!(c_names %in% schema$phi)], 0))
        }
      }
    }
  }
  stop ("Undefined behavior. eval_excess.numeric failed to match input. Unknown error.")
}

validate_field_sources <- function(data,
                                   schema,
                                   ...#,
                                   #valid_fun = function(x) !(is.null(x) || (length(x) == 1 && is.na(x)))
                                    ) {
  args <- list(...)

  out <- setNames(vector("list", length(schema)), names(schema))

  for (nm in names(schema)) {
    schema_name <- schema[[nm]]

    as_column <- !is.null(schema_name) &&
      is.character(schema_name) &&
      length(schema_name) >= 1L &&
      all(schema_name %in% names(data)) #&&
      #valid_fun(data[schema_name])

    as_parameter <- nm %in% names(args)# && valid_fun(args[[nm]])

    out[[nm]] <- c(
      as_column = as_column,
      as_parameter = as_parameter
    )
  }

  out
}

resolve_field_name <- function(field) {
  expr <- rlang::enexpr(field)

  if (rlang::is_string(expr)) {
    return(rlang::as_string(expr))
  }

  if (rlang::is_symbol(expr)) {
    return(rlang::as_string(expr))
  }

  stop("field must be either a string or an unquoted name")
}

resolve_field_value <- function(field, data, schema, source_info, ...) {
  q <- rlang::enquo(field)
  field <- resolve_field_name(q)
  args <- list(...)

  if (isTRUE(source_info[[field]]["as_column"])) {
    schema_name <- schema[[field]]

    # x may be multi-column
    if (length(schema_name) > 1L) {
      return(data[, schema_name, drop = FALSE])
    }

    return(data[[schema_name]])
  }

  if (isTRUE(source_info[[field]]["as_parameter"])) {
    return(args[[field]])
  }

  NULL
}
resolve_field_value_chr <- function(field, data, schema, source_info, ...) {
  args <- list(...)

  if (isTRUE(source_info[[field]]["as_column"])) {
    schema_name <- schema[[field]]

    # x may be multi-column
    if (length(schema_name) > 1L) {
      return(data[, schema_name, drop = FALSE])
    }

    return(data[[schema_name]])
  }

  if (isTRUE(source_info[[field]]["as_parameter"])) {
    return(args[[field]])
  }

  NULL
}

require_sources <- function(source_info, fields) {
  ok <- vapply(source_info[fields], any, FUN.VALUE = logical(1))
  all(ok)
}

new_excursion_path <- function(data,
                               mapping = NULL,
                               schema = default_excursion_schema(),
                               ...,
                               u = NULL,
                               t = NULL,
                               excess = NULL,
                               radial_fun = L1_radial,
                               angle_fun = L1_angle,
                               metadata = list()) {
  UseMethod("new_excursion_path")
}

new_excursion_path.data.frame <- function(data,
                                          mapping = NULL,
                                          schema = default_excursion_schema(),
                                          ...,
                                          u = NULL,
                                          t = NULL,
                                          excess = NULL,
                                          radial_fun = L1_radial,
                                          angle_fun = L1_angle,
                                          metadata = list()) {
  stopifnot(nrow(data) >= 1)

  mapping <- resolve_mapping(mapping)
  def_schema <- default_excursion_schema()
  schema <- merge_schema_mapping(def_schema, schema)
  schema <- merge_schema_mapping(schema, mapping)

  src <- validate_field_sources(data, schema, u=u, t=t, excess=excess, radial_fun=radial_fun, angle_fun=angle_fun, metadata=metadata,...)

  if (!require_sources(src, names(def_schema))) {
    print(src)
    stop ("Missing required values in default schema.")
  }
  resolved <- list()
  for (nm in names(schema)){
    resolved[[nm]] <- resolve_field_value_chr(nm, data, schema, src, u=u, t=t, excess=excess, radial_fun=radial_fun, angle_fun=angle_fun, metadata=metadata, ...)
  }

  resolved$X <- as.matrix(resolved$X)
  storage.mode(resolved$X) <- "double"

  obs <- nrow(resolved$X)
  stopifnot(length(resolved$u) == 1 || length(resolved$u) == obs)
  if(is.null(resolved$phi) || (dim(resolved$phi)[1] %||% length(resolved$phi)) < obs) resolved$phi <- angle_fun(resolved$X)
  if(is.null(resolved$R) || length(resolved$R) < obs) resolved$R <- radial_fun(resolved$X)
  if(is.null(resolved$excess) || length(resolved$excess) < obs) resolved$excess <- pmax(resolved$R - resolved$u, 0)
  null_eval <- is.null(resolved$t)
  na_eval <- if(null_eval) FALSE else is.na(resolved$t)
  if(null_eval || length(resolved$t) != obs || any(na_eval)){
    if(!null_eval) warning(paste("Time parameter t has", sum(na_eval), "NA values out of", length(resolved$t), "for", obs,"observations."))
    row_names <- suppressWarnings(as.numeric(row.names(data)))
    row_names_valid <- length(row_names) == obs && !any(is.na(row_names))
    resolved$t <- if (row_names_valid){
      row_names
    } else {
      warning("Defaulting to sequential unit time distance observations.")
      seq_len(obs)
    }
  }


  stopifnot(
    length(resolved$R) == obs,
    (dim(resolved$phi)[1] %||% length(resolved$phi)) == obs,
    length(resolved$excess) == obs,
    length(resolved$t) == obs
  )

  structure(
    list(
      X = resolved$X,
      u = resolved$u,
      t = resolved$t,
      radial_fun = resolved$radial_fun,
      angle_fun = resolved$angle_fun,
      R = resolved$R,
      phi = resolved$phi,
      excess = resolved$excess,
      cluster_id = resolved$cluster_id,
      schema = schema,
      metadata = metadata
    ),
    class = "excursion_path"
  )
}

new_excursion_path.matrix <- function(data,
                                      mapping = NULL,
                                      schema = default_excursion_schema(),
                                      ...,
                                      u = NULL,
                                      t = NULL,
                                      excess = NULL,
                                      radial_fun = L1_radial,
                                      angle_fun = L1_angle,
                                      metadata = list()) {
  stopifnot(nrow(data) >= 1)

  mapping <- resolve_mapping(mapping)
  def_schema <- default_excursion_schema()
  schema <- merge_schema_mapping(def_schema, schema)
  schema <- merge_schema_mapping(schema, mapping)


  src <- validate_field_sources(data, schema, ..., u, t, excess, radial_fun, angle_fun, metadata)

  if (!require_sources(src, names(def_schema))) {
    stop ("Missing required values in default schema.")
  }
  resolved <- list()
  for (nm in names(schema)){
    resolved[[nm]] <- resolve_field_value(!!rlang::sym(nm), data, schema, src, ..., u, t, excess, radial_fun, angle_fun, metadata)
  }

  obs <- nrow(resolved$X)
  stopifnot(length(resolved$u) == 1 || length(resolved$u) == obs)
  if(is.null(resolved$phi) || (dim(resolved$phi)[1] %||% length(resolved$phi)) < obs) resolved$phi <- angle_fun(resolved$X)
  if(is.null(resolved$R) || length(resolved$R) < obs) resolved$R <- radial_fun(resolved$X)
  if(is.null(resolved$excess) || length(resolved$excess) < obs) resolved$excess <- pmax(resolved$R - resolved$u, 0)
  null_eval <- is.null(resolved$t)
  na_eval <- if(null_eval) FALSE else is.na(resolved$t)
  if(null_eval || length(resolved$t) != obs || any(na_eval)){
    if(!null_eval) warning(paste("Time parameter t has", sum(na_eval), "NA values out of", length(resolved$t), "for", obs,"observations."))
    row_names <- suppressWarnings(as.numeric(row.names(data)))
    row_names_valid <- length(row_names) == obs && !any(is.na(row_names))
    resolved$t <- if (row_names_valid){
      row_names
    } else {
      warning("Defaulting to sequential unit time distance observations.")
      seq_len(obs)
    }
  }


  stopifnot(
    length(resolved$R) == obs,
    (dim(resolved$phi)[1] %||% length(resolved$phi)) == obs,
    length(resolved$excess) == obs,
    length(resolved$t) == obs
  )

  structure(
    list(
      X = resolved$X,
      u = resolved$u,
      t = resolved$t,
      radial_fun = resolved$radial_fun,
      angle_fun = resolved$angle_fun,
      R = resolved$R,
      phi = resolved$phi,
      excess = resolved$excess,
      cluster_id = resolved$cluster_id,
      schema = schema,
      metadata = metadata
    ),
    class = "excursion_path"
  )
}

print.excursion_path <- function(x, ...) {
  cat("excursion_path\n")
  cat("  points    :", nrow(x$X), "\n")
  cat("  dimension :", ncol(x$X), "\n")
  cat("  threshold :", x$u, "\n")
  cat("  max excess:", max(x$excess, na.rm = TRUE), "\n")
  invisible(x)
}

# ============================================================
# 3. Crossing helper for original-space interpolation
#
# Find s in [0,1] such that:
#   g((1-s)Xa + s Xb) = u
# using bisection, assuming one endpoint is below and the other above.

find_threshold_crossing_original <- function(Xa, Xb, u, radial_fun,
                                             tol = 1e-8, max_iter = 100) {
  Xa <- as.matrix(Xa)
  Xb <- as.matrix(Xb)

  if (nrow(Xa) != 1L) Xa <- matrix(as.numeric(Xa), nrow = 1)
  if (nrow(Xb) != 1L) Xb <- matrix(as.numeric(Xb), nrow = 1)

  x_names <- colnames(Xa)
  if (is.null(x_names)) {
    x_names <- colnames(Xb)
  }
  f <- function(s) {
    xs <- (1 - s) * as.numeric(Xa[1, ]) + s * as.numeric(Xb[1, ])
    Xs <- matrix(xs, nrow = 1)
    if (!is.null(x_names)) colnames(Xs) <- x_names
    radial_fun(Xs) - u
  }

  fa <- f(0)
  fb <- f(1)

  if (is.na(fa) || is.na(fb)) return(NULL)
  if (fa == 0) return(0)
  if (fb == 0) return(1)
  if (fa * fb > 0) return(NULL)

  lo <- 0
  hi <- 1
  flo <- fa

  for (iter in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    fmid <- f(mid)

    if (abs(fmid) < tol || (hi - lo) < tol) return(mid)

    if (flo * fmid <= 0) {
      hi <- mid
    } else {
      lo <- mid
      flo <- fmid
    }
  }

  0.5 * (lo + hi)
}

# ============================================================
# 4. Segment builders
#
# Two interpolation modes:
#   original-space:
#      interpolate in X, then transform
#   angular-radial:
#      transform endpoints, then interpolate in (phi, Y)
#
# Returned segment objects are sampled densely enough that we can
# recover an envelope on a fixed angular grid robustly.

build_segment_original_space <- function(Xa, Xb, u, radial_fun, angle_fun,
                                         n_sub = 100L) {
  s <- seq(0, 1, length.out = n_sub)
  Xs <- vapply(
    s,
    function(si) as.numeric((1 - si) * Xa + si * Xb),
    numeric(length(Xa))
  )
  Xs <- t(Xs)

  tr <- transform_points(Xs, radial_fun, angle_fun, (1-s)*u[1]+s*u[2])

  data.frame(
    s = s,
    phi = tr$phi,
    Y = tr$Y
  )
}

build_segment_angular_radial <- function(phi_a, Y_a, phi_b, Y_b,
                                         n_sub = 100L) {
  s <- seq(0, 1, length.out = n_sub)
  data.frame(
    s = s,
    phi = (1 - s) * phi_a + s * phi_b,
    Y   = pmax((1 - s) * Y_a + s * Y_b, 0)
  )
}

# ============================================================
# 5. Closure points
#
# We allow closure in two ways:
#   A) explicit baseline angles phi0 / phi_end provided by user
#   B) infer from neighboring sub-threshold points if available
#
# For a skeleton, we support optional external points X_prev, X_next.
# If unavailable, we default to vertical closure in angle-space:
#   start baseline at phi_1, end baseline at phi_m
#
# Vertical closure is simple and stable, though less physically rich.

infer_closure_original <- function(X_prev, X1, X_last, X_next,
                                    u, radial_fun, angle_fun) {
  out <- list(phi0 = NULL, phi_end = NULL)

  X1 <- as.numeric(X1)
  X_last <- as.numeric(X_last)

  if (!is.null(X_prev)) {
    X_prev <- as.numeric(X_prev)
    s0 <- find_threshold_crossing_original(X_prev, X1, u[1], radial_fun)
    if (!is.null(s0)) {
      Xc <- matrix((1 - s0) * X_prev + s0 * X1, nrow = 1)
      out$phi0 <- angle_fun(Xc)
    }
  }

  if (!is.null(X_next)) {
    X_next <- as.numeric(X_next)
    s1 <- find_threshold_crossing_original(X_last, X_next, u[length(u)], radial_fun)
    if (!is.null(s1)) {
      Xc <- matrix((1 - s1) * X_last + s1 * X_next, nrow = 1)
      out$phi_end <- angle_fun(Xc)
    }
  }

  out
}
