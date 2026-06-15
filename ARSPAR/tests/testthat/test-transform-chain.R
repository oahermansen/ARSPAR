test_that("compiled transform step enforces output length invariants", {
  ctx <- new_spar_transform_context(
    current = data.frame(x = 1:4, check.names = FALSE),
    original = data.frame(x = 1:4, check.names = FALSE),
    steps = list()
  )

  step <- spar_compile_transform_step(
    exprs = list(x = quote(x + 1)),
    name = "ok",
    context = ctx
  )

  out <- spar_eval_transform_step(step, current = ctx$current)
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(4L, 1L))
  expect_identical(colnames(out), "x")

  bad_step <- spar_compile_transform_step(
    exprs = list(x = quote(c(1, 2))),
    name = "bad",
    context = ctx
  )

  expect_error(
    spar_eval_transform_step(bad_step, current = ctx$current),
    "has length"
  )
})

test_that("step diagnostics report runtime dependencies", {
  X <- cbind(x = 1:5, y = 6:10)

  chain <- spar_build_transform_chain(
    X_original = X,
    spar_step_mutate(
      x = x - mean(x),
      y = y - mean(y),
      .name = "center"
    ),
    spar_step_mutate(
      z = x + .y,
      .name = "combine"
    ),
    name = "demo"
  )

  step <- chain$steps[[2]]
  diag <- spar_transform_diagnostics(step)

  expect_s3_class(diag, "spar_transform_diagnostics")
  expect_identical(diag$name, "combine")
  expect_false(isTRUE(diag$needs_step_cache))
  expect_identical(diag$per_output[[1]]$runtime_steps, list())
  expect_identical(diag$per_output[[1]]$runtime_current, "x")
  expect_identical(diag$per_output[[1]]$runtime_original, "y")
})

test_that("compile chain ignores removed input_space spec field", {
  X <- cbind(x = 1:3)

  chain <- spar_compile_transform_chain(
    step_specs = list(
      list(
        name = "s1",
        exprs = list(x = quote(x + 1)),
        input_space = "original"
      )
    ),
    X_original = X,
    name = "legacy-spec"
  )

  expect_s3_class(chain, "spar_transform_chain")
  out <- spar_run_transform_chain(chain, X)
  expect_equal(drop(out$current), c(2, 3, 4))
})

test_that("multi-step chain matches manual computation on long random multivariate data", {
  set.seed(20260317)

  n <- 3000L
  i <- seq_len(n)

  a <- rnorm(n, mean = sin(i / 150), sd = 0.6 + i / n)
  b <- rexp(n, rate = 0.4 + i / (3 * n))
  c <- rnorm(n, mean = i / n, sd = 1 + 0.2 * cos(i / 80))

  X <- cbind(a = a, b = b, c = c)

  chain <- spar_build_transform_chain(
    X_original = X,
    spar_step_mutate(
      c1 = (a - mean(a)) / sd(a),
      c2 = log1p(b),
      c3 = c - median(c),
      c4 = a + b + c,
      .name = "s1"
    ),
    spar_step_mutate(
      d1 = c1 + 0.5 * .a,
      d2 = c2 - mean(c2) + .c,
      d3 = mean(c4) + sd(.b),
      .name = "s2"
    ),
    spar_step_mutate(
      e1 = d1 + d2 + d3,
      e2 = d1 - mean(d1) + s1.c1,
      .name = "s3"
    ),
    name = "long-random"
  )

  out <- spar_run_transform_chain(chain, X)

  c1 <- (a - mean(a)) / sd(a)
  c2 <- log1p(b)
  c3 <- c - median(c)
  c4 <- a + b + c

  d1 <- c1 + 0.5 * a
  d2 <- c2 - mean(c2) + c
  d3 <- rep(mean(c4) + sd(b), n)

  e1 <- d1 + d2 + d3
  e2 <- d1 - mean(d1) + c1

  expected <- cbind(e1 = e1, e2 = e2)

  expect_identical(colnames(out$current), c("e1", "e2"))
  expect_equal(out$current, expected, tolerance = 1e-12)

  expect_identical(names(chain$cache_requirements$forward), c("s1", "s2", "s3"))
  expect_identical(chain$cache_requirements$forward$s1, c("c1"))
  expect_identical(chain$cache_requirements$forward$s2, character(0))
  expect_identical(chain$cache_requirements$forward$s3, character(0))

  expect_identical(names(out$step_cache$s1), c("c1"))
  expect_identical(out$step_cache$s2, list())
  expect_identical(out$step_cache$s3, list())

  diag_s2 <- spar_transform_diagnostics(chain$steps[[2]])
  expect_identical(diag_s2$per_output[[3]]$output, "d3")
  expect_identical(diag_s2$per_output[[3]]$runtime_steps, list())
  expect_identical(diag_s2$per_output[[3]]$captured_count, 2L)
})

test_that("summary on step-qualified field is captured while reusing field names", {
  set.seed(20260317)

  n <- 2500L
  i <- seq_len(n)

  a <- rnorm(n, mean = cos(i / 90), sd = 0.4 + i / (2 * n))
  b <- rexp(n, rate = 0.5 + i / (4 * n))
  c <- rnorm(n, mean = i / (2 * n), sd = 0.8 + 0.1 * sin(i / 50))

  X <- cbind(a = a, b = b, c = c)

  chain <- spar_build_transform_chain(
    X_original = X,
    spar_step_mutate(
      c1 = (a - mean(a)) / sd(a),
      c2 = log1p(b),
      c3 = c - median(c),
      c4 = a + b + c,
      .name = "s1"
    ),
    spar_step_mutate(
      c1 = c1 + 0.1 * .a,
      c2 = c2 - mean(c2),
      c3 = c3 + .b,
      c4 = c4 + c1,
      .name = "s2"
    ),
    spar_step_mutate(
      c1 = c1 + s1.c1,
      c2 = c2 + sd(s2.c4),
      c3 = c3 + mean(s1.c3),
      c4 = c4,
      .name = "s3"
    ),
    name = "same-field-names"
  )

  out <- spar_run_transform_chain(chain, X)

  s1_c1 <- (a - mean(a)) / sd(a)
  s1_c2 <- log1p(b)
  s1_c3 <- c - median(c)
  s1_c4 <- a + b + c

  s2_c1 <- s1_c1 + 0.1 * a
  s2_c2 <- s1_c2 - mean(s1_c2)
  s2_c3 <- s1_c3 + b
  s2_c4 <- s1_c4 + s1_c1

  s3_c1 <- s2_c1 + s1_c1
  s3_c2 <- s2_c2 + sd(s2_c4)
  s3_c3 <- s2_c3 + mean(s1_c3)
  s3_c4 <- s2_c4

  expected <- cbind(c1 = s3_c1, c2 = s3_c2, c3 = s3_c3, c4 = s3_c4)

  expect_equal(out$current, expected, tolerance = 1e-12)

  expect_identical(chain$cache_requirements$forward$s1, c("c1"))
  expect_identical(chain$cache_requirements$forward$s2, character(0))
  expect_identical(chain$cache_requirements$forward$s3, character(0))

  expect_identical(names(out$step_cache$s1), c("c1"))
  expect_identical(out$step_cache$s2, list())
  expect_identical(out$step_cache$s3, list())

  expr_c2 <- chain$steps[[3]]$exprs$c2
  expect_identical(expr_c2$runtime_dependencies$steps, list())
  expect_identical(length(expr_c2$params), 1L)
  expect_equal(unname(unlist(expr_c2$params)), sd(s2_c4), tolerance = 1e-12)
})

test_that(".data reference supports identity and row-wise matrix operations", {
  set.seed(20260317)

  n <- 1500L
  x <- rnorm(n, mean = 0.2, sd = 1.1)
  y <- rexp(n, rate = 0.7)
  X <- cbind(x = x, y = y)

  chain <- spar_build_transform_chain(
    X_original = X,
    spar_step_mutate(
      .data = .data,
      .name = "s1"
    ),
    spar_step_mutate(
      .data = .data / rowSums(.data),
      .name = "s2"
    ),
    name = "row-ref"
  )

  out <- spar_run_transform_chain(chain, X)

  s1 <- X
  expected <- s1 / rowSums(s1)

  expect_equal(out$current, expected, tolerance = 1e-12)
  expect_identical(colnames(out$current), c("x", "y"))

  deps <- chain$steps[[1]]$exprs[[1]]$runtime_dependencies
  expect_true(isTRUE(deps$data))
  expect_identical(deps$steps, list())
})
