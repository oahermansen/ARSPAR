test_that("representation can attach and run transform chain", {
  testthat::skip_if_not_installed("ellipsis", "0.3.2")

  set.seed(20260317)

  n <- 1200L
  i <- seq_len(n)
  x1 <- rnorm(n, mean = i / n, sd = 0.7 + i / (3 * n))
  x2 <- rexp(n, rate = 0.6 + i / (5 * n))
  x3 <- rnorm(n, mean = sin(i / 60), sd = 0.9)

  X <- cbind(x1 = x1, x2 = x2, x3 = x3)
  obj <- spar_representation(X)

  chain <- spar_build_transform_chain(
    X_original = X,
    spar_step_mutate(
      x1 = (x1 - mean(x1)) / sd(x1),
      x2 = log1p(x2),
      x3 = x3 - median(x3),
      .name = "s1"
    ),
    spar_step_mutate(
      x1 = x1 + .x1,
      x2 = x2 + sd(x3),
      x3 = x3,
      .name = "s2"
    ),
    name = "repr-chain"
  )

  obj <- spar_set_transform_chain(
    x = obj,
    chain = chain,
    run = TRUE,
    keep_step_cache = TRUE
  )

  s1_x1 <- (x1 - mean(x1)) / sd(x1)
  s1_x2 <- log1p(x2)
  s1_x3 <- x3 - median(x3)

  expected <- cbind(
    x1 = s1_x1 + x1,
    x2 = s1_x2 + sd(s1_x3),
    x3 = s1_x3
  )

  expect_equal(obj$transform$transformed_data, expected, tolerance = 1e-12)
  expect_identical(obj$schema$transformed_names, c("x1", "x2", "x3"))

  transformed <- spar_data(obj, space = "transformed", format = "matrix")
  expect_equal(transformed, expected, tolerance = 1e-12)

  expect_length(obj$transform$steps, 2)
  expect_s3_class(obj$transform$chain, "spar_transform_chain")
  expect_identical(obj$transform$step_cache$s1, list())
  expect_identical(obj$transform$step_cache$s2, list())
})

test_that("representation default transform is identity", {
  testthat::skip_if_not_installed("ellipsis", "0.3.2")

  X <- cbind(a = rnorm(50), b = rexp(50, 0.8))
  obj <- spar_representation(X)

  expect_equal(obj$transform$transformed_data, X, tolerance = 1e-12)
  expect_identical(obj$schema$transformed_names, c("a", "b"))

  Xt <- spar_data(obj, space = "transformed", format = "matrix")
  expect_equal(Xt, X, tolerance = 1e-12)
})

test_that("representation convenience builder compiles and runs chain", {
  testthat::skip_if_not_installed("ellipsis", "0.3.2")

  X <- cbind(a = 1:6, b = 11:16)
  obj <- spar_representation(X)

  obj <- spar_build_representation_transform(
    obj,
    spar_step_mutate(
      a = a - mean(a),
      b = b - mean(b),
      .name = "center"
    ),
    spar_step_mutate(
      a = a + .a,
      b = b,
      .name = "combine"
    ),
    name = "inline",
    run = TRUE
  )

  expected <- cbind(
    a = ((X[, "a"] - mean(X[, "a"])) + X[, "a"]),
    b = (X[, "b"] - mean(X[, "b"]))
  )

  expect_equal(obj$transform$transformed_data, expected, tolerance = 1e-12)
  expect_identical(obj$schema$transformed_names, c("a", "b"))
})
