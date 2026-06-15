test_that("manual threshold assignment updates per-observation and excess", {
  X <- cbind(x = c(1, 2, 3), y = c(1, 1, 1))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) M[, 1] / rowSums(abs(M)),
    domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 1),
    source = "original",
    run = TRUE
  )

  obj <- spar_set_threshold(obj, u = 2.5, name = "fixed", set_active = TRUE, compute_excess = TRUE)

  expect_identical(obj$threshold$active, "fixed")
  expect_equal(obj$threshold$per_observation, rep(2.5, 3), tolerance = 1e-12)

  expected_excess <- pmax(obj$angular$R - 2.5, 0)
  expect_equal(obj$excess$value, expected_excess, tolerance = 1e-12)
  expect_identical(obj$excess$is_exceedance, expected_excess > 0)
})

test_that("functional threshold can use angular fields", {
  X <- cbind(x = c(1, 2, 3, 4), y = c(2, 2, 2, 2))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) M[, 1] / rowSums(abs(M)),
    domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 1),
    source = "original",
    run = TRUE
  )

  obj <- spar_set_threshold(
    obj,
    u = function(R, phi) 0.2 * R + 2 * phi,
    name = "phi_adaptive",
    set_active = TRUE,
    compute_excess = TRUE
  )

  expected_u <- 0.2 * obj$angular$R + 2 * obj$angular$phi
  expect_equal(obj$threshold$per_observation, expected_u, tolerance = 1e-12)
  expect_equal(obj$excess$value, pmax(obj$angular$R - expected_u, 0), tolerance = 1e-12)
})

test_that("registered estimator can be applied", {
  X <- cbind(x = c(1, 2, 4), y = c(1, 3, 2))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) atan2(M[, 2], M[, 1]),
    domain = new_spar_angle_domain(type = "cyclical", lower = -pi, upper = pi),
    source = "original",
    run = TRUE
  )

  est <- list(offset = 1.0)
  pred <- function(estimator, x) estimator$offset + abs(x$angular$phi)

  obj <- spar_register_threshold_estimator(
    obj,
    name = "custom",
    estimator = est,
    predict_fun = pred,
    set_active = TRUE
  )

  obj <- spar_apply_threshold(obj, compute_excess = TRUE)

  expected_u <- 1 + abs(obj$angular$phi)
  expect_equal(obj$threshold$per_observation, expected_u, tolerance = 1e-12)
  expect_equal(obj$excess$value, pmax(obj$angular$R - expected_u, 0), tolerance = 1e-12)
})

test_that("evgam ALD threshold fitting registers estimator when available", {
  testthat::skip_if_not_installed("evgam")

  set.seed(20260317)
  X <- cbind(x = rnorm(300), y = rnorm(300))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) atan2(M[, 2], M[, 1]),
    domain = new_spar_angle_domain(type = "cyclical", lower = -pi, upper = pi),
    source = "original",
    run = TRUE
  )

  obj <- spar_fit_threshold_ald_evgam(
    obj,
    name = "ald",
    formula = stats::as.formula("R ~ s(phi, k = 8)"),
    tau = 0.9,
    set_active = TRUE
  )

  obj <- spar_apply_threshold(obj)

  expect_true(all(is.finite(obj$threshold$per_observation)))
  expect_identical(length(obj$threshold$per_observation), nrow(X))
})

test_that("evgam ALD defaults use cyclic spline for cyclical domains", {
  testthat::skip_if_not_installed("evgam")

  set.seed(20260317)
  X <- cbind(x = rnorm(200), y = rnorm(200))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) atan2(M[, 2], M[, 1]),
    domain = new_spar_angle_domain(type = "cyclical", lower = -pi, upper = pi),
    source = "original",
    run = TRUE
  )

  obj <- spar_fit_threshold_ald_evgam(
    obj,
    name = "ald_default",
    formula = NULL,
    k = 9,
    tau = 0.9,
    set_active = TRUE
  )

  ftxt <- paste(deparse(obj$threshold$estimators$ald_default$formula[[1]]), collapse = "")
  expect_match(ftxt, "bs = \"cc\"")
  expect_match(ftxt, "k = 9")

  obj <- spar_apply_threshold(obj)
  expect_true(all(is.finite(obj$threshold$per_observation)))
  expect_identical(length(obj$threshold$per_observation), nrow(X))

  grid <- spar_predict_threshold(obj, n = 120, name = "ald_default")
  expect_identical(nrow(grid), 120L)
  expect_true(all(is.finite(grid$threshold)))
  expect_true(all(grid$phi >= -pi & grid$phi <= pi))
})

test_that("high-level threshold fit wrapper works with evgam", {
  testthat::skip_if_not_installed("evgam")

  set.seed(20260317)
  X <- cbind(x = rnorm(180), y = rnorm(180))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) atan2(M[, 2], M[, 1]),
    domain = new_spar_angle_domain(type = "cyclical", lower = -pi, upper = pi),
    source = "original",
    run = TRUE
  )

  obj <- spar_fit_threshold(
    obj,
    method = "evgam_ald",
    name = "ald_wrap",
    tau = 0.8,
    k = 8,
    trace = 0,
    verbose = FALSE,
    apply = TRUE,
    compute_excess = TRUE
  )

  expect_identical(obj$threshold$active, "ald_wrap")
  expect_identical(length(obj$threshold$per_observation), nrow(X))
  expect_true(all(is.finite(obj$threshold$per_observation)))
  expect_identical(length(obj$excess$value), nrow(X))
})
