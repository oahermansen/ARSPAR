test_that("angular mutators compute and store representation fields", {
  X <- cbind(a = c(1, 2, 3), b = c(2, 1, 4))
  obj <- spar_representation(X)

  radial_fun <- function(M) rowSums(abs(M))
  angle_fun <- function(M) M[, 1] / rowSums(abs(M))
  dom <- new_spar_angle_domain(type = "interval", lower = 0, upper = 1)

  obj <- spar_set_angular_map(
    x = obj,
    radial_fun = radial_fun,
    angle_fun = angle_fun,
    domain = dom,
    name = "l1-ratio",
    source = "original",
    run = TRUE
  )

  expected_R <- c(3, 3, 7)
  expected_phi <- c(1 / 3, 2 / 3, 3 / 7)

  expect_equal(obj$angular$R, expected_R, tolerance = 1e-12)
  expect_equal(obj$angular$phi, expected_phi, tolerance = 1e-12)
  expect_identical(obj$angular$active, "l1-ratio")
  expect_identical(obj$angular$source, "original")

  A <- spar_data(obj, space = "angular", format = "matrix")
  expect_equal(A[, 1], expected_R, tolerance = 1e-12)
  expect_equal(A[, 2], expected_phi, tolerance = 1e-12)
})

test_that("angular mutators can use transformed source", {
  X <- cbind(x = c(1, 2, 4), y = c(2, 3, 5))
  obj <- spar_representation(X)

  obj <- spar_build_representation_transform(
    obj,
    spar_step_mutate(
      x = x - mean(x),
      y = y - mean(y),
      .name = "center"
    ),
    run = TRUE
  )

  radial_fun <- function(M) sqrt(rowSums(M^2))
  angle_fun <- function(M) atan2(M[, 2], M[, 1])

  obj <- spar_set_angular_map(
    x = obj,
    radial_fun = radial_fun,
    angle_fun = angle_fun,
    domain = new_spar_angle_domain(type = "cyclical", lower = -pi, upper = pi),
    source = "transformed",
    run = FALSE
  )

  obj <- spar_apply_angular_map(obj)

  Xt <- obj$transform$transformed_data
  expected_R <- sqrt(rowSums(Xt^2))
  expected_phi <- atan2(Xt[, 2], Xt[, 1])

  expect_equal(obj$angular$R, expected_R, tolerance = 1e-12)
  expect_equal(obj$angular$phi, expected_phi, tolerance = 1e-12)
  expect_identical(obj$angular$source, "transformed")
})
