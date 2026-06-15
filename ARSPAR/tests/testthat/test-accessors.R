test_that("spar_data returns original matrix with schema names", {
  X <- cbind(a = 1:3, b = 4:6)
  obj <- spar_representation(X_original = X)

  out <- spar_data(obj, space = "original", format = "matrix",
                   include_time = FALSE, include_observation_id = FALSE)

  expect_true(is.matrix(out))
  expect_identical(colnames(out), c("a", "b"))
})

test_that("spar_data returns angular names from schema", {
  X <- cbind(Hs = 1:3, Tm = 4:6)
  obj <- spar_representation(X_original = X, validate = FALSE)
  obj$angular$R <- c(1, 2, 3)
  obj$angular$phi <- c(0.1, 0.2, 0.3)
  obj$schema$radial_name <- "radius"
  obj$schema$angle_name <- "theta"

  out <- spar_data(obj, space = "angular", format = "data.frame",
                   include_time = FALSE, include_observation_id = FALSE)

  expect_identical(names(out), c("radius", "theta"))
})

test_that("spar_data includes schema-respecting metadata names", {
  X <- cbind(Hs = 1:3, Tm = 4:6)
  tt <- as.POSIXct(c("2026-01-01 00:00:00", "2026-01-01 01:00:00", "2026-01-01 02:00:00"), tz = "UTC")

  obj <- spar_representation(X_original = X, time = tt)
  obj$schema$observation_id_name <- "obs_id"
  obj$schema$time_name <- "timestamp"

  out <- spar_data(obj, format = "data.frame")

  expect_identical(names(out), c("obs_id", "timestamp", "Hs", "Tm"))
})
