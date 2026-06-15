test_that("spar exceedance extraction works from representation state", {
  X <- cbind(a = 1:20, b = rep(1, 20))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) M[, 1] / rowSums(abs(M)),
    domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 1),
    source = "original",
    run = TRUE
  )

  obj <- spar_set_threshold(obj, u = 15, name = "fixed", set_active = TRUE, compute_excess = TRUE)

  idx <- spar_extract_exceedance_index(obj)

  expect_s3_class(idx, "exceedance_index")
  expect_true(length(idx$t) > 0L)
  expect_true(all(idx$excess > 0))
  expect_identical(length(idx$R), length(idx$phi))
})

test_that("spar excursion path-group builder works and can store outputs", {
  X <- cbind(a = 1:20, b = rep(1, 20))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) M[, 1] / rowSums(abs(M)),
    domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 1),
    source = "original",
    run = TRUE
  )

  obj <- spar_set_threshold(obj, u = 15, name = "fixed", set_active = TRUE, compute_excess = TRUE)

  grp <- spar_build_excursion_path_group(obj, gap_rule = 1, space = "original", store = FALSE)
  expect_s3_class(grp, "excursion_path_group")
  expect_true(length(grp$paths) >= 1L)

  obj2 <- spar_build_excursion_path_group(obj, gap_rule = 1, space = "original", store = TRUE)
  expect_s3_class(obj2, "spar_representation")
  expect_s3_class(obj2$excursions$paths, "excursion_path_group")
  expect_s3_class(obj2$excursions$clusters, "clustered_excursion_spans")
  expect_s3_class(obj2$excursions$pointwise, "exceedance_index")
})

test_that("declustering and envelope diagnostics work on representation", {
  set.seed(20260317)
  X <- cbind(a = rnorm(400), b = rnorm(400))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) atan2(M[, 2], M[, 1]),
    domain = new_spar_angle_domain(type = "cyclical", lower = -pi, upper = pi),
    source = "original",
    run = TRUE
  )

  obj <- spar_set_threshold(obj, u = quantile(obj$angular$R, 0.8), name = "fixed", set_active = TRUE)
  obj <- spar_decluster_excursions(obj, gap_rule = 2, store = TRUE)

  expect_s3_class(obj$excursions$pointwise, "exceedance_index")
  expect_s3_class(obj$excursions$clusters, "clustered_excursion_spans")
  expect_true(all(c("lphi", "uphi") %in% names(extract_exceedance_spans(obj$excursions$pointwise, angle_domain = obj$angular$domain)$spans)))
  expect_true("support_ranges" %in% names(obj$excursions$clusters$metadata))

  obj <- spar_build_excursion_path_group(obj, gap_rule = 2, store = TRUE)
  expect_s3_class(obj$excursions$paths, "excursion_path_group")

  obj <- spar_build_excursion_envelope_sample(obj, n_phi = 120, lambda = 0, store = TRUE)
  expect_true(is.list(obj$excursions$envelopes))
  expect_s3_class(obj$excursions$envelopes$sample, "envelope_sample")

  cal <- spar_compute_threshold_calibration(obj, store = FALSE)
  expect_s3_class(cal, "threshold_calibration")
  expect_true(length(cal$phi_grid) >= 1L)

  ad <- spar_compute_angular_density(obj, sample = obj$excursions$envelopes$sample, time_span = nrow(X), store = FALSE)
  expect_s3_class(ad, "angular_density_diagnostic")
  expect_true(length(ad$phi_grid) >= 2L)
  expect_true(all(is.finite(ad$declustered_rate)))

  reg_decl <- spar_extract_region_declustered_maxima(
    obj,
    region_phi = c(-0.5, 0.5),
    gap_rule = 2
  )
  expect_true(is.data.frame(reg_decl))
  expect_true(all(c("cluster_id", "max_excess") %in% names(reg_decl)))

  reg_env <- spar_extract_region_envelope_maxima(
    obj,
    region_phi = c(-0.5, 0.5),
    sample = obj$excursions$envelopes$sample
  )
  expect_true(is.data.frame(reg_env))
  expect_true(all(c("cluster_id", "max_excess") %in% names(reg_env)))

  obj <- spar_build_upper_excursion_paths(obj, store = TRUE)
  expect_true(is.list(obj$excursions$upper_paths$paths))
  expect_true(length(obj$excursions$upper_paths$paths) >= 1L)

  reg_upper <- spar_extract_region_upper_path_maxima(
    obj,
    region_phi = c(-0.5, 0.5),
    upper_paths = obj$excursions$upper_paths
  )
  expect_true(is.data.frame(reg_upper))
  expect_true(all(c("cluster_id", "max_excess") %in% names(reg_upper)))
})

test_that("envelope polyline update respects cyclical seam", {
  phi_grid <- seq(-2, 2, length.out = 201)
  H <- rep(0, length(phi_grid))
  dom <- new_spar_angle_domain(type = "cyclical", lower = -2, upper = 2)

  H2 <- update_envelope_from_polyline(
    H = H,
    phi_grid = phi_grid,
    phi_vec = c(1.9, -1.9),
    Y_vec = c(1, 1),
    angle_domain = dom
  )

  idx_pos <- which(H2 > 0)
  expect_true(length(idx_pos) > 0)

  phi_pos <- phi_grid[idx_pos]
  expect_true(any(phi_pos <= -1.8))
  expect_true(any(phi_pos >= 1.8))
  expect_false(any(phi_pos > -1.5 & phi_pos < 1.5))
})

test_that("optimized envelope method matches brute-force and fast-paths by lambda", {
  X <- cbind(a = 1:40, b = rep(1, 40))
  obj <- spar_representation(X)

  obj <- spar_set_angular_map(
    obj,
    radial_fun = function(M) rowSums(abs(M)),
    angle_fun = function(M) M[, 1] / rowSums(abs(M)),
    domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 1),
    source = "original",
    run = TRUE
  )

  obj <- spar_set_threshold(obj, u = 15, name = "fixed", set_active = TRUE, compute_excess = TRUE)
  grp <- spar_build_excursion_path_group(obj, gap_rule = 1, space = "original", store = FALSE)
  path <- grp$paths[[1]]

  dom <- new_spar_angle_domain(type = "cyclical", lower = -2, upper = 2)
  phi_grid <- spar_angle_grid(n = 160, domain = dom)

  env_bf0 <- build_cluster_envelopes(
    path = path,
    phi_grid = phi_grid,
    lambda = 0,
    angle_domain = dom,
    n_sub = 120,
    method = "brute_force"
  )

  env_opt0 <- build_cluster_envelopes(
    path = path,
    phi_grid = phi_grid,
    lambda = 0,
    angle_domain = dom,
    n_sub = 120,
    method = "optimized"
  )

  expect_null(env_opt0$angular_radial)
  expect_s3_class(env_opt0$original, "excursion_envelope")
  expect_equal(env_opt0$mixed$H, env_bf0$mixed$H, tolerance = 1e-10)

  env_bf1 <- build_cluster_envelopes(
    path = path,
    phi_grid = phi_grid,
    lambda = 1,
    angle_domain = dom,
    n_sub = 120,
    method = "brute_force"
  )

  env_opt1 <- build_cluster_envelopes(
    path = path,
    phi_grid = phi_grid,
    lambda = 1,
    angle_domain = dom,
    n_sub = 120,
    method = "optimized"
  )

  expect_null(env_opt1$original)
  expect_s3_class(env_opt1$angular_radial, "excursion_envelope")
  expect_equal(env_opt1$mixed$H, env_bf1$mixed$H, tolerance = 1e-10)
})

test_that("upper excursion paths use path transitions instead of angle ordering", {
  path <- structure(
    list(
      X = matrix(c(1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE),
      t = 1:3,
      phi = c(0, 1, 2),
      excess = c(2, 1, 2),
      R = c(2, 1, 2),
      metadata = list(cluster_id = 1L)
    ),
    class = "excursion_path"
  )

  up <- extract_upper_excursion_path(
    path,
    angle_domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 2),
    interpolation = "angular"
  )

  expect_equal(up$idx, 1:3)
})

test_that("upper excursion angular interpolation only uses actual transitions", {
  path <- structure(
    list(
      X = matrix(c(1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE),
      t = 1:3,
      phi = c(0, 2, 1),
      excess = c(2, 2, 1),
      R = c(2, 2, 1),
      metadata = list(cluster_id = 1L)
    ),
    class = "excursion_path"
  )

  up <- extract_upper_excursion_path(
    path,
    angle_domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 2),
    interpolation = "angular"
  )

  expect_equal(up$idx, c(1L, 2L))
})

test_that("upper excursion transformed interpolation uses bivariate ray intersections", {
  path <- structure(
    list(
      X = matrix(c(2, 0, 0, 2, 0.5, 0.5), ncol = 2, byrow = TRUE),
      t = 1:3,
      phi = c(0, 1, 0.5),
      excess = c(2, 2, 1),
      R = c(2, 2, 1),
      metadata = list(cluster_id = 1L)
    ),
    class = "excursion_path"
  )

  up <- extract_upper_excursion_path(
    path,
    angle_domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 1),
    interpolation = "transformed"
  )

  expect_equal(up$idx, c(1L, 2L))
})

test_that("dominated upper excursion observations still contribute transitions", {
  path <- structure(
    list(
      X = matrix(c(2, 0, 0.5, 0.5, 0, 2, 2, 0, 0.25, 0.5), ncol = 2, byrow = TRUE),
      t = 1:5,
      phi = c(0, 0.5, 1, 0, 0.75),
      excess = c(2, 1, 2, 2, 0.5),
      R = c(2, 1, 2, 2, 0.5),
      metadata = list(cluster_id = 1L)
    ),
    class = "excursion_path"
  )

  up <- extract_upper_excursion_path(
    path,
    angle_domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 1),
    interpolation = "angular"
  )

  expect_equal(up$idx, c(1L, 3L, 4L))
})

test_that("upper excursion paths do not bridge non-consecutive spans", {
  path <- structure(
    list(
      X = matrix(c(1, 0, 0, 1, 0.5, 0.5), ncol = 2, byrow = TRUE),
      t = c(1, 3, 4),
      phi = c(0, 2, 1),
      excess = c(2, 2, 1),
      R = c(2, 2, 1),
      metadata = list(cluster_id = 1L)
    ),
    class = "excursion_path"
  )

  up <- extract_upper_excursion_path(
    path,
    angle_domain = new_spar_angle_domain(type = "interval", lower = 0, upper = 2),
    interpolation = "angular"
  )

  expect_equal(up$idx, 1:3)
})
