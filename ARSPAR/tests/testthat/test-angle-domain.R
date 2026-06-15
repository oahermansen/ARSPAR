test_that("spar angle-domain utilities handle cyclical wrapping", {
  dom <- new_spar_angle_domain(type = "cyclical", lower = -2, upper = 2)

  x <- c(-2.2, -1.5, 1.9, 2.3)
  xn <- spar_normalize_angle(x, dom)
  expect_equal(xn, c(1.8, -1.5, 1.9, -1.7), tolerance = 1e-12)

  d <- spar_angle_delta(c(1.9, -1.9), phi0 = c(-1.9, 1.9), domain = dom)
  expect_equal(d, c(-0.2, 0.2), tolerance = 1e-12)

  span <- spar_angle_span(c(1.9, -1.9, -1.8), dom)
  expect_equal(span, 0.3, tolerance = 1e-12)
})

test_that("legacy angle_domain is adapted to spar_angle_domain", {
  legacy <- new_angle_domain(type = "circular", min = -2, max = 2, period = 4)
  dom <- as_spar_angle_domain(legacy)

  expect_s3_class(dom, "spar_angle_domain")
  expect_identical(dom$type, "cyclical")
  expect_equal(dom$lower, -2)
  expect_equal(dom$upper, 2)
})

test_that("angle extraction interpolation works with cyclical domain", {
  df <- data.frame(
    cluster_id = c(1L, 1L, 2L, 2L),
    phi = c(-1.9, 1.9, -1.0, 1.0),
    excess = c(1, 3, 2, 2),
    weight = c(1, 1, 1, 1)
  )

  sample <- new_envelope_sample(
    df,
    angle_domain = new_spar_angle_domain(type = "cyclical", lower = -2, upper = 2)
  )

  out <- extract_angle_distribution(sample, phi0 = -2, mode = "interpolate")

  expect_equal(nrow(out), 2)
  expect_true(all(is.finite(out$excess)))
})

test_that("smallest angular span handles cyclical seam crossing", {
  dom <- new_spar_angle_domain(type = "cyclical", lower = -2, upper = 2)
  sp <- spar_smallest_angular_span(1.9, -1.9, dom)

  expect_true(isTRUE(sp$crosses_seam))
  expect_identical(nrow(sp$segments), 2L)
  expect_equal(sp$delta, 0.2, tolerance = 1e-12)
})
