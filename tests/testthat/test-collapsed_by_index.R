test_that("return correct collapsed measures w/o 1 or nfreq + 1", {
  spec <- matrix(c(
    rep(rep(1:2, each = 10), 5), rep(rep(3:4, each = 10), 5)
  ), nrow = 20, ncol = 10)
  endpoints_index <- 11
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index), expected)
})

test_that("return correct collapsed measures w/o 1 but with nfreq + 1", {
  spec <- matrix(c(
    rep(rep(1:2, each = 10), 5), rep(rep(3:4, each = 10), 5)
  ), nrow = 20, ncol = 10)
  endpoints_index <- c(11, 21)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index), expected)
})

test_that("return correct collapsed measures w/0 nfreq + 1 but with 1", {
  spec <- matrix(c(
    rep(rep(1:2, each = 10), 5), rep(rep(3:4, each = 10), 5)
  ), nrow = 20, ncol = 10)
  endpoints_index <- c(1, 11)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index), expected)
})

test_that("return correct collapsed measures w 1 and nfreq + 1", {
  spec <- matrix(c(
    rep(rep(1:2, each = 10), 5), rep(rep(3:4, each = 10), 5)
  ), nrow = 20, ncol = 10)
  endpoints_index <- c(1, 11, 21)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index), expected)
})

test_that("catch lower endpoints out of bounds", {
  spec <- matrix(rnorm(1000 * 10), nrow = 1000)
  endpoints_index <- c(0, 300, 500, 700, 1000)
  expect_error(rep_collapsed_by_index(spec, endpoints_index),
               "endpoints must be between")
})

test_that("catch upper endpoints out of bounds", {
  spec <- matrix(rnorm(1000 * 10), nrow = 1000)
  endpoints_index <- c(1, 300, 500, 700, 1003)
  expect_error(rep_collapsed_by_index(spec, endpoints_index),
               "endpoints must be between")
})

test_that("one (column) spec provided", {
  spec <- rep(1:2, each = 10)
  endpoints_index <- c(11, 21)
  expected <- matrix(1:2, nrow = 1, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index), expected)
})

test_that("return correct collapsed measures w/o 1 or nfreq + 1; multi-clust", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- rep(1:2, each = 5)
  endpoints_index <- matrix(c(11, 21), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index, labels), expected)
})

test_that("return correct collapsed measures w/o 1 but w/ nfreq + 1; multi-clust", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- rep(1:2, each = 5)
  endpoints_index <- matrix(c(11, 31, 21, 31), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index, labels), expected)
})

test_that("return correct collapsed measures w/ 1 but w/o nfreq + 1; multi-clust", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- rep(1:2, each = 5)
  endpoints_index <- matrix(c(1, 11, 1, 21), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index, labels), expected)
})

test_that("return correct collapsed measures w/ 1 and nfreq + 1; multi-clust", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- rep(1:2, each = 5)
  endpoints_index <- matrix(c(1, 11, 31, 1, 21, 31), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_equal(rep_collapsed_by_index(spec, endpoints_index, labels), expected)
})

test_that("multiple endpoints but no labels", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  endpoints_index <- matrix(c(1, 11, 1, 21), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_warning(rep_collapsed_by_index(spec, endpoints_index), "no labels")
})

test_that("not enough sets of endpoints", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3)
  endpoints_index <- matrix(c(1, 11, 1, 21), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_error(rep_collapsed_by_index(spec, endpoints_index, labels),
               "not enough sets of frequency bands")
})

test_that("too many labels", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
  endpoints_index <- matrix(c(1, 11, 1, 21), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_error(rep_collapsed_by_index(spec, endpoints_index, labels),
               "too many labels provided")
})

test_that("not enough labels", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- c(1, 1, 1, 1, 1, 3, 3, 3)
  endpoints_index <- matrix(c(1, 11, 1, 21), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_error(rep_collapsed_by_index(spec, endpoints_index, labels),
               "not enough labels provided")
})

test_that("labels outside valid range", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- rep(c(1, 4), each = 5)
  endpoints_index <- matrix(c(1, 11, 1, 21), nrow = 2, byrow = T)
  expected <- matrix(c(
    rep(1:2, 5), rep(3:4, 5)
  ), nrow = 10, byrow = T)
  expect_error(rep_collapsed_by_index(spec, endpoints_index, labels),
               "some labels are outside the range 1:nclust")
})

test_that("return correct avg collapsed measures w/ 1 and nfreq + 1; multi-clust", {
  spec <- matrix(c(
    rep(rep(1:2, c(10, 20)), 5), rep(rep(3:4, c(20, 10)), 5)
  ), nrow = 30, ncol = 10)
  labels <- rep(1:2, each = 5)
  endpoints_index <- matrix(c(1, 11, 31, 1, 21, 31), nrow = 2, byrow = T)
  expected <- matrix(1:4, nrow = 2, byrow = T)
  expect_equal(avg_collapsed_by_index(spec, endpoints_index, labels), expected)
})
