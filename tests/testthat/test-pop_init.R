test_that("correct dimensions", {
  spec <- matrix(1, nrow = 100, ncol = 10)
  popsize <- 100
  nclust <- 2; nbands <- 4
  pop <- pop_init(popsize, nclust, nbands, spec)
  for (i in 1:popsize) {
    expect_equal(dim(pop[[i]]), c(nclust, 2 * nbands - 1))
  }
})

test_that("equal endpoints across clusters", {
  spec <- matrix(1, nrow = 100, ncol = 10)
  popsize <- 100
  nclust <- 2; nbands <- 4
  pop <- pop_init(popsize, nclust, nbands, spec)
  for (i in 1:popsize) {
    endpoints <- pop[[i]][, 1:(nbands - 1), drop = FALSE]
    mean <- colMeans(endpoints)
    for (j in 1:nclust) {
      expect_equal(pop[[i]][j, 1:(nbands - 1)], mean)
    }
  }
})

test_that("expected output with nclust, nbands > 1", {
  spec <- matrix(1, nrow = 100, ncol = 10)
  popsize <- 100
  nclust <- 2; nbands <- 4
  pop <- pop_init(popsize, nclust, nbands, spec)
  for (i in 1:popsize) {
    expect_equal(pop[[i]][, nbands:(nbands + nbands - 1)],
                 matrix(1, nrow = nclust, ncol = nbands))
    expect_in(pop[[i]][1, 1:(nbands - 1)], 1:100)
  }
})

test_that("expected output with nclust = 1, nbands > 1", {
  spec <- matrix(1, nrow = 100, ncol = 10)
  popsize <- 100
  nclust <- 1; nbands <- 4
  pop <- pop_init(popsize, nclust, nbands, spec)
  for (i in 1:popsize) {
    expect_equal(pop[[i]][, nbands:(nbands + nbands - 1), drop = FALSE],
                 matrix(1, nrow = nclust, ncol = nbands))
    expect_in(pop[[i]][1, 1:(nbands - 1)], 1:100)
  }
})

test_that("expected output with large s", {
  spec <- matrix(1, nrow = 100, ncol = 100)
  popsize <- 100
  nclust <- 6; nbands <- 4
  pop <- pop_init(popsize, nclust, nbands, spec, s = 100)
  for (i in 1:popsize) {
    expect_equal(pop[[i]][, nbands:(nbands + nbands - 1), drop = FALSE],
                 matrix(1, nrow = nclust, ncol = nbands))
    expect_in(pop[[i]][1, 1:(nbands - 1)], 1:100)
  }
})
