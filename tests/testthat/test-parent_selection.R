test_that("SUS gives expected ratios", {
  # probability of selection is proportional to index of the value
  # for example, item 10 in items should be selected about 10 times as often
  # as item 1, about 5 times as often as item 2, and about twice as often as
  # item 5.
  for (i in 1:10) {
    items <- letters[1:10]; x <- 1:10; cdf <- cumsum(x / sum(x))
    selected <- stochastic_universal_sampling(cdf, 1000000)
    ratios <- table(selected) / length(selected)
    exp_ratio <- seq(from = 10, to = 100, by = 10)
    expect_true(all(abs(ratios / ratios[10] * 100 - exp_ratio) < 0.01))
  }
})
