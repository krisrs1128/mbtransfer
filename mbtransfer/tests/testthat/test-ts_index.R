test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

data(sim_ts)
test_that("Dimension of appended interventions makes sense", {
  new_inter <- matrix(0, ncol = 5, nrow = 1)
  sim_ts2 <- replace_inter(sim_ts, new_inter)
  for (i in seq_along(sim_ts)) {
    expect_equal(
      ncol(interventions(sim_ts2[[i]])),
      ncol(interventions(sim_ts[[i]])) + 1
    )
  }
})

test_that("Number of samples returned makes sense", {
  expect_equal(length(sample_ts(ts, 10)), 10)
  expect_equal(length(sample_ts(ts, 1)), 1)
})

test_that("Patch length across samples makes sense", {
  for (j in c(4, 10, 1)) {
    ts_star <- sample_ts(ts, 10, j)
    for (i in seq_along(ts_star)) {
      expect_equal(ncol(ts_star[[i]]), j)
    }
  }
})