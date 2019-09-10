context("count_n_mb_events")

test_that("no event", {

  newick <- "((A:1, B:1):1);"
  tree <- ape::read.tree(text = newick)
  n_events <- count_n_mb_events(ape::branching.times(tree))
  expect_equal(0, n_events)
})

test_that("one event", {

  newick <- "((A:1, B:1):1, (C:1, D:1):1);"
  tree <- ape::read.tree(text = newick)
  n_events <- count_n_mb_events(ape::branching.times(tree))
  expect_equal(1, n_events)
})
