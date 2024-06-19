library(testthat)
library(gnomesims)

test_that("power function works", {
  expect_equal(round(gnome_power(df = 1, ncp = 10), 2), 0.89)
})
