library(testthat)
library(gnomesims)

test_that("effect function works mz", {
  expect_equal(round(gnome_effect(a = sqrt(.4), c = sqrt(.3), e = sqrt(.3), g = 0, b = sqrt(.05))$mz, 2), 0.33)
})

test_that("effect function works dz", {
  expect_equal(round(gnome_effect(a = sqrt(.4), c = sqrt(.3), e = sqrt(.3), g = 0, b = sqrt(.05))$dz, 2), 0.19)
})
