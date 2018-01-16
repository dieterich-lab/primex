context("plot and transformations")

test_that("unify difference works", {
  x <- c(10,20)
  y <- c(30,40)
  expect_equal(unifyDiff(x,y), list(c(1,2), c(3,4)))
  
  x <- c(10,20)
  y <- c(20,40)
  expect_equal(unifyDiff(x,y), list(c(1,2), c(2,3)))
  
  x <- c(10,20,20,21,20,19)
  y <- c(20,15,16,30,19,40)
  expect_equal(unifyDiff(x,y), 
               list(c(1,5,5,6,5,4), c(5,2,3,7,4,8)))
})

test_that("normalisation works", {
  x <- c(10,20,20,21,20,19)
  y <- c(20,15,16,30,19,40)
  result <-
    list(a = data.frame(start = c(1, 5, 5, 6, 5, 4), end = c(5, 2, 3, 7, 4, 8)))
  expect_equal(normaliseData(a = data.frame(start = x, end = y)), result)
               
})

test_that("plotting without errors", {
  
})