context("settings")

test_that("checks for user input", {
  expect_error(p3Settings() %>% primerTm(100,50,60))
  expect_error(p3Settings() %>% primerTm(10,50,20))
  expect_silent(p3Settings() %>% primerTm(10,20,30))
  expect_silent(p3Settings() %>% primerTm(10,20,20))
  expect_silent(p3Settings() %>% primerTm(10,10,20))
  
  expect_error(p3Settings() %>% productSize(1:3))
  expect_error(p3Settings() %>% productSize(2:1))
  expect_silent(p3Settings() %>% productSize(1:2))
  
  expect_error(p3Settings() %>% primerSize(100,10,20))
  expect_error(p3Settings() %>% primerSize(10,40,30))
  expect_silent(p3Settings() %>% primerSize(10,20,30))
  expect_silent(p3Settings() %>% primerSize(10,20,20))
  expect_silent(p3Settings() %>% primerSize(10,10,20))
})