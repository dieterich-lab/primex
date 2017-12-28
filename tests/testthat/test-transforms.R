context("coordinate transforms")


test_that("toGenome works properly", {
  downExon <- GRanges(
    seqnames = 'a',
    strand = Rle('+', 1),
    ranges = IRanges(start = 1, end = 10)
  )
  upExon <- GRanges(
    seqnames = 'a',
    strand = Rle('+', 1),
    ranges = IRanges(start = 1001, end = 1010)
  )
  # inside to genome coords
  expect_equal(1004, primex:::.toGenome(4, c(upExon, downExon)))
  expect_equal(4, primex:::.toGenome(14, c(upExon, downExon)))
  # swap and minus strand
  strand(downExon) <- "-"
  strand(upExon)   <- "-"
  expect_equal(7, primex:::.toGenome(4, c(downExon, upExon)))
  expect_equal(1007, primex:::.toGenome(14, c(downExon, upExon)))
})

test_that("coordinate extracted correctly", {
  result <- data.frame(start = 10, width = 20)
  expect_equal(extractCoords("10,20"), result)
})

test_that("left and right treated correctly", {
  p <- list(PRIMER_LEFT = c("0,1", "9,10"),
            PRIMER_RIGHT = c("9,1", "99,10"))
  res <- list(
    data.frame(
      start = c(1, 10),
      width = c(1, 1),
      direction = c("forward", "reverse")
    ),
    data.frame(
      start = c(10, 91),
      width = c(10, 10),
      direction = c("forward", "reverse")
    )
  )
  expect_equivalent(primersToPairs(p), res)
})
