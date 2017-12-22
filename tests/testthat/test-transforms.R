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