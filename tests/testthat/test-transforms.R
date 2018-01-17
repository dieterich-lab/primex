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

test_that("splitting works", {
  pair <- data.frame(
    start = c(1, 100),
    width = c(10, 10),
    direction = c("forward", "reverse")
  )
  expect_equivalent(splitPrimer(pair, 0), pair)
  expect_equivalent(splitPrimer(pair, 109), pair)
  res <-  data.frame(
    start = c(1, 6, 100),
    width = c(5, 5, 10),
    direction = c(rep("forward", 2), "reverse")
  )
  expect_equivalent(splitPrimer(pair, 5), res)
  res <-  data.frame(
    start = c(1, 2, 100),
    width = c(1, 9, 10),
    direction = c(rep("forward", 2), "reverse")
  )
  expect_equivalent(splitPrimer(pair, 1), res)
  res <-  data.frame(
    start = c(1, 100, 101),
    width = c(10, 1, 9),
    direction = c("forward", rep("reverse", 2))
  )
  expect_equivalent(splitPrimer(pair, 100), res)
})

test_that("pair to genome works + strand", {
  pair <- data.frame(
    start = c(1, 11, 1001 + 10),
    width = c(10, 5, 10),
    direction = c("forward", "reverse")[c(1,1,2)]
  )
  upExon <- GRanges(
    seqnames = 'a',
    strand = Rle('+', 1),
    ranges = IRanges(start = 1, end = 10)
  )
  downExon <- GRanges(
    seqnames = 'a',
    strand = Rle('+', 1),
    ranges = IRanges(start = 1001, end = 10000)
  )
  exons <- c(upExon, downExon)
  primerRanges <- GRanges(
    seqnames = "a",
    ranges = IRanges(start = c(1, 1001, 2001),
                     end = c(10, 1005, 2010)),
    strand = "+")
  mcols(primerRanges) <- pair["direction"]
  res <- pairToGenome(pair, exons)
  expect_equal(res, primerRanges)
})

test_that("pair to genome works - strand", {
  pair <- data.frame(
    start = c(1, 11, 1001 + 10),
    width = c(10, 5, 10),
    direction = c("forward", "reverse")[c(1,1,2)]
  )
  upExon <- GRanges(
    seqnames = 'a',
    strand = Rle('-', 1),
    ranges = IRanges(start = 2991, end = 3000)
  )
  downExon <- GRanges(
    seqnames = 'a',
    strand = Rle('-', 1),
    ranges = IRanges(start = 1, end = 2000)
  )
  exons <- c(upExon, downExon)
  primerRanges <- GRanges(
    seqnames = "a",
    ranges = IRanges(start = c(2991, 1996, 991),
                     end = c(3000, 2000, 1000)),
    strand = "-")
  mcols(primerRanges) <- pair["direction"]
  
  res <- pairToGenome(pair, exons)
  expect_equal(res, primerRanges)
})