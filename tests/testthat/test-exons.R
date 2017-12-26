context("exon selection")


test_that("filter by distance works", {
  expairs <- list(
    tx1 = list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100), width = 10),
                        exon_id = c("e1", "e2"))), 
    tx2 = list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100), width = 15), 
                        exon_id = c("e3", "e4"))),
    tx3 = list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100), c(15,200)), 
                        exon_id = c("e5", "e6")))
  )
  minDiff <- 0
  res <- expairs
  expect_equal(filterByDistance(expairs, minDiff), res)
  minDiff <- 5
  res[[2]] <- res[[3]][integer(0)]
  res[[3]] <- res[[3]][integer(0)]
  expect_equal(filterByDistance(expairs, minDiff), res)
  minDiff <- 6
  res[[1]] <- res[[1]][integer(0)]
  expect_equal(filterByDistance(expairs, minDiff), res)
})

test_that("asserion works", {
  expair <- GenomicRanges::GRanges(
    seqnames = "a",
    ranges = IRanges::IRanges(c(1, 100), width = 10),
    exon_id = c("e1", "e2")
  )
  expect_silent(assertColumns(expair, "exon_id"))
  expect_error (assertColumns(expair, "exon_ip"))
})

test_that("single exon tx are correctly processed", {
  singleExon <-  GenomicRanges::GRanges(
    seqnames = "a",
    ranges = IRanges::IRanges(1, width = 15),
    exon_id = c("single")
  )
  expect_equal(list(), splitToPairs(singleExon))
})

test_that("splitting works", {
  threeExons <- GenomicRanges::GRanges(
    seqnames = "a",
    ranges = IRanges::IRanges(c(1, 100, 200), width = 10),
    exon_id = c("e1", "e2", "e3")
  )
  expect_equal(splitToPairs(threeExons),
               list(c(threeExons[1], threeExons[2]),
                    c(threeExons[2], threeExons[3])))
})

test_that("exonSeqs input checks work", {
  expairs <- list(
    tx1 = list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100), width = 10),
                        exon_id = c("e1", "e2"))), 
    tx2 = list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100), width = 15), 
                        exon_id = c("e3", "e4"))),
    tx3 = list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100), c(15,200)), 
                        exon_id = c("e5", "e6")))
  )
  expect_error(addSeq(expairs, list()))
  dna <- mockDNAStringSet("a", 200)
  expect_error(addSeq(list(), dna))
  expect_error(addSeq(list(), iris))
  expect_silent(addSeq(expairs[[1]], dna))
  expairs[[1]][[1]]$exon_id <- NULL
  expect_error(addSeq(expairs[[1]], dna))
})

test_that("check for exonsByTx class works", {
  x <- 1
  expect_error(extractExonPairs(list(x)))
  class(x) <- "GRanges"
  expect_silent(extractExonPairs(list()))
})