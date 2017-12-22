context("exon selection")


test_that("filter by distance works", {
  expairs <-   list(
    tx1 =  list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100),
                                         width = 10),
                        exon_id = c("e1", "e2")
                        )), 
    tx2 =  list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100),
                                         width = 15), 
                        exon_id = c("e3", "e4"))),
    tx3 =  list(GenomicRanges::GRanges(seqnames = "a",
                        ranges = IRanges::IRanges(c(1, 100),
                                          c(15,200)), 
                        exon_id = c("e5", "e6"))
                ))
    
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
