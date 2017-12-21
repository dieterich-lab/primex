
.toGenome <- function(x, strand, upExon, downExon) {
  vapply(x, function(y) {
    if (strand == '+') {
      if (y <= IRanges::width(upExon)) {
        res <- y + IRanges::start(upExon) - 1
      } else {
        res <- y - IRanges::width(upExon) + IRanges::start(downExon) - 1
      }
      res
    } else {
      if (y <= IRanges::width(upExon)) {
        res <- IRanges::end(upExon) - y + 1
      } else {
        res <-  IRanges::end(downExon) - (y - IRanges::width(upExon)) + 1
      }
      res
    }
  }, numeric(1))
}
