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

primerCoords <- function(coords) {
  coords <- do.call(
    rbind,
    lapply(coords,
           function(x) {
             x <- strsplit(x, ",")[[1]]
             as.numeric(x)
           }))
  colnames(coords) <- c("start", "width")
  as.data.frame(coords)
}

primersToRanges <- function(res) {
  lefts <- do.call(IRanges::IRanges, primerCoords(res$PRIMER_LEFT))
  rights <- do.call(IRanges::IRanges, primerCoords(res$PRIMER_RIGHT))
  lapply(seq_along(lefts), function(i) {
    c(lefts[i], rights[i])
  })
}

pairToGRanges <- function(pair, seqs) {
  seqstrand <- as.character(strand(seqs))[1]
  res <- lapply(
    list(start = IRanges::start(pair), end = IRanges::end(pair)),
    .toGenome,
    upExon = seqs[1],
    downExon = seqs[2],
    strand = seqstrand
  )
  res <- do.call(rbind, Map(range, res$start, res$end))
  res <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(
      unique(GenomicRanges::seqnames(seqs[1])), lengths = 2),
    IRanges::IRanges(start = res[, 1], end = res[, 2]),
    strand = S4Vectors::Rle(unique(GenomicRanges::strand(seqs[1])),2)
  )
}